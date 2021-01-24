###################################################################################################
# Author:         Enrico Mattea (@unifr.ch), inspired by the IDL version by Matthias Huss.        #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the main loop and instructions.                              #
###################################################################################################

#### Load from files or reboot file ####
params_file_write <- FALSE                # Save .RData file with run parameters, for faster reload. The file is saved AFTER data loading since some parameters depend on the loaded grids (e.g. cell size).
params_file_read  <- FALSE                # Load .RData file with run parameters, instead of setting new run parametrs.
params_file_name  <- "params_file.RData"  # Name of the .RData run parameters file.

boot_file_write   <- FALSE                # Save .RData file with the input data, for faster reload.
boot_file_read    <- TRUE                 # Load .RData file with the input data, instead of loading input files.
boot_file_name    <- "boot_file.RData"    # Name of the .RData input data file.


#### Include libraries ####
library(raster)
library(spatialEco)   # curvature()
library(scales)       # rescale()
library(Rcpp)         # avalanche function implemented in C++ for performance
library(gstat)        # IDW of snow probing data
library(Rfast)        # rowSort() of the stake cells indices

# Experimental plots.
library(ggplot2)
library(RStoolbox)

library(profvis)


# library(bitops)  # cksum(), for debugging.



#### Load function definitions ####
invisible(sapply(paste("functions/", list.files("functions/", pattern = "\\.R$"), sep=""), source))
sourceCpp("functions/func_avalanche_gruber.cpp", cacheDir = "functions/") # Remove cacheDir option to force reload of the C++ code (useful after changing computer or editing the source file).


#### Load or set run parameters ####
if (params_file_read) {
  load(params_file_name) } else {
  run_params <- func_set_params()
}


#### Load input data from sources or reboot file ####
if (boot_file_read) {
  load(boot_file_name)
} else {
  data_weather               <-   func_load_weather(run_params)
  data_dems                  <-   func_load_elevation_grids(run_params, "dem")
  data_dhms                  <-   func_load_elevation_grids(run_params, "dhm")
  data_surftype              <-   func_load_surftype_grids(run_params)
  data_radiation             <-   func_load_radiation_grids(run_params)
  data_massbalance_annual    <-   func_load_massbalance_measurements(run_params, "annual")
  data_massbalance_winter    <-   func_load_massbalance_measurements(run_params, "winter")
}

# Assign global grid parameters to run_params.
# We do it explicitly here since we need to have
# loaded the elevation grids to get their parameters.
# These overwrite the loaded run parameters, to avoid
# the possibility of dangerous inconsistencies with the loaded data.
run_params$grid_nrow       <- nrow(data_dhms$elevation[[1]])
run_params$grid_ncol       <- ncol(data_dhms$elevation[[1]])
run_params$grid_cell_size  <- xres(data_dhms$elevation[[1]])
run_params$grid_ncells     <- run_params$grid_nrow * run_params$grid_ncol

# Compute fixed grids to speed up processing, but
# not if these were loaded from the boot file.
if (!boot_file_read) {
  grids_avalanche            <-   func_compute_avalanche_fixed_grids(run_params, data_dhms)
  grids_snowdist_topographic <-   func_compute_snowdist_topographic(run_params, data_dhms, data_dems)
}

# Save reboot file.
if (boot_file_write) {
  save(list = c(apropos("^data_"), apropos("^grids_")), file = boot_file_name)
}
# Save parameters file.
if (params_file_write) {
  save(list = "run_params", file = params_file_name)
}

# Cleanup memory (temporary variables during loading!)
gc()


#### Main loop ####
# for (year_id in 1:run_params$n_years) {
  year_id <- 1 # TESTING!!

  # Select current year.
  year_cur <- run_params$years[year_id]
  year_cur_params <- func_load_year_params(run_params, year_cur)
  year_cur_params$rad_fact_firn <- (year_cur_params$rad_fact_ice + year_cur_params$rad_fact_snow) / 2 # As per IDL implementation.
  
  # Select grids of the current year.
  grid_id <- data_dhms$grid_year_id[year_id]
  
  grids_avalanche_cur <- sapply(grids_avalanche, `[[`, grid_id)
  
  
  #------------------------- Compute image for background of daily SWE plots -------------------------#
  surf_r <- subs(data_surftype[[grid_id]], data.frame(from = c(0, 1, 4, 5), to = c(100, 170, 0, 70)))
  surf_g <- subs(data_surftype[[grid_id]], data.frame(from = c(0, 1, 4, 5), to = c(150, 213, 0, 20)))
  surf_b <- subs(data_surftype[[grid_id]], data.frame(from = c(0, 1, 4, 5), to = c(200, 255, 0, 20)))
  surf_base <- ggRGB(stack(surf_r, surf_g, surf_b), r = 1, g = 2, b = 3, ggLayer = TRUE)
  #---------------------------------------------------------------------------------------------------#
  
  
  # Select mass balance measurements of the current year.
  massbal_annual_ids <- func_select_year_measurements(data_massbalance_annual, year_cur)
  nstakes_annual <- length(massbal_annual_ids)
  massbal_winter_ids <- func_select_year_measurements(data_massbalance_winter, year_cur)
  nstakes_winter <- length(massbal_winter_ids)
  massbal_annual_cur <- data_massbalance_annual[massbal_annual_ids,]
  massbal_winter_cur <- data_massbalance_winter[massbal_winter_ids,] # Empty if we have no winter stakes for the year.
  
  # Should we make a winter run to optimize the precipitation correction?
  # Only if we have some measurements of winter snow cover, else we can't.
  process_winter <- (nstakes_winter > 0)
  
  if (process_winter) {
    dist_probes_idw           <- func_snow_probes_idw(run_params, massbal_winter_cur, data_dhms)
    dist_probes_idw           <- clamp(dist_probes_idw, lower = 0, upper = Inf)
    dist_probes_idw_norm      <- dist_probes_idw / mean(dist_probes_idw[data_dems$glacier_cell_ids[[grid_id]]])
  } else {
    # No winter probes to work with, so uniform distribution for the probes component.
    dist_probes_idw_norm      <- setValues(data_dhms$elevation[[1]], 1.0)
  }
  dist_probes_norm_values     <- getValues(dist_probes_idw_norm) # For the accumulation model.
  dist_probes_norm_mean       <- mean(dist_probes_norm_values)
  dist_probes_norm_values_red <- dist_probes_norm_mean + run_params$accum_probes_red_fac * (dist_probes_norm_values - dist_probes_norm_mean)
  
  
  
  #### .  INITIAL SNOW COVER ####
  snowdist_init <- func_compute_initial_snow_cover(run_params,
                                                   data_dhms,
                                                   data_dems,
                                                   grids_snowdist_topographic,
                                                   grids_avalanche_cur,
                                                   dist_probes_idw_norm,
                                                   grid_id,
                                                   massbal_winter_cur)
  
  
  #### .  MODELING PERIOD BOUNDS ####
  # Four POSIXct time objects: start and end of the annual
  # modeling period, and same for the winter modeling period.
  model_time_bounds   <- func_compute_modeling_periods(run_params,
                                                       massbal_annual_cur,
                                                       massbal_winter_cur,
                                                       year_cur)
  
  
  #### .  WINTER MASS BALANCE ####
  if (process_winter)  {
    model_winter_bounds <- model_time_bounds[3:4]
    # TODO:
      # find grid cells of winter stakes
      # run massbal_model (we have to run the whole model - not only on the stake grid cells - due to avalanches: we want the whole swe maps!)
  }
    
  
  #### .  ANNUAL MASS BALANCE ####
  # Here instead do the annual processing.
  # Find grid cells corresponding to the annual stakes.
  # We sort them to enable vectorized bilinear filtering.
  annual_stakes_cells <- rowSort(fourCellsFromXY(data_dhms$elevation[[grid_id]], as.matrix(massbal_annual_cur[,4:5])))
  
  # Time specification of the annual run.
  model_annual_bounds <- model_time_bounds[1:2]
  
  # Select weather series period.
  # If the weather series time specification is wrong
  # this step is where it all falls apart.
  weather_series_cur <- data_weather[which(data_weather$timestamp == model_time_bounds[1]):(which(data_weather$timestamp == model_time_bounds[2]) - 1),]
  model_days_n <- length(weather_series_cur[,1])
  
  dist_topographic_values      <- getValues(grids_snowdist_topographic[[grid_id]])
  dist_topographic_values_mean <- mean(dist_topographic_values)
  dist_topographic_values_red  <- dist_topographic_values_mean + run_params$accum_snow_dist_red_fac * (dist_topographic_values - dist_topographic_values_mean)
  
  
  #### . .  RUN MASS BALANCE MODEL ####
  mb_model_output <- func_massbal_model(run_params,
                                        year_cur_params,
                                        getValues(data_dhms$elevation[[grid_id]]),
                                        data_dems$glacier_cell_ids[[grid_id]],
                                        getValues(data_surftype[[grid_id]]),
                                        getValues(snowdist_init),
                                        data_radiation,
                                        weather_series_cur,
                                        dist_topographic_values_red,
                                        dist_probes_norm_values_red,
                                        grids_avalanche_cur)

  
  # Find indices of the days corresponding to the stake measurements.
  # We match w.r.t. weather_series_cur whose index is off by ~0.5 with the
  # mass balance (mb_model_out$gl_massbal_cumul[1] is the initial condition
  # (i.e. 0.0) at 00:00 of the first day, then the index of the weather series
  # corresponds to the full following 24 hours, then gl_massbal_cumul[2] is the
  # cumulative mass balance by the end of that same day.
  # So it would be equally correct to also shift all the day indices by one (little to no change).
  annual_stakes_start_ids <- pmatch(massbal_annual_cur$start_date,
                                    weather_series_cur$timestamp,
                                    duplicates.ok = TRUE)
  annual_stakes_end_ids   <- pmatch(massbal_annual_cur$end_date,
                                    weather_series_cur$timestamp,
                                    duplicates.ok = TRUE)

  
  
  
  
  #### Extract the whole modeled series for all stakes. ####
  
  stakes_series_mod_all <- matrix(NA, nrow = model_days_n + 1, ncol = nstakes_annual) # One row per day, one column per stake
  # We do a manual bilinear interpolation of the four
  # cells surrounding each stake, with a (fast) loop on stakes.
  # We have verified that it corresponds exactly to raster::extract(..., method = "bilinear"), but much faster.
  # First find (vectorized) the distance of each stake from the 4 surrounding cell centers.
  # dx1 = x distance from the two cells to the left,
  # dy2 = y distance from the two cells above.
  dx1 <- (massbal_annual_cur$x - (extent(data_dhms$elevation[[grid_id]])[1] - (run_params$grid_cell_size / 2))) %% run_params$grid_cell_size
  dx2 <- run_params$grid_cell_size - dx1
  dy1 <- (massbal_annual_cur$y - (extent(data_dhms$elevation[[grid_id]])[3] - (run_params$grid_cell_size / 2))) %% run_params$grid_cell_size
  dy2 <- run_params$grid_cell_size - dy1
  
  # Now extract the values.
  for (stake_id in 1:nstakes_annual) {
    
    # Cells are ordered like this:
    # 1 2
    # 3 4
    # with the stake somewhere in the middle.
    # Repeated cells (i.e. if the stake lies at
    # the same x and/or y as a cell) work fine.
    cell_series1 <- mb_model_output$vec_massbal_cumul[annual_stakes_cells[stake_id, 1] + seq(0,length(mb_model_output$vec_massbal_cumul)-1,run_params$grid_ncells)]
    cell_series2 <- mb_model_output$vec_massbal_cumul[annual_stakes_cells[stake_id, 2] + seq(0,length(mb_model_output$vec_massbal_cumul)-1,run_params$grid_ncells)]
    cell_series3 <- mb_model_output$vec_massbal_cumul[annual_stakes_cells[stake_id, 3] + seq(0,length(mb_model_output$vec_massbal_cumul)-1,run_params$grid_ncells)]
    cell_series4 <- mb_model_output$vec_massbal_cumul[annual_stakes_cells[stake_id, 4] + seq(0,length(mb_model_output$vec_massbal_cumul)-1,run_params$grid_ncells)]
    
    stakes_series_mod_all[, stake_id] <- (cell_series1 * dx2[stake_id] * dy1[stake_id] +
                                          cell_series2 * dx1[stake_id] * dy1[stake_id] +
                                          cell_series3 * dx2[stake_id] * dy2[stake_id] +
                                          cell_series4 * dx1[stake_id] * dy2[stake_id]) / (run_params$grid_cell_size^2)
    
  }
    
  
  # Find start date for stakes with NA (i.e. mass balance minimum of previous year):
    # Find the minimum of the modeled series between the beginning and February 28 (or customizable parameter!)
    # Set the found minimum as the stake period start.
  annual_stakes_start_ids_corr <- annual_stakes_start_ids  # We leave the original set unaltered, it will serve during optimization.
  stakes_start_unknown_ids <- which(is.na(annual_stakes_start_ids))
  
  # Latest possible day for the search of the stake start at the mass balance minimum.
  stakes_start_latest_id <- which(format(weather_series_cur$timestamp, "%m/%d") == run_params$stakes_unknown_latest_start)
  
  for (stake_cur_id in stakes_start_unknown_ids) {
    # cat("Finding start date for stake", stake_cur_id, "...\n")
    annual_stakes_start_ids_corr[stake_cur_id] <- which.min(stakes_series_mod_all[1:stakes_start_latest_id, stake_cur_id])
  }
  
  # Cumulative mass balance of each stake
  # over the measurement period (numeric vector).
  stakes_mb_mod  <- as.numeric(stakes_series_mod_all)[((1:nstakes_annual)-1)*(model_days_n+1) + annual_stakes_end_ids] -
                    as.numeric(stakes_series_mod_all)[((1:nstakes_annual)-1)*(model_days_n+1) + annual_stakes_start_ids_corr]
  # Corresponding measurement.
  stakes_mb_meas <- massbal_annual_cur$dh_cm * massbal_annual_cur$density * 10 # 10: cm w.e. to mm w.e.
  
  stakes_bias <- stakes_mb_mod - stakes_mb_meas

  # Compute mean BIAS and RMS: mean(bias_vec) and sqrt(mean(bias_vec^2))
  # Then we're ready to optimize the model parameters.
  
  
  
# }
