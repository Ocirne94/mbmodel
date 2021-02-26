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
boot_file_read    <- FALSE                 # Load .RData file with the input data, instead of loading input files.
boot_file_name    <- "boot_file_gries.RData"    # Name of the .RData input data file.


#### Include libraries ####
library(raster)
library(spatialEco)   # curvature()
library(scales)       # rescale()
library(Rcpp)         # avalanche function implemented in C++ for performance
library(gstat)        # IDW of snow probing data
library(Rfast)        # rowSort() of the stake cells indices
library(stats)        # uniroot()
library(sf)           # geom_sf(), to plot the glacier outline in the maps
library(metR)         # geom_text_contour()
library(ggplot2)
library(RStoolbox)    # For the surface type basemap under the SWE plots.

# library(profvis) # profiling.
# library(bitops)  # cksum(), for debugging.



#### Load function definitions ####
invisible(sapply(paste("functions/", list.files("functions/", pattern = "\\.R$"), sep=""), source))
source("func_set_params.R")
sourceCpp("functions/func_avalanche_gruber.cpp", cacheDir = "functions/") # Remove cacheDir option to force reload of the C++ code (useful after changing computer or editing the source file).


#### Load or set fixed run parameters ####
if (params_file_read) {
  load(params_file_name)
  } else {
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
  data_outlines              <-   func_load_outlines(run_params)
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
  grids_ice_albedo_fact      <-   func_compute_variable_ice_albedo(run_params, data_dhms)
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
invisible(gc())

# This will be set to TRUE after the first modeled year,
# to enable re-using of the modeled SWE as starting condition.
swe_prev_available <- FALSE

#### Main loop ####
for (year_id in 1:run_params$n_years) {
  # year_id <- 1 # TESTING!!

  # Select current year.
  year_cur <- run_params$years[year_id]
  year_cur_params <- func_load_year_params(run_params, year_cur)
  
  cat("\n\n\n\n============  STARTING NEW YEAR:", year_cur, " ============\n")
  
  # Select grids of the current year.
  # We may be using different ids for elevation and surface type
  # since elevation grids could be interpolated (unlike surface type).
  # The fixed avalanche grids use the same indices as the elevation ones.
  elevation_grid_id <- data_dhms$grid_year_id[year_id]
  surftype_grid_id  <- data_surftype$grid_year_id[year_id]
  outline_id        <- data_outlines$outline_year_id[year_id]
  
  # Extract pre-computed avalanche grids for this year.
  grids_avalanche_cur <- sapply(grids_avalanche, `[[`, elevation_grid_id)
  
  # Compute reduced-intensity base topographic distribution of solid precipitation.
  dist_topographic_values      <- getValues(grids_snowdist_topographic[[elevation_grid_id]])
  dist_topographic_values_mean <- mean(dist_topographic_values)
  dist_topographic_values_red  <- dist_topographic_values_mean + run_params$accum_snow_dist_red_fac * (dist_topographic_values - dist_topographic_values_mean)

  # Extract ice albedo factor grid for this year.
  grid_ice_albedo_fact_cur_values <- getValues(grids_ice_albedo_fact[[elevation_grid_id]])
  
  # Select mass balance measurements of the current year.
  massbal_annual_ids <- func_select_year_measurements(data_massbalance_annual, year_cur)
  nstakes_annual <- length(massbal_annual_ids)
  massbal_winter_ids <- func_select_year_measurements(data_massbalance_winter, year_cur)
  nstakes_winter <- length(massbal_winter_ids)
  massbal_annual_meas_cur <- data_massbalance_annual[massbal_annual_ids,]
  massbal_winter_meas_cur <- data_massbalance_winter[massbal_winter_ids,] # Empty if we have no winter stakes for the year.
  

  # Find (vectorized) the distance of each annual (and then winter) stake from the 4 surrounding cell centers.
  # We will use this later to extract the modeled series for each stake.
  # dx1 = x distance from the two cells to the left,
  # dy2 = y distance from the two cells above.
  dx1_annual <- (massbal_annual_meas_cur$x - (extent(data_dhms$elevation[[elevation_grid_id]])[1] - (run_params$grid_cell_size / 2))) %% run_params$grid_cell_size
  dx2_annual <- run_params$grid_cell_size - dx1_annual
  dy1_annual <- (massbal_annual_meas_cur$y - (extent(data_dhms$elevation[[elevation_grid_id]])[3] - (run_params$grid_cell_size / 2))) %% run_params$grid_cell_size
  dy2_annual <- run_params$grid_cell_size - dy1_annual
  
  dx1_winter <- (massbal_winter_meas_cur$x - (extent(data_dhms$elevation[[elevation_grid_id]])[1] - (run_params$grid_cell_size / 2))) %% run_params$grid_cell_size
  dx2_winter <- run_params$grid_cell_size - dx1_winter
  dy1_winter <- (massbal_winter_meas_cur$y - (extent(data_dhms$elevation[[elevation_grid_id]])[3] - (run_params$grid_cell_size / 2))) %% run_params$grid_cell_size
  dy2_winter <- run_params$grid_cell_size - dy1_winter
  
  
  
  # Should we make a winter run to optimize the precipitation correction?
  # Only if we have some measurements of winter snow cover, else we can't.
  process_winter <- (nstakes_winter > 0)
  
  if (process_winter) {
    dist_probes_idw           <- func_snow_probes_idw(run_params, massbal_winter_meas_cur, data_dhms)
    dist_probes_idw           <- clamp(dist_probes_idw, lower = 0, upper = Inf)
    dist_probes_idw_norm      <- dist_probes_idw / mean(dist_probes_idw[data_dems$glacier_cell_ids[[elevation_grid_id]]])
  } else {
    # No winter probes to work with, so uniform distribution for the probes component.
    dist_probes_idw_norm      <- setValues(data_dhms$elevation[[1]], 1.0)
  }
  dist_probes_norm_values     <- getValues(dist_probes_idw_norm) # For the accumulation model.
  dist_probes_norm_mean       <- mean(dist_probes_norm_values, na.rm = T)
  dist_probes_norm_values_red <- dist_probes_norm_mean + run_params$accum_probes_red_fac * (dist_probes_norm_values - dist_probes_norm_mean)
  
  
  #### .  MODELING PERIOD BOUNDS ####
  # Vector of four Date objects: start and end of the annual
  # modeling period, and same for the winter modeling period.
  model_time_bounds   <- func_compute_modeling_periods(run_params,
                                                       massbal_annual_meas_cur,
                                                       massbal_winter_meas_cur,
                                                       year_cur,
                                                       year_cur_params)
  
  
  #### .  INITIAL SNOW COVER ####
  # The initial snow cover can be either (1) a single, constant estimated map
  # (from topography, avalanches and user-defined snow line elevation)
  # or (2) the result of the previous year of modeling at the starting date
  # of the simulation.
  # In case (2), we distinguish between initial snow cover for winter
  # and for annual modeling, since the two can have different
  # dates, depending on the dates of the annual and winter stakes.
  # Case (2) is obviously not applicable to the first year of modeling
  # (there is no previous result for it).
  if (run_params$initial_snow_dist_from_model && swe_prev_available) {
    
    # NOTE: here weather_series_annual_cur and mod_output_annual_cur are
    # still the weather series and modeled series of the PREVIOUS year!
    swe_prev_annual_day_id <- which.min(abs(weather_series_annual_cur$timestamp - model_time_bounds[1]))
    snowdist_init_annual <- setValues(data_dhms$elevation[[elevation_grid_id]], mod_output_annual_cur$vec_swe_all[(swe_prev_annual_day_id - 1) * run_params$grid_ncells + 1:run_params$grid_ncells])
    
    if (process_winter) {
      swe_prev_winter_day_id <- which.min(abs(weather_series_annual_cur$timestamp - model_time_bounds[3]))
      snowdist_init_winter <- setValues(data_dhms$elevation[[elevation_grid_id]], mod_output_annual_cur$vec_swe_all[(swe_prev_winter_day_id - 1) * run_params$grid_ncells + 1:run_params$grid_ncells])
    }
    
  # Here instead estimate the initial snow cover from snow line elevation,
  # topography, avalanches and snow probes if available.
  } else {
    snowdist_init_annual <- func_compute_initial_snow_cover(run_params,
                                                            data_dhms,
                                                            data_dems,
                                                            grids_snowdist_topographic,
                                                            grids_avalanche_cur,
                                                            dist_probes_idw_norm,
                                                            elevation_grid_id,
                                                            massbal_winter_meas_cur)
    snowdist_init_winter <- snowdist_init_annual
  }
  
  
  #### .  WINTER MASS BALANCE, to optimize the precipitation correction ####
  # We set this here so that there is no correction
  # if we don't do the winter optimization.
  corr_fact_winter      <- 0
  # We set this to NULL to have it defined (for the
  # extraction functions) in case we don't do winter processing.
  mod_output_winter_cur <- NULL
  if (process_winter)  {
    
    winter_stakes_cells <- rowSort(fourCellsFromXY(data_dhms$elevation[[elevation_grid_id]], as.matrix(massbal_winter_meas_cur[,4:5])))
    
    model_winter_bounds <- model_time_bounds[3:4]
    
    # Select weather series period.
    # If the weather series time specification is
    # wrong this step is where it all falls apart.
    weather_series_winter_cur <- data_weather[which(data_weather$timestamp == model_winter_bounds[1]):(which(data_weather$timestamp == model_winter_bounds[2])),]
    model_winter_days_n <- nrow(weather_series_winter_cur)
    
    # This leaves the result of the last optimization
    # iteration in mod_output_annual_cur.
    # The NA is for the corr_fact_winter (which we are
    # determining here, so we don't use a previous value:
    # it is ignored).
    optim_corr_winter <- func_optimize_mb("winter", NA,
                                          run_params, year_cur_params, elevation_grid_id, surftype_grid_id,
                                          data_dhms, data_dems, data_surftype, snowdist_init_winter, data_radiation,
                                          weather_series_winter_cur, dist_topographic_values_red,
                                          dist_probes_norm_values_red, grids_avalanche_cur,
                                          grid_ice_albedo_fact_cur_values,
                                          dx1_winter, dx2_winter, dy1_winter, dy2_winter,
                                          nstakes_winter, model_winter_days_n, massbal_winter_meas_cur,
                                          winter_stakes_cells)
    
    # Save the correction factor to re-use it
    # during the annual optimization.
    # We divide by the original prec_corr
    # since the corr_fact is relative
    # (it gets multiplied again during
    # optimization, inside func_optim_worker()).
    corr_fact_winter <- optim_corr_winter$prec_corr / year_cur_params$prec_corr
    
    # Free some memory after processing.
    invisible(gc())
  }
    
  
  #### .  ANNUAL MASS BALANCE ####
  # Here do the annual processing.
  # Find grid cells corresponding to the annual stakes.
  # We sort them to enable vectorized bilinear filtering.
  annual_stakes_cells <- rowSort(fourCellsFromXY(data_dhms$elevation[[elevation_grid_id]], as.matrix(massbal_annual_meas_cur[,4:5])))
  
  # Time specification of the annual run.
  model_annual_bounds <- model_time_bounds[1:2]
  
  # Select weather series period.
  # If the weather series time specification is
  # wrong this step is where it all falls apart.
  weather_series_annual_cur <- data_weather[which(data_weather$timestamp == model_annual_bounds[1]):(which(data_weather$timestamp == model_annual_bounds[2])),]
  model_annual_days_n <- nrow(weather_series_annual_cur)
  
  
  # This leaves the result of the last optimization
  # iteration in mod_output_annual_cur.
  optim_corr_annual <- func_optimize_mb("annual", corr_fact_winter,
                                        run_params, year_cur_params, elevation_grid_id, surftype_grid_id,
                                        data_dhms, data_dems, data_surftype, snowdist_init_annual, data_radiation,
                                        weather_series_annual_cur, dist_topographic_values_red,
                                        dist_probes_norm_values_red, grids_avalanche_cur,
                                        grid_ice_albedo_fact_cur_values,
                                        dx1_annual, dx2_annual, dy1_annual, dy2_annual,
                                        nstakes_annual, model_annual_days_n, massbal_annual_meas_cur,
                                        annual_stakes_cells)
  # Free some memory after processing.
  invisible(gc())
  
  # After an annual model run we have SWE information
  # suitable for use as starting condition of the next
  # year, if the user decides to use it.
  swe_prev_available <- TRUE
  
  
  #### . EXTRACT CUMULATIVE MASS BALANCE AT DATES OF INTEREST ####
  # We extract three maps of cumulative annual mass balances:
  # (1) "hydro":       hydrological year (1 October <Year-1> - 30 September <Year>)
  # (2) "meas_period": measurement period, defined as (earliest annnual stake start - latest annual stake end)
  # (3) "fixed":       user-defined fixed period.
  massbal_annual_maps_data <- func_extract_massbal_maps_annual(run_params,
                                                               year_cur_params,
                                                               weather_series_annual_cur,
                                                               mod_output_annual_cur,
                                                               data_dhms)
  massbal_annual_maps <- massbal_annual_maps_data$massbal_maps
  massbal_annual_meas_period <- massbal_annual_maps_data$meas_period
  
  # We also extract two winter mass balances:
  # (1) "fixed":       user-defined fixed period.
  # (2) "meas_period": measurement period, defined as (earliest winter stake start - latest winter stake end).
  # If process_winter is FALSE, the list contains only (1).
  massbal_winter_maps_data <- func_extract_massbal_maps_winter(run_params,
                                                               year_cur_params,
                                                               weather_series_annual_cur,
                                                               mod_output_annual_cur,
                                                               process_winter,
                                                               mod_output_winter_cur,
                                                               data_dhms)
  massbal_winter_maps <- massbal_winter_maps_data$massbal_maps
  if (process_winter) {
    massbal_winter_meas_period <- massbal_winter_maps_data$meas_period
  } else {
    massbal_winter_meas_period <- NA
  }
  
  
  #### . CORRECT ANNUAL MASS BALANCE IN ELEVATION BANDS ####
  massbal_annual_maps$meas_period_corr <- func_correct_massbal_elebands(year_cur_params,
                                                                        data_dems,
                                                                        massbal_annual_meas_cur,
                                                                        mod_output_annual_cur,
                                                                        massbal_annual_maps)
  
  massbal_annual_values <- sapply(massbal_annual_maps, cellStats, stat = "mean", na.rm = TRUE)
  massbal_winter_values <- sapply(massbal_winter_maps, cellStats, stat = "mean", na.rm = TRUE)
  
  
  
  #### . PLOTS OF THE EXTRACTED MASS BALANCE MAPS ####
  func_plot_year_mb_maps(run_params,
                         year_cur,
                         data_dems,
                         data_outlines,
                         elevation_grid_id,
                         outline_id,
                         massbal_annual_maps,
                         massbal_winter_maps,
                         massbal_annual_values,
                         massbal_winter_values,
                         massbal_annual_meas_period,
                         massbal_winter_meas_period,
                         process_winter)
  
  
    
  
  
  
  
  
  #### .  DAILY PLOTS (SLOW!) ####
  # Compute surface type image, to use as background for the daily SWE plots.
  # surf_r    <- subs(data_surftype$grids[[surftype_grid_id]], data.frame(from = c(0, 1, 4, 5), to = c(100, 170, 0, 70)))
  # surf_g    <- subs(data_surftype$grids[[surftype_grid_id]], data.frame(from = c(0, 1, 4, 5), to = c(150, 213, 0, 20)))
  # surf_b    <- subs(data_surftype$grids[[surftype_grid_id]], data.frame(from = c(0, 1, 4, 5), to = c(200, 255, 0, 20)))
  # surf_base <- ggRGB(stack(surf_r, surf_g, surf_b), r = 1, g = 2, b = 3, ggLayer = TRUE)
  # func_plot_daily_maps(run_params,
  #                      year_cur,
  #                      weather_series_annual_cur,
  #                      data_dems,
  #                      data_outlines,
  #                      mod_output_annual_cur,
  #                      surf_base,
  #                      elevation_grid_id,
  #                      outline_id)
  
  
}
