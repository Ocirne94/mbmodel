###################################################################################################
# Author:         Enrico Mattea (@unifr.ch), inspired by the IDL version by Matthias Huss.        #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the main loop and instructions.                              #
###################################################################################################

Sys.setlocale(category = "LC_TIME", locale = "en_US.UTF-8") # Set English language for dates (in the plots).

#### Load from files or reboot file ####
params_file_write <- FALSE                # Save .RData file with run parameters, for faster reload. The file is saved AFTER data loading since some parameters depend on the loaded grids (e.g. cell size).
params_file_read  <- FALSE                # Load .RData file with run parameters, instead of setting new run parametrs.
params_file_name  <- "params_file.RData"  # Name of the .RData run parameters file.

boot_file_write   <- TRUE                # Save .RData file with the input data, for faster reload.
boot_file_read    <- FALSE                 # Load .RData file with the input data, instead of loading input files.
boot_file_name    <- "boot_file_barkrak.RData"    # Name of the .RData input data file.


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
library(ggplot2)      # Base plotting library
library(ggpubr)       # Additional plotting functions
library(grid)         # Additional plotting functions
library(cowplot)      # Additional plotting functions
library(ggplotify)    # Additional plotting functions
library(ggpattern)    # Additional plotting functions (install with remotes::install_github("coolbutuseless/ggpattern"))
library(reshape2)     # melt()
library(RStoolbox)    # For the surface type basemap under the SWE plots


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
  func_load_data_all(run_params) # This calls all the data loading routines.
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

# Create output directory.
dir.create(run_params$output_dirname, recursive = TRUE)

# Cleanup memory (temporary variables during loading!)
invisible(gc())

# This will be set to TRUE after the first modeled year,
# to enable re-using of the modeled SWE as starting condition.
swe_prev_available <- FALSE

# In this data frame we put all the annual results
# for the overview document: glacier-wide mass balance
# (annual and winter, all versions), ELA, AAR, RMSE,
# optimized parameters, cumulative mass balance.
df_overview <- data.frame(year                = run_params$years,
                          mb_annual_meas_corr = NA,
                          mb_annual_meas      = NA,
                          mb_annual_hydro     = NA,
                          mb_annual_fixed     = NA,
                          mb_winter_meas      = NA, # This stays NA unless winter measurements are available.
                          mb_winter_fixed     = NA,
                          ela                 = NA,
                          aar                 = NA,
                          rmse                = NA,
                          melt_factor         = NA,
                          rad_fact_ice        = NA,
                          rad_fact_snow       = NA,
                          prec_corr           = NA,
                          mb_cumul            = NA)

# Here we will put just the final mass balance for each
# year, to produce the overview_areaplot multi-page PDF file.
overview_areaplots <- list()

#### Main loop ####
for (year_id in 1:run_params$n_years) {
  # year_id <- 1 # TESTING

  # Select current year and the corresponding parameters.
  year_cur <- run_params$years[year_id]
  year_cur_params <- func_load_year_params(run_params, year_cur)
  
  cat("\n\n\n\n============  STARTING NEW YEAR:", year_cur, " ============\n")
  
  # Select grids of the current year from the list of available grids.
  # We could be using different ids for elevation and surface type
  # because the elevation grids can also be interpolated annually (unlike surface type).
  # So we have different _id variables.
  # The fixed avalanche grids use the same indices as the elevation ones.
  elevation_grid_id <- data_dhms$grid_year_id[year_id]
  surftype_grid_id  <- data_surftype$grid_year_id[year_id]
  outline_id        <- data_outlines$outline_year_id[year_id]
  
  # Extract avalanche grids for this year
  # (pre-computed before the start of the loop).
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
    weather_series_winter_cur <- data_weather[which(data_weather$timestamp == model_winter_bounds[1]):(which(data_weather$timestamp == model_winter_bounds[2])),]
    model_winter_days_n <- nrow(weather_series_winter_cur)
    
    # This leaves the result of the last optimization
    # iteration in mod_output_annual_cur.
    # The NA is for the optimized corr_fact_winter (which we are
    # determining here, so we don't use a previous value: it is ignored).
    optim_corr_winter <- func_optimize_mb("winter", NA,
                                          run_params, year_cur_params, elevation_grid_id, surftype_grid_id,
                                          data_dhms, data_dems, data_surftype, snowdist_init_winter, data_radiation,
                                          weather_series_winter_cur, dist_topographic_values_red,
                                          dist_probes_norm_values_red, grids_avalanche_cur,
                                          grid_ice_albedo_fact_cur_values,
                                          dx1_winter, dx2_winter, dy1_winter, dy2_winter,
                                          nstakes_winter, model_winter_days_n, massbal_winter_meas_cur,
                                          winter_stakes_cells)
    
    # Save the correction factor to re-use it during the annual optimization.
    # We divide by the original prec_corr since the corr_fact is relative
    # (it gets multiplied again during optimization, inside func_optim_worker()).
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
  weather_series_annual_cur <- data_weather[which(data_weather$timestamp == model_annual_bounds[1]):(which(data_weather$timestamp == model_annual_bounds[2])),]
  model_annual_days_n <- nrow(weather_series_annual_cur)
  
  
  # This returns the optimized parameters
  # and also leaves the result of the
  # last optimization iteration in
  # a variable called mod_output_annual_cur.
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
  
  # Extract the grid values since we will use them several times below.
  mb_meas_period_corr_values <- getValues(massbal_annual_maps$meas_period_corr)
  
  
  #### . COMPUTE ELA and AAR ####
  ela_aar <- func_compute_ela_aar(run_params, mb_meas_period_corr_values, data_dems)
  
  
  #### . SAVE OVERVIEW VALUES for the year ####
  df_overview$mb_annual_meas_corr[year_id] <- massbal_annual_values[["meas_period_corr"]] / 1e3
  df_overview$mb_annual_meas[year_id]      <- massbal_annual_values[["meas_period"]] / 1e3
  df_overview$mb_annual_hydro[year_id]     <- massbal_annual_values[["hydro"]] / 1e3
  df_overview$mb_annual_fixed[year_id]     <- massbal_annual_values[["fixed"]] / 1e3
  if (process_winter) {df_overview$mb_winter_meas[year_id] <- massbal_winter_values[["meas"]] / 1e3}
  df_overview$mb_winter_fixed[year_id]     <- massbal_winter_values [["fixed"]] / 1e3
  df_overview$ela[year_id]                 <- ela_aar[["ela"]]
  df_overview$aar[year_id]                 <- ela_aar[["aar"]] * 100
  df_overview$rmse[year_id]                <- mod_output_annual_cur$global_rms / 1e3
  df_overview$melt_factor[year_id]         <- year_cur_params$melt_factor + optim_corr_annual$melt_factor
  df_overview$rad_fact_ice[year_id]        <- year_cur_params$rad_fact_ice + optim_corr_annual$rad_fact_ice
  df_overview$rad_fact_snow[year_id]       <- year_cur_params$rad_fact_snow + optim_corr_annual$rad_fact_ice * year_cur_params$rad_fact_ratio_snow_ice
  df_overview$prec_corr[year_id]           <- year_cur_params$prec_corr + optim_corr_annual$prec_corr
  
  
  #### . PLOT THE MASS BALANCE MAPS ####
  # This returns a list with the (5 or 6, depending on whether we have winter measurements)
  # mass balance maps for the current year.
  # Then we will append to this list also the
  # other plots of the year (time series,
  # vertical distributions and so on).
  plots_year <- func_plot_year_mb_maps(run_params,
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
  
  
  #### . PLOT THE DAILY TIME SERIES OF GLACIER-WIDE MASS BALANCE ####
  plots_mb_cumul <- as.ggplot(func_plot_massbal_cumul(run_params,
                                                      process_winter,
                                                      massbal_annual_meas_cur,
                                                      massbal_winter_meas_cur,
                                                      mod_output_annual_cur,
                                                      model_annual_bounds))
  plots_year <- append(plots_year, list(plots_mb_cumul))
    
  
  #### . PLOT MASS BALANCE VERSUS ELEVATION BANDS ####
  ele_bands_plot_values <- getValues(data_dems$elevation_bands_plot[[elevation_grid_id]])
  ele_bands_plot_min <- min(ele_bands_plot_values, na.rm = T)
  ele_bands_plot_max <- max(ele_bands_plot_values, na.rm = T)
  ele_bands_plot_df <- data.frame(ele                 = seq(ele_bands_plot_min, ele_bands_plot_max, run_params$ele_bands_plot_size),
                                  ncells              = NA,
                                  mb_annual_meas_corr = NA,
                                  mb_annual_meas      = NA,
                                  mb_annual_hydro     = NA,
                                  mb_annual_fixed     = NA,
                                  mb_winter_fixed     = NA,
                                  mb_winter_meas      = NA)
  for (band_id in 1:length(ele_bands_plot_df[,1])) {
    band_cell_ids <- which(ele_bands_plot_values == ele_bands_plot_df$ele[band_id])
    ele_bands_plot_df$ncells[band_id] <- length(band_cell_ids)
    ele_bands_plot_df$mb_annual_meas_corr[band_id] <- mean(mb_meas_period_corr_values[band_cell_ids])
    ele_bands_plot_df$mb_annual_meas[band_id]      <- mean(getValues(massbal_annual_maps$meas_period)[band_cell_ids])
    ele_bands_plot_df$mb_annual_hydro[band_id]     <- mean(getValues(massbal_annual_maps$hydro)[band_cell_ids])
    ele_bands_plot_df$mb_annual_fixed[band_id]     <- mean(getValues(massbal_annual_maps$fixed)[band_cell_ids])
    ele_bands_plot_df$mb_winter_fixed[band_id]     <- mean(getValues(massbal_winter_maps$fixed)[band_cell_ids])
    if (process_winter) {
      ele_bands_plot_df$mb_winter_meas[band_id]     <- mean(getValues(massbal_winter_maps$meas_period)[band_cell_ids])
    }
  }
  
  ele_bands_plot_df_melt <- na.omit(melt(ele_bands_plot_df, id.vars = c("ele", "ncells"))) # This also removes the empty mb_winter_meas values if we don't have winter measurements.
  ele_bands_plot_df_melt$variable <- factor(ele_bands_plot_df_melt$variable, levels = c("mb_annual_fixed", "mb_annual_meas", "mb_annual_hydro",
                                                                                        "mb_winter_fixed", "mb_winter_meas",
                                                                                        "mb_annual_meas_corr"))
  
  base_size <- 16 # For the plot
  theme_elebands_plot <- theme_bw(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(face = "bold"),
          panel.grid = element_blank(),
          legend.position = c(0.4,0.85),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.title = element_blank())
  
  dat_ncells <- data.frame(ele    = c(ele_bands_plot_df$ele[1] - rep(25,2), rep(ele_bands_plot_df$ele, each = 2) + 25),
                           ncells = c(0, rep(ele_bands_plot_df$ncells, each = 2), 0))
  
  # Below: use geom_rect() instead of geom_polygon_pattern()
  # if there is a problem with package ggpattern.
  plot_elebands <-  ggplot(ele_bands_plot_df_melt) +
    # geom_rect(data = ele_bands_plot_df, aes(xmin = ele - run_params$ele_bands_plot_size/2, xmax = ele + run_params$ele_bands_plot_size/2, ymin = min(ele_bands_plot_df_melt$value) / 1e3, ymax = ncells * (max(ele_bands_plot_df_melt$value) - min(ele_bands_plot_df_melt$value)) / (1e3 * 4 * max(ncells)) + min(ele_bands_plot_df_melt$value) / 1e3)) +
    geom_polygon_pattern(data = dat_ncells,
                         aes(x = ele, y = ncells * (max(ele_bands_plot_df_melt$value) - min(ele_bands_plot_df_melt$value)) / (1e3 * 4 * max(ncells)) + min(ele_bands_plot_df_melt$value) / 1e3),
                         fill = "#FFFFFF", color = "#000000",
                         pattern_fill = "#000000", pattern_colour = "#000000",
                         pattern_angle = 35, pattern_size = 0.1, pattern_spacing = 0.02, pattern_density = 0.05) +
    geom_hline(yintercept = 0) +
    geom_line(aes(x = ele, y = value / 1e3, color = variable), size = 1) +
    scale_color_manual(breaks = c("mb_annual_fixed", "mb_annual_meas", "mb_annual_hydro", "mb_winter_fixed", "mb_winter_meas", "mb_annual_meas_corr"),
                       values = c("#8C00D4", "#FF0000", "#FF9000", "#0000FF", "#8080FF", "#000000"),
                       labels = c("Annual, fixed dates", "Annual, measurement period", "Annual, hydrological year",
                                  "Winter, fixed dates", "Winter, measurement period", "Annual, final")) +
    scale_y_continuous(breaks = pretty(ele_bands_plot_df_melt$value / 1e3), expand = expansion(0,0),
                       sec.axis = sec_axis(~ (. - min(ele_bands_plot_df_melt$value/1e3)) * 4 * max(ele_bands_plot_df$ncells) / ((max(ele_bands_plot_df_melt$value) - min(ele_bands_plot_df_melt$value))/1e3)  )) +
    scale_x_continuous(expand = expansion(0,0)) +
    coord_flip() +
    xlab("Elevation [m a.s.l.]") +
    ylab("Mass balance [m w.e.]") +
    theme_elebands_plot
  plots_year <- append(plots_year, list(plot_elebands))
  
  
  
  # Write multi-page PDF for the current year.
  plots_year_out <- ggarrange(plotlist = plots_year, ncol = 1, nrow = 1, align = "hv")
  ggexport(plots_year_out, filename = file.path(run_params$output_dirname, paste0(year_cur, ".pdf")), width = 21 * run_params$size_mult, height = 29.7 * run_params$size_mult)
  
  
  
  # Save the plot of the final mass balance of the year,
  # to be put in a PDF file with 1 plot per year.
  overview_areaplots[[year_id]] <- plots_year[[3]]
  

  
  
  
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

# save.image("zz.RData")
# load("zz.RData")

#### Overview plots ####
df_overview$mb_cumul <- cumsum(df_overview$mb_annual_meas_corr)
func_plot_overview(df_overview)

# Save to file the maps of final annual mass balance.
overview_areaplot <- ggarrange(plotlist = overview_areaplots, ncol = 1, nrow = 1, align = "hv")
ggexport(overview_areaplot, filename = file.path(run_params$output_dirname, "overview_areaplot.pdf"), width = 21 * run_params$size_mult, height = 29.7 * run_params$size_mult)
