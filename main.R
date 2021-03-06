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

boot_file_write   <- FALSE                # Save .RData file with the input data, for faster reload.
boot_file_read    <- TRUE                 # Load .RData file with the input data, instead of loading input files.
boot_file_name    <- "boot_file_barkrak.RData"    # Name of the .RData input data file.


#### Load function definitions and R modules ####
invisible(sapply(paste("functions/", list.files("functions/", pattern = "\\.R$"), sep=""), source))
source("func_set_params.R")
# Remove cacheDir option to force recompilation of the C++ code
# (useful after changing computer or editing the source file).
sourceCpp("functions/func_avalanche_gruber.cpp", cacheDir = "functions/")
func_load_libraries()


source(file.path("procedures", "pro_load_data_parameters.R"))    # Load data and parameters.
source(file.path("procedures", "pro_compute_grid_parameters.R")) # Compute grid parameters.
source(file.path("procedures", "pro_compute_all_fixed_grids.R")) # Compute fixed grids.
source(file.path("procedures", "pro_save_boot_files.R"))         # Save boot files.
source(file.path("procedures", "pro_setup_loop.R"))              # Prepare variables before main loop.


# Create output directory.
dir.create(file.path(run_params$output_dirname, "annual_results"), recursive = TRUE)

#### Main loop ####
for (year_id in 1:run_params$n_years) {

  #### . Select current year, parameters, data. ####
  year_cur <- run_params$years[year_id]
  year_cur_params <- func_load_year_params(run_params, year_cur)
  cat("\n\n\n\n============  STARTING NEW YEAR:", year_cur, " ============\n")
  source(file.path("procedures", "pro_select_year_data.R"))
  
  source(file.path("procedures", "pro_find_stake_dxdy.R")) # Find stake offsets on the grid.
  
  
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
  
  
  #### .  COMPUTE INITIAL SNOW COVER ####
  source(file.path("procedures", "pro_compute_initial_snow_cover.R"))
  
  
  #### .  WINTER MASS BALANCE, to optimize the precipitation correction ####
  # We set this here so that there is no correction
  # if we don't do the winter optimization.
  corr_fact_winter      <- 0
  # We set this to NULL to have it defined (for the
  # extraction functions) in case we don't do winter processing.
  mod_output_winter_cur <- NULL
  if (process_winter)  {
    
    # We ask duplicates = FALSE, else the bilinear filtering in func_extract_modeled_stakes()
    # can fail when a stake is exactly at the same (X and/or Y) coordinate as a cell center.
    # duplicates = FALSE returns four different cells. In case we have a stake exactly
    # aligned with a cell center, unless we are at the lower raster border (which we should
    # always avoid!) the additional cells returned with duplicates = FALSE (cells which would
    # not be part of the actual adjacent cells) have higher index than the "true" adjacent cells.
    winter_stakes_cells <- rowSort(fourCellsFromXY(data_dhms$elevation[[elevation_grid_id]], as.matrix(massbal_winter_meas_cur[,4:5]), duplicates = FALSE))
    
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
  annual_stakes_cells <- rowSort(fourCellsFromXY(data_dhms$elevation[[elevation_grid_id]], as.matrix(massbal_annual_meas_cur[,4:5]), duplicates = FALSE))
  
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
  
  
  #### . EXTRACT MAPS OF CUMULATIVE MASS BALANCE AT DATES OF INTEREST ####
  source(file.path("procedures", "pro_extract_massbalance.R"))
  
  
  #### . CORRECT ANNUAL MASS BALANCE IN ELEVATION BANDS ####
  massbal_annual_maps$meas_period_corr <- func_correct_massbal_elebands(year_cur_params,
                                                                        data_dems,
                                                                        massbal_annual_meas_cur,
                                                                        mod_output_annual_cur,
                                                                        massbal_annual_maps)
  
  massbal_annual_values <- sapply(massbal_annual_maps, cellStats, stat = "mean", na.rm = TRUE)
  massbal_winter_values <- sapply(massbal_winter_maps, cellStats, stat = "mean", na.rm = TRUE)
  
  # Extract the grid values of the "final" mass balance, since we will use them several times below.
  mb_meas_period_corr_values <- getValues(massbal_annual_maps$meas_period_corr)
  
  # Compute time series of glacier-wide mass balance,
  # including the bias correction in elevation bands.
  # We assign the correction to the melt component,
  # accumulation stays the same.
  id_measperiod_start <- min(mod_output_annual_cur$stakes_start_ids_corr)
  id_measperiod_end   <- max(mod_output_annual_cur$stakes_end_ids)
  mb_band_bias <- massbal_annual_values[["meas_period"]] - massbal_annual_values[["meas_period_corr"]]
  mb_band_corr_fact <- (mod_output_annual_cur$gl_melt_cumul[id_measperiod_end] - mod_output_annual_cur$gl_melt_cumul[id_measperiod_start] + mb_band_bias) / (mod_output_annual_cur$gl_melt_cumul[id_measperiod_end] - mod_output_annual_cur$gl_melt_cumul[id_measperiod_start])
  mod_output_annual_cur$gl_melt_cumul_bandcorr <- mod_output_annual_cur$gl_melt_cumul * mb_band_corr_fact
  mod_output_annual_cur$gl_massbal_cumul_bandcorr <- mod_output_annual_cur$gl_accum_cumul - mod_output_annual_cur$gl_melt_cumul_bandcorr
  
  
  
  #### . COMPUTE ELA and AAR ####
  ela_aar <- func_compute_ela_aar(run_params, mb_meas_period_corr_values, data_dems)
  
  
  #### . COMPUTE STANDARDIZED STAKE MEASUREMENTS ####
  massbal_annual_meas_cur$massbal_standardized <- func_compute_stake_mb_standardized(mod_output_annual_cur,
                                                                                     massbal_annual_meas_cur,
                                                                                     nstakes_annual)


  #### . Save overview values for the year ####
  source(file.path("procedures", "pro_save_overview_values.R"))
  
  
  #### . Produce all plots for the year ####
  source(file.path("procedures", "pro_plot_year.R"))
  
  
  #### . WRITE MODEL OUTPUT TO FILES ####
  func_write_year_output(run_params,
                         year_cur,
                         massbal_annual_maps,
                         massbal_winter_maps,
                         mod_output_annual_cur,
                         massbal_annual_meas_cur,
                         model_time_bounds,
                         ele_bands_plot_df,
                         process_winter)
  
  
  
  if (max(abs((extract(massbal_annual_maps$meas_period, cbind(massbal_annual_meas_cur$x, massbal_annual_meas_cur$y), method = "bilinear") - massbal_annual_meas_cur$massbal_standardized) - (mod_output_annual_cur$stakes_mb_mod - mod_output_annual_cur$stakes_mb_meas))) > 1e-5) {
    stop("VERY BAD ERROR: the recomputed stake mass balance biases over the stake period and over the single \"measurement period\" do not match. Probably an issue with the manual bilinear filtering of the stakes series. Check if there are stakes coordinates exactly aligned with cell centers, they are likely the cause.")
    Sys.sleep(1e9)
  }
  
}


#### Overview plots ####
df_overview$mb_cumul <- cumsum(df_overview$mb_annual_meas_corr)
func_plot_overview(df_overview)
df_overview_out <- data.frame(year = df_overview$year,
                              apply(df_overview[,2:7], 2, sprintf, fmt="%.3f"),
                              df_overview[,8],
                              sprintf("%.1f", df_overview[,9]),
                              apply(df_overview[,10:13], 2, sprintf, fmt="%.3f"),
                              df_overview[,14],
                              sprintf("%.1f", df_overview[,15]))
names(df_overview_out) <- names(df_overview)
write.csv(df_overview_out,
          file.path(run_params$output_dirname, "overview.csv"),
          quote = FALSE,
          row.names = FALSE)

# Save to file the maps of final annual mass balance.
overview_areaplot <- ggarrange(plotlist = overview_areaplots, ncol = 1, nrow = 1, align = "hv")
suppressMessages(ggexport(overview_areaplot,
         filename = file.path(run_params$output_dirname, "overview_areaplot.pdf"),
         width = 21 * run_params$size_mult,
         height = 29.7 * run_params$size_mult))
