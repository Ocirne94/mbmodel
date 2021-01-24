###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the routine to call the mass balance model and compare       #
#                 its output to the measured stakes, computing RMS and BIAS of the mass balance.  #
#                 This routine is set up to be run in the optim() function (i.e., it does not     #
#                 save or return the actual mass balance grids or series, just the values needed  #
#                 for optimization).                                                              #
###################################################################################################

func_run_model_year_optim <- function(year_param_multipliers,
                                      run_params, year_cur_params, grid_id, data_dhms, data_dems, data_surftype,
                                      snowdist_init, data_radiation, weather_series_cur, dist_topographic_values_red,
                                      dist_probes_norm_values_red, grids_avalanche_cur, dx1, dx2, dy1, dy2,
                                      nstakes, model_days_n, massbal_meas_cur) {
  
  
  #### . .  APPLY MULTIPLIERS FOR OPTIMIZATION ####
  year_cur_params_multiplied <- year_cur_params
  year_cur_params_multiplied$melt_factor <- year_cur_params$melt_factor * year_param_multipliers[1]
  year_cur_params_multiplied$rad_fact_ice <- year_cur_params$rad_fact_ice * year_param_multipliers[2]
  year_cur_params_multiplied$rad_fact_snow <- year_cur_params$rad_fact_snow * year_param_multipliers[3]
  year_cur_params_multiplied$prec_corr <- year_cur_params$prec_corr * year_param_multipliers[4]
  year_cur_params_multiplied$prec_summer_fact <- year_cur_params$prec_summer_fact * year_param_multipliers[5]
  year_cur_params_multiplied$prec_elegrad <- year_cur_params$prec_elegrad * year_param_multipliers[6]
  
  
  #### . .  RUN MASS BALANCE MODEL ####
  mb_model_output <- func_massbal_model(run_params,
                                        year_cur_params_multiplied,
                                        getValues(data_dhms$elevation[[grid_id]]),
                                        data_dems$glacier_cell_ids[[grid_id]],
                                        getValues(data_surftype[[grid_id]]),
                                        getValues(snowdist_init),
                                        data_radiation,
                                        weather_series_cur,
                                        dist_topographic_values_red,
                                        dist_probes_norm_values_red,
                                        grids_avalanche_cur)
  
  
  #### . .  COMPARE TO STAKE MEASUREMENTS ####
  # Extract the whole modeled series for all stakes.
  stakes_series_mod_all <- func_extract_modeled_stakes(run_params,
                                                       dx1, dx2, dy1, dy2,
                                                       mb_model_output$vec_massbal_cumul,
                                                       nstakes,
                                                       model_days_n)
  
  
  # Find indices of the days corresponding to the stake measurements.
  # We match w.r.t. weather_series_cur whose index is off by ~0.5 with the
  # mass balance (mb_model_out$gl_massbal_cumul[1] is the initial condition
  # (i.e. 0.0) at 00:00 of the first day, then the index of the weather series
  # corresponds to the full following 24 hours, then gl_massbal_cumul[2] is the
  # cumulative mass balance by the end of that same day.
  # So it would be equally correct to also shift all the day indices by one (little to no change).
  # We also find the start date for stakes set to NA (i.e. start date = date of mass balance minimum).
  annual_stakes_start_ids <- pmatch(massbal_meas_cur$start_date,
                                    weather_series_cur$timestamp,
                                    duplicates.ok = TRUE)
  annual_stakes_end_ids   <- pmatch(massbal_meas_cur$end_date,
                                    weather_series_cur$timestamp,
                                    duplicates.ok = TRUE)
  annual_stakes_start_ids_corr <- func_compute_unknown_stakes_start_ids(run_params, annual_stakes_start_ids, weather_series_cur, stakes_series_mod_all)
  
  
  # Cumulative mass balance of each stake
  # over the measurement period (numeric vector).
  stakes_mb_mod  <- as.numeric(stakes_series_mod_all)[((1:nstakes_annual)-1)*(model_days_n+1) + annual_stakes_end_ids] -
                    as.numeric(stakes_series_mod_all)[((1:nstakes_annual)-1)*(model_days_n+1) + annual_stakes_start_ids_corr]
  # Corresponding measurement.
  stakes_mb_meas <- massbal_meas_cur$dh_cm * massbal_meas_cur$density * 10 # 10: cm w.e. to mm w.e.
  
  stakes_bias <- stakes_mb_mod - stakes_mb_meas
  
  year_bias <- mean(stakes_bias)
  year_rms  <- sqrt(mean(stakes_bias^2))
  
  cat("BIAS:", year_bias, "mm w.e.\n")
  cat("RMS:", year_rms, "mm w.e.\n")
  
  # score_cur <- (abs(year_bias) + 500) * year_rms # Minimizing this should minimize both BIAS and RMS. 500 (mm w.e.) is arbitrary!
  score_cur <- year_rms
  
  return(score_cur)
  
}
