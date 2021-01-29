###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the routines to optimize the annual mass balance (currently  #
#                 acting on the melt factor and the radiation factors).                           #
#                 NOTE: several functions are defined in this file, to call the annual simulation #
#                 while optimizing either the bias or the rms.                                    #
################################################################################################### 

# Annual optimization:
# cancel the bias by altering the melt factor
# and radiation factor together, by the same amount.

# We have a list of additive corrections (see func_run_simulation_single)
# which are applied to the year_cur_params before the
# model run, so that we can tune parameters for the optimization.
# The variable specifying the amount of correction is
# corr_fact_cur, which multiplies the original parameter value to
# get the additive correction. So a corr_fact_cur of -1 means
# that the melt factor and radiation factors become 0,
# a corr_fact_cur of +1 means that they get doubled.
# This approach works seamlessly for the optimization
# of winter and annual mass balance, and allows easy
# optimization of other parameters.
# For the moment we just do the same optimization
# as the original IDL implementation, on Barkrak glacier
# there is almost no RMS gain (< 2 mm w.e. out of 500)
# when optimizing independently the radiation factors.


#### CONVENIENCE FUNCTION TO RUN uniroot() ####
func_optim_annual_worker <- function(corr_fact_cur,
                                     run_params, year_cur_params, elevation_grid_id, surftype_grid_id,
                                     data_dhms, data_dems, data_surftype,
                                     snowdist_init, data_radiation, weather_series_cur, dist_topographic_values_red,
                                     dist_probes_norm_values_red, grids_avalanche_cur,
                                     dx1, dx2, dy1, dy2,
                                     nstakes, model_days_n, massbal_meas_cur, stakes_cells) {
  
  annual_corrections_cur <- list(melt_factor  = corr_fact_cur * year_cur_params$melt_factor,
                                 rad_fact_ice = corr_fact_cur * year_cur_params$rad_fact_ice)
  
  # We already set this in the global environment
  # so that we don't have to re-run the model with
  # the optimized parameters to get the actual
  # model output.
  mod_output_annual_cur <<- func_run_simulation_single(annual_corrections_cur,
                                           run_params, year_cur_params, elevation_grid_id, surftype_grid_id,
                                           data_dhms, data_dems, data_surftype,
                                           snowdist_init, data_radiation, weather_series_cur, dist_topographic_values_red,
                                           dist_probes_norm_values_red, grids_avalanche_cur,
                                           dx1, dx2, dy1, dy2,
                                           nstakes, model_days_n, massbal_meas_cur, stakes_cells)
  bias <- mod_output_annual_cur$global_bias
  
  return(bias * (abs(bias) >= run_params$optim_bias_threshold))
}


#### ACTUAL OPTIMIZATION FUNCTION ####
func_optimize_mb_annual <- function(run_params, year_cur_params, elevation_grid_id, surftype_grid_id,
                                    data_dhms, data_dems, data_surftype,
                                    snowdist_init, data_radiation, weather_series_cur, dist_topographic_values_red,
                                    dist_probes_norm_values_red, grids_avalanche_cur, dx1, dx2, dy1, dy2,
                                    nstakes, model_days_n, massbal_meas_cur, stakes_cells) {

  # extendInt lets the root-finding work even if the original search interval is wrong.
  # "downX" is based on the fact that the bias is becoming more negative as the additive
  # correction to the melt and radiation factors becomes more positive.
  optim_res <- uniroot(func_optim_annual_worker,
                       
                       run_params = run_params, year_cur_params = year_cur_params,
                       elevation_grid_id = elevation_grid_id, surftype_grid_id = surftype_grid_id,
                       data_dhms = data_dhms, data_dems = data_dems, data_surftype = data_surftype,
                       snowdist_init = snowdist_init, data_radiation = data_radiation,
                       weather_series_cur = weather_series_cur, dist_topographic_values_red = dist_topographic_values_red,
                       dist_probes_norm_values_red = dist_probes_norm_values_red,
                       grids_avalanche_cur = grids_avalanche_cur, dx1 = dx1, dx2 = dx2, dy1 = dy1, dy2 = dy2,
                       nstakes = nstakes_annual, model_days_n = model_days_n, massbal_meas_cur = massbal_meas_cur,
                       stakes_cells = stakes_cells,
                       
                       lower = -min(year_cur_params$melt_factor, year_cur_params$rad_fact_ice),
                       upper = run_params$optim_max_corr_fact,
                       tol = 0.01,
                       extendInt = "downX")
  
  annual_corrections_best <- list(melt_factor  = optim_res$root * year_cur_params$melt_factor,
                                  rad_fact_ice = optim_res$root * year_cur_params$rad_fact_ice)
  
  return(annual_corrections_best)
  
}
