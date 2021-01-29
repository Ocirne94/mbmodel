###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the routines to optimize the annual mass balance (currently  #
#                 acting on the melt factor and the radiation factors).                           #
#                 Optimization is performed by computing the bias derivative w.r.t. the           #
#                 correction factor, since the bias is actually quasi-linear (weak albedo         #
#                 feedback). Then we can converge quickly to zero bias.                           #
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


#### CONVENIENCE FUNCTION TO RUN THE OPTIMIZATION ON ####
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
  
  return(bias)
}




#### ACTUAL OPTIMIZATION FUNCTION ####
func_optimize_mb_annual <- function(run_params, year_cur_params, elevation_grid_id, surftype_grid_id,
                                    data_dhms, data_dems, data_surftype,
                                    snowdist_init, data_radiation, weather_series_cur, dist_topographic_values_red,
                                    dist_probes_norm_values_red, grids_avalanche_cur, dx1, dx2, dy1, dy2,
                                    nstakes, model_days_n, massbal_meas_cur, stakes_cells) {
  
  cat("** Mass balance optimization **\n")
  cat("\n* Model run # 1\n")
  corr_fact_prev <- 0
  bias_prev <- func_optim_annual_worker(corr_fact_prev,
                                        run_params, year_cur_params, elevation_grid_id, surftype_grid_id,
                                        data_dhms, data_dems, data_surftype,
                                        snowdist_init, data_radiation, weather_series_cur, dist_topographic_values_red,
                                        dist_probes_norm_values_red, grids_avalanche_cur,
                                        dx1_annual, dx2_annual, dy1_annual, dy2_annual,
                                        nstakes_annual, model_days_n, massbal_annual_meas_cur, annual_stakes_cells)
  
  cat("\n* Model run # 2\n")
  # This 0.01 increment is arbitrary, we just need
  # a small interval to approximate the bias
  # derivative with a finite difference.
  # A very small value is safer in case the starting
  # value of the factors was very low
  # (we don't want to go to the negatives!).
  corr_fact_cur <- 0.01
  bias_cur <- func_optim_annual_worker(corr_fact_cur,
                                       run_params, year_cur_params, elevation_grid_id, surftype_grid_id,
                                       data_dhms, data_dems, data_surftype,
                                       snowdist_init, data_radiation, weather_series_cur, dist_topographic_values_red,
                                       dist_probes_norm_values_red, grids_avalanche_cur,
                                       dx1_annual, dx2_annual, dy1_annual, dy2_annual,
                                       nstakes_annual, model_days_n, massbal_annual_meas_cur, annual_stakes_cells)
  
  niter <- 2
  while (abs(bias_cur) > run_params$optim_bias_threshold) {
    bias_slope <- (bias_cur - bias_prev) / (corr_fact_cur - corr_fact_prev)
    bias_prev <- bias_cur
    corr_fact_prev <- corr_fact_cur
    corr_fact_cur <- corr_fact_cur - (bias_cur / bias_slope) # Apply linear correction with the computed derivative.
    niter <- niter + 1
    cat("\n* Model run #", niter, "\n")
    bias_cur <- func_optim_annual_worker(corr_fact_cur,
                                         run_params, year_cur_params, elevation_grid_id, surftype_grid_id,
                                         data_dhms, data_dems, data_surftype,
                                         snowdist_init, data_radiation, weather_series_cur, dist_topographic_values_red,
                                         dist_probes_norm_values_red, grids_avalanche_cur,
                                         dx1_annual, dx2_annual, dy1_annual, dy2_annual,
                                         nstakes_annual, model_days_n, massbal_annual_meas_cur, annual_stakes_cells)
  }
  
  
  annual_corrections_best <- list(melt_factor  = corr_fact_cur * year_cur_params$melt_factor,
                                  rad_fact_ice = corr_fact_cur * year_cur_params$rad_fact_ice)
  
  return(annual_corrections_best)
  
  
}
