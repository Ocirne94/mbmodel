###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the code to post-process the simulated mass balance.         #
#                 Specifically, here are performed:                                               #
#                   (1) mass balance correction in elevation bands                                #
#                   (2) computation of ELA and AAR                                                #
#                   (3) standardization of stake measurements to the whole measurement period     #
###################################################################################################

# NOTE: this code is source()'d as part of main.R.
# We put code here just to make it more organized.


#### Correct mass balance bias according to user-defined elevation bands ####
massbal_annual_maps$meas_period_corr <- func_correct_massbal_elebands(year_cur_params,
                                                                      data_dems,
                                                                      massbal_annual_meas_cur,
                                                                      mod_output_annual_cur,
                                                                      massbal_annual_maps,
                                                                      dem_grid_id)

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



#### Compute ELA and AAR ####
ela_aar <- func_compute_ela_aar(run_params,
                                mb_meas_period_corr_values,
                                data_dems,
                                dem_grid_id)


#### Compute standardized stake measurements ####
massbal_annual_meas_cur$massbal_standardized <- func_compute_stake_mb_standardized(mod_output_annual_cur,
                                                                                   massbal_annual_meas_cur,
                                                                                   nstakes_annual) 
