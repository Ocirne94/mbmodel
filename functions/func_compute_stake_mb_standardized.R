###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the routine to compute a "standardized" stake mass balance,  #
#                 which is the sum of the measured mass balance and the mass balance simulated    #
#                 at the stake over the "measurement period" but outside of the stake observation #
#                 period. The result is the best estimate for the mass balance at the stake       #
#                 over the "measurement period". The resulting "standardized" mass balances       #
#                 are comparable even if the stakes were measured at different points in time.    #
###################################################################################################

# For each stake we compute its cumulative mass balance over the whole "measurement period"
# (i.e. earliest stake start - latest stake end). We do this by adding to each stake measurement
# the MODELED mass balance at the stake on the days which are within the "measurement period"
# but before the individual stake start / after the individual stake end.
# This is important to have a standardized measured value, which can be compared to the simulated value
# over the "measurement period"; for each year, this period includes the earliest stake start and the latest stake end.
# The standardized value is the one which is plotted in the "measurement period" maps
# and also in the scatterplot of vertical mass balance distribution.
# If all the stakes starts are the same, and all the stake ends are the same, the standardized stake measurement
# is the same as the original stake measurement.

func_compute_stake_mb_standardized <- function(mod_output_annual_cur,
                                               massbal_annual_meas_cur,
                                               nstakes_annual) {
  
  id_measperiod_start <- min(mod_output_annual_cur$stakes_start_ids_corr)
  id_measperiod_end   <- max(mod_output_annual_cur$stakes_end_ids)
  
  massbal_standardized <- rep(NA, nstakes_annual)
  for (stake_id in 1:nstakes_annual) {
    stake_mod_mb_measperiod_start <- mod_output_annual_cur$stakes_series_mod_all[id_measperiod_start,stake_id]
    stake_mod_mb_stake_start <- mod_output_annual_cur$stakes_series_mod_all[mod_output_annual_cur$stakes_start_ids_corr[stake_id],stake_id]
    stake_mod_mb_stake_end <- mod_output_annual_cur$stakes_series_mod_all[mod_output_annual_cur$stakes_end_ids[stake_id],stake_id]
    stake_mod_mb_measperiod_end <- mod_output_annual_cur$stakes_series_mod_all[id_measperiod_end,stake_id]
    
    stake_standard_corr <- (stake_mod_mb_stake_start - stake_mod_mb_measperiod_start) + (stake_mod_mb_measperiod_end - stake_mod_mb_stake_end)
    massbal_standardized[stake_id] <- massbal_annual_meas_cur$massbal[stake_id] + stake_standard_corr
  }
  
  return(massbal_standardized)
  
}
