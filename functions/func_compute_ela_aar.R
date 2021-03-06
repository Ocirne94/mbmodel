###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the routine to compute the equilibrium line altitude and     #
#                 the accumulation-area ratio (ELA and AAR) after modeling one year.              #
################################################################################################### 


# The equilibrium line altitude is computed by classifying
# the glacier grid into elevation bands (with user-defined
# vertical extent) and then taking the band whose mean mass
# balance over the corrected measured period is closest to 0.
# NOTE: these bands (typically 10 vertical meters wide)
# are NOT the same as for the correction of modeled mass
# balance in elevation bands (glacier-dependent bands),
# and NOT the same as for the plot of modeled mass balance
# vs elevation bands (typically 40 vertical meters wide).
func_compute_ela_aar <- function(run_params,
                                 mb_meas_period_corr_values,
                                 data_dems,
                                 dem_grid_id) {
  
  ele_bands_values <- getValues(data_dems$elevation_bands_ela[[dem_grid_id]])
  ele_bands_min <- min(ele_bands_values, na.rm = T)
  ele_bands_max <- max(ele_bands_values, na.rm = T)
  ele_bands_df <- data.frame(ele = seq(ele_bands_min, ele_bands_max, run_params$ele_bands_ela_size),
                             mb_corr = NA)
  for (band_id in 1:length(ele_bands_df[,1])) {
    ele_bands_df$mb_corr[band_id] <- mean(mb_meas_period_corr_values[ele_bands_values == ele_bands_df$ele[band_id]], na.rm=T)
  }
  ela_band_id <- which.min(abs(ele_bands_df$mb_corr))
  ela <- ele_bands_df$ele[ela_band_id]
  
  aar <- length(which(mb_meas_period_corr_values >= 0)) / length(data_dems$glacier_cell_ids[[dem_grid_id]])
  
  return(c(ela = ela, aar = aar))
  
}
