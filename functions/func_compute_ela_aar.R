###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the routine to compute the equilibrium line altitude and     #
#                 the accumulation-area ratio (ELA and AAR) after modeling one year.              #
################################################################################################### 

func_compute_ela_aar <- function(run_params,
                                 massbal_annual_maps,
                                 data_dems) {
  
  mb_meas_period_corr_values <- getValues(massbal_annual_maps$meas_period_corr)
  ele_bands_values <- getValues(data_dems$elevation_bands[[elevation_grid_id]])
  ele_bands_min <- min(ele_bands_values, na.rm = T)
  ele_bands_max <- max(ele_bands_values, na.rm = T)
  ele_bands_df <- data.frame(ele = seq(ele_bands_min, ele_bands_max, run_params$ele_bands_size),
                             mb_corr = NA)
  for (band_id in 1:length(ele_bands_df[,1])) {
    ele_bands_df$mb_corr[band_id] <- mean(mb_meas_period_corr_values[ele_bands_values == ele_bands_df$ele[band_id]], na.rm=T)
  }
  ela_band_id <- which.min(abs(ele_bands_df$mb_corr))
  ela <- ele_bands_df$ele[ela_band_id]
  
  aar <- length(which(mb_meas_period_corr_values >= 0)) / length(data_dems$glacier_cell_ids[[elevation_grid_id]])
  
  return(c(ela = ela, aar = aar))
  
}
