###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the computation of a grid of relative (normalized)           #
#                 snow distribution from elevation, slope and curvature.                          #
# Latest update:  2021.1.7                                                                        #
###################################################################################################


func_compute_snowdist_topographic <- function(run_params, dem, dhm) {
  
  #### Curvature factor (lower accumulation on convex surfaces, higher on concave) ####
  # In the original IDL implementation, the curvature multiplication factor
  # is (1 - x), with x linearly dependent on the terrain curvature up to a cutoff:
  # curvature values at or exceeding ±threshold produce an x = ±0.5.
  # threshold is currently the smaller of the two curvature extremes (abs(max) and abs(min))
  # reduced by a tunable factor (default 1.2), so that both 1.5 and 0.5 are reached (even with some margin!)
  # in the curvature multiplication factor (arbitrary!).
  # IMPORTANT: terrain curvature in IDL is computed manually by taking cells at distance 3, 5 and 6
  # from the focal cell. This is a bit arbitrary (taken from some old ArcGIS documentation)
  # and produces a smoothed-out curvature.
  
  # We use a smoothed DEM to compute curvature because it is very sensitive to DEM noise.
  # The window size used for the smoothing is automatically computed from the smoothing amount.
  # Then for the snow distribution we use the same approach as the IDL version.
  dhm_smooth <- raster.gaussian.smooth(dhm,
                                       run_params$curvature_dhm_smooth,
                                       run_params$dhm_smooth_windowsize,
                                       type = "sum")
  dhm_curvature <- curvature(dhm_smooth, "total")
  
  # Compute curvature cutoff.
  dhm_curvature_bound <- min(abs(max(dhm_curvature[], na.rm=T)), abs(min(dhm_curvature[], na.rm=T))) / run_params$curvature_cutoff_fact
  
  # Apply corrected cutoff to curvature.
  dhm_curvature_cut <- setValues(dhm_curvature,
                                 pmin(dhm_curvature_bound, pmax(-dhm_curvature_bound, dhm_curvature[])))
  
  # Compute final curvature factor for snow distribution.
  snowdist_curv_mult <- 1 - (dhm_curvature_cut * run_params$curvature_effect_limit / dhm_curvature_bound)
  
  
  #### Elevation factor (lower accumulation at very high elevations) ####
  snowdist_ele_mult <- setValues(dhm, 1)
  ele_scaling_fac <- (max(dhm[], na.rm=T) - run_params$elevation_effect_threshold) / run_params$elevation_effect_fact
  dhm_ids_high <- which(dhm[] > run_params$elevation_effect_threshold)
  
  snowdist_ele_mult[dhm_ids_high] <- 2 - 2^( (dhm[dhm_ids_high] - run_params$elevation_effect_threshold) / ele_scaling_fac )
  
  
  #### Slope factor (avalanches) ####
  # Here we follow Gruber (2007), doi:10.1029/2006WR004868.
  dhm_slope <- terrain(dhm, "slope", unit = "degrees")
  
  
  snowdist_topographic <- snowdist_curv_mult * snowdist_ele_mult * snowdist_slope_mult
  snowdist_topographic_norm <- snowdist_topographic / mean(snowdist_topographic[which(!is.na(dem[]))])
  
  return(snowdist_topographic_norm)
  
}
