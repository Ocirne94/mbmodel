###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the routine to compute the initial snow cover for a year.    #
################################################################################################### 


# Algorithm:
  # Compute avalanche on the topographic snow distribution (with appropriate multiplier for max deposition)
  # Reduce effect of the computed grid (custom reduction parameter, in IDL it is 0.5) and normalize the resulting grid over the glacier surface
  # Compute probes_idw, normalize the probes idw output
  # Multiply the two grids (probes idw and avalanched topographic dist) and normalize again over the glacier surface
  # Compute grid with snowgrad and snow line elevation
  # Combine with the distribution grid and return

func_compute_initial_snow_cover <- function(run_params,
                                            data_dhms,
                                            data_dems,
                                            grids_snowdist_topographic,
                                            grids_avalanche,
                                            grid_id,
                                            data_massbal_winter) {
  
  grids_avalanche_cur <- sapply(grids_avalanche, `[[`, grid_id)
  
  dist_initial <- grids_snowdist_topographic[[grid_id]]
  
  # We set the mass deposition limit to 1 (the input
  # snowdist_topographic grid is normalized, so it has
  # values close to 1; deposition of 1 leads to a realistic
  # avalanche deposit for Barkrak glacier).
  dist_post_avalanche <- func_avalanche(grids_avalanche_cur, dist_initial, 1 / run_params$deposition_mass_lim)
  
  writeRaster(dist_post_avalanche, "dist_post_avalanche.tif", overwrite = T)
  
  # Reduce importance of topographic/avalanche distribution variability.
  dist_red <- 1 + run_params$initial_snow_dist_red_fac * (dist_post_avalanche - 1)
  
  # Not sure we should renormalize this, as it would kill the effect of avalanches
  # reaching the glacier from ice-free slopes! Plus the grids used so far were already normalized.
  # dist_red <- dist_red / mean(dist_red[data_dems$glacier_cell_ids[[grid_id]]])
  
  writeRaster(dist_red, "dist_red.tif", overwrite = T)
  
  dist_cur <- dist_red
  
  # If we have any winter stakes for the year,
  # use them to correct the 
  if(length(data_massbal_winter[,1]) > 0) {
    
    dist_probes_idw <- func_snow_probes_idw(run_params, data_massbal_winter, data_dhms)
    dist_probes_idw_norm <- dist_probes_idw / mean(dist_probes_idw[data_dems$glacier_cell_ids[[grid_id]]])
    
    dist_corr <- dist_cur * dist_probes_idw_norm
    
    # Not sure we should renormalize this one!
    # The probes are already normalized, and for the
    # topographic/avalanche grids see the comment after
    # the computation of dist_red.
    # dist_corr <- dist_corr / mean(dist_corr[data_dems$glacier_cell_ids[[grid_id]]])
    
    dist_cur <- dist_corr
    
  }
  
  # Distribution from snow line elevation and snow gradient.
  dist_snl <- setValues(dist_cur,
                        pmax(0,
                             getValues(data_dems$elevation[[grid_id]] - run_params$initial_snowline_elevation) * run_params$initial_snow_gradient / 100))

  writeRaster(dist_snl, "dist_snl.tif", overwrite = T)
  
  # Compute final distribution.
  dist_final <- dist_snl * dist_cur
  
  writeRaster(dist_final, "dist_final.tif", overwrite = T)
  
  return(dist_final)
  
}
