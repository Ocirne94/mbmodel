###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the loading routine for the DEM(s) and DHM(s). As output     #
#                 we get a list with one raster per year (selected from the closest input         #
#                 available, or interpolated if asked to do so).                                  #
#                 NOTE: we should never modify the list elements, rather work on copies.          #
# Latest update:  2021.1.6                                                                        #
###################################################################################################


# ALGORITHM:
# First select either the DEMs or the DHMs (we re-use the code, exact same processing): load_which can be "dem" or "dhm".
# Then: if we have only one dem just use it,
# else,
#   if dem_interpolate is false, for each modeled year find the closest dem year and use its dem.
#   If dem_interpolate is true, for each modeled year look for the year's dem, if you have it use it,
#     if not look for the two closest enclosing years,
#       if both found do linear interpolation,
#       if only one found (i.e. modeling a year before the earliest dem or after the last one) just use it.

func_load_elevation_grids <- function(run_params, load_which) {
  
  # Here we will put the output.
  grids_out <- list()
  
  
  # Work on DEMs or DHMs?
  if (load_which == "dem") {
    # Set variables to DEMs.
    grid_paths <- paste(run_params$dir_data_dem,
                        run_params$filename_dem_prefix,
                        run_params$dem_years,
                        run_params$filename_dem_suffix,
                        sep = "")
    grid_interpolate <- run_params$dem_interpolate
    grid_years <- run_params$dem_years
    
  } else {
    # Set variables to DHMs.
    grid_paths <- paste(run_params$dir_data_dhm,
                        run_params$filename_dhm_prefix,
                        run_params$dhm_years,
                        run_params$filename_dhm_suffix,
                        sep = "")
    grid_interpolate <- run_params$dhm_interpolate
    grid_years <- run_params$dhm_years
  }
  
  
  # Do we have a single DEM/DHM? If so just use it every year.
  if (length(grid_years) == 1) {
    for (year_cur_id in 1:run_params$n_years) {
      grids_out[[year_cur_id]] <- raster(grid_paths[1])
    }
    
    # We have more than a single DEM/DHM!
  } else {
    
    if (grid_interpolate == FALSE) {
      # For each modeled year find the closest DEM year and use its DEM.
      for (year_cur_id in 1:run_params$n_years) {
        
        year_cur <- run_params$years[year_cur_id]
        grid_year_closest_id <- which.min(abs(grid_years - year_cur))
        grids_out[[year_cur_id]] <- raster(grid_paths[grid_year_closest_id])
        
      }
      
      # Here the case grid_interpolate == TRUE
    } else {
      # For each modeled year look for a DEM exactly from that year,
      # if found use it,
      # if not look for the two closest enclosing years,
      #     if both found do linear interpolation,
      #     if only one found (i.e. modeling a year before the earliest DEM or after the last one) just use it.
      for (year_cur_id in 1:run_params$n_years) {
        
        year_cur <- run_params$years[year_cur_id]
        grid_id_cur <- which(grid_years == year_cur)
        
        # Found DEM for the current year!
        if (length(grid_id_cur) != 0) { # This should be only 0 or 1, i.e. two DEMs for a single year are not allowed.
          
          grids_out[[year_cur_id]] <- raster(grid_paths[grid_id_cur])
          
        # No DEM for the current year, we have to interpolate if we can.  
        } else {
          
          year_dist <- grid_years - year_cur
          
          # Check if we have two enclosing years, or if instead we are outside the range of the DEM years.
          if (max(year_dist) * min(year_dist) < 0) {
            # We can interpolate! Find the enclosing DEMs/DHMs.
            grid_year_earlier_id <- which.max(year_dist[which(year_dist < 0)]) # year_dist is a sorted vector, so this indexing should work.
            grid_year_later_id <- grid_year_earlier_id + 1
            grid_year_earlier <- grid_years[grid_year_earlier_id]
            grid_year_later <- grid_years[grid_year_later_id]
            
            # Here interpolate between the two grids of grid_year_earlier_id and grid_year_later_id
            grid_earlier <- raster(grid_paths[grid_year_earlier_id])
            grid_later <- raster(grid_paths[grid_year_later_id])
            grid_interpolated <- grid_earlier + (grid_later - grid_earlier) * (year_cur - grid_year_earlier) / (grid_year_later - grid_year_earlier)
            grids_out[[year_cur_id]] <- grid_interpolated
            
          # Else: we are outside the range of the DEM years.
          } else {
            grid_year_closest_id <- which.min(abs(year_dist))
            grids_out[[year_cur_id]] <- raster(grid_paths[grid_year_closest_id])
          }
        }
      }
    }
  }
  
  return(grids_out)
  
}
