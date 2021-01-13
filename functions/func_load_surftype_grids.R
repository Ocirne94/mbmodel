###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the loading routine for the grid(s) of surface type.         #
#                 As output we get a list with one raster per year (selected from the             #
#                 closest input available).                                                       #
#                 NOTE: we should never modify the list elements, rather work on copies.          #
###################################################################################################


func_load_surftype_grids <- function(run_params) {
  
  # Here we will put the output.
  grids_out <- list()
  
  grid_paths <- paste(run_params$dir_data_surftype,
                      run_params$filename_surftype_prefix,
                      run_params$surftype_years,
                      run_params$filename_surftype_suffix,
                      sep = "")
  
  grid_years <- run_params$surftype_years
  
  # For each modeled year find the closest grid year and use its grid.
  for (year_cur_id in 1:run_params$n_years) {
    
    year_cur <- run_params$years[year_cur_id]
    grid_year_closest_id <- which.min(abs(grid_years - year_cur))
    grids_out[[year_cur_id]] <- readAll(raster(grid_paths[grid_year_closest_id]))
    crs(grids_out[[year_cur_id]]) <- run_params$grids_crs
    
  }
  
  return(grids_out)
  
}
