###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the loading routine for the grid(s) of daily potential       #
#                 radiation sum.                                                                  #
#                 As output we get a list with one raster per day of year (366 elements,          #
#                 last two equal)                                                                 #
###################################################################################################


func_load_radiation_grids <- function(run_params) {
  
  # Here we will put the output.
  grids_out <- list()
  
  grid_paths <- paste(run_params$dir_data_radiation,
                      run_params$filename_radiation_prefix,
                      sprintf("%03d", 1:365),
                      run_params$filename_radiation_suffix,
                      sep = "")
  
  # Actual loading happens here.
  for (doy in 1:365) {
    
    cat("\rLoading radiation files...", doy, "/", 365)
    grids_out[[doy]] <- readAll(raster(grid_paths[doy]))
    crs(grids_out[[doy]]) <- run_params$grids_crs
    
  }
  cat("\n")
  
  grids_out[[366]] <- grids_out[[365]]
  
  return(grids_out)
  
}
