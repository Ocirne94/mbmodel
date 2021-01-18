###################################################################################################
# Author:         Enrico Mattea (@unifr.ch), inspired by the IDL version by Matthias Huss.        #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the main loop and instructions.                              #
###################################################################################################

#### Load from files or reboot file ####
params_file_write <- FALSE                # Save .RData file with run parameters, for faster reload. The file is saved AFTER data loading since some parameters depend on the loaded grids (e.g. cell size).
params_file_read  <- FALSE                # Load .RData file with run parameters, instead of setting new run parametrs.
params_file_name  <- "params_file.RData"  # Name of the .RData run parameters file.

boot_file_write   <- FALSE                # Save .RData file with the input data, for faster reload.
boot_file_read    <- TRUE                 # Load .RData file with the input data, instead of loading input files.
boot_file_name    <- "boot_file.RData"    # Name of the .RData input data file.


#### Include libraries ####
library(raster)
library(spatialEco)   # curvature()
library(scales)       # rescale()
library(Rcpp)         # avalanche function implemented in C++ for performance
library(gstat)        # IDW of snow probing data

# Experimental plots.
library(ggplot2)
library(RStoolbox)


# library(bitops)  # cksum(), for debugging.



#### Load function definitions ####
invisible(sapply(paste("functions/", list.files("functions/", pattern = "\\.R$"), sep=""), source))
sourceCpp("functions/func_avalanche_gruber.cpp", cacheDir = "functions/") # Remove cacheDir option to force reload of the C++ code (useful after changing computer or editing the source file).


#### Load or set run parameters ####
if (params_file_read) {
  load(params_file_name) } else {
  run_params <- func_set_params()
}


#### Load input data from sources or reboot file ####
if (boot_file_read) {
  load(boot_file_name)
} else {
  data_weather               <-   func_load_weather(run_params)
  data_dems                  <-   func_load_elevation_grids(run_params, "dem")
  data_dhms                  <-   func_load_elevation_grids(run_params, "dhm")
  data_surftype              <-   func_load_surftype_grids(run_params)
  data_radiation             <-   func_load_radiation_grids(run_params)
  data_massbalance_annual    <-   func_load_massbalance_measurements(run_params, "annual")
  data_massbalance_winter    <-   func_load_massbalance_measurements(run_params, "winter")
}

# Assign global grid parameters to run_params.
# We do it explicitly here since we need to have
# loaded the elevation grids to get their parameters.
# These overwrite the loaded run parameters, to avoid
# the possibility of dangerous inconsistencies with the loaded data.
run_params$grid_nrow       <- nrow(data_dhms$elevation[[1]])
run_params$grid_ncol       <- ncol(data_dhms$elevation[[1]])
run_params$grid_cell_size  <- xres(data_dhms$elevation[[1]])
run_params$grid_ncells     <- run_params$grid_nrow * run_params$grid_ncol

# Compute fixed grids to speed up processing, but
# not if these were loaded from the boot file.
if (!boot_file_read) {
  grids_avalanche            <-   func_compute_avalanche_fixed_grids(run_params, data_dhms)
  grids_snowdist_topographic <-   func_compute_snowdist_topographic(run_params, data_dhms, data_dems)
}

# Save reboot file.
if (boot_file_write) {
  save(list = c(apropos("^data_"), apropos("^grids_")), file = boot_file_name)
}
# Save parameters file.
if (params_file_write) {
  save(list = "run_params", file = params_file_name)
}

# Cleanup memory (temporary variables during loading!)
gc()



#### TESTING, for experimental plots ####
surf_r <- subs(data_surftype[[1]], data.frame(from = c(0, 1, 4, 5), to = c(100, 170, 0, 70)))
surf_g <- subs(data_surftype[[1]], data.frame(from = c(0, 1, 4, 5), to = c(150, 213, 0, 20)))
surf_b <- subs(data_surftype[[1]], data.frame(from = c(0, 1, 4, 5), to = c(200, 255, 0, 20)))
surf_base <- ggRGB(stack(surf_r, surf_g, surf_b), r = 1, g = 2, b = 3, ggLayer = TRUE)


#### Main loop ####
# for (year_id in 1:run_params$n_years) {
  year_id <- 1 # TESTING!!

  # Select current year.
  year_cur <- run_params$years[year_id]
  year_cur_params <- func_load_year_params(run_params, year_cur)
  
  # Select grids of the current year.
  grid_id <- data_dhms$grid_year_id[year_id]
  
  # Select mass balance measurements of the current year.
  massbal_annual_ids <- func_select_year_measurements(data_massbalance_annual, year_cur)
  massbal_winter_ids <- func_select_year_measurements(data_massbalance_winter, year_cur)
  massbal_annual_cur <- data_massbalance_annual[massbal_annual_ids,]
  massbal_winter_cur <- data_massbalance_winter[massbal_winter_ids,] # Empty if we have no winter stakes for the year.
  
  # Should we make a winter run to optimize the precipitation correction?
  # Only if we have some measurements of winter snow cover, else we can't.
  process_winter <- length(massbal_winter_ids) > 0
  
  # Compute an initial snow distribution for the current year.
  snowdist_init <- func_compute_initial_snow_cover(run_params,
                                                   data_dhms,
                                                   data_dems,
                                                   grids_snowdist_topographic,
                                                   grids_avalanche,
                                                   grid_id,
                                                   massbal_winter_cur)
  
  
  # Four POSIXct time objects: start and end of the annual
  # modeling period, and same for the winter modeling period.
  model_time_bounds   <- func_compute_modeling_periods(run_params, massbal_annual_cur, massbal_winter_cur, year_cur)
  

  if (process_winter)  {
    model_winter_bounds <- model_time_bounds[3:4]
    # TODO:
      # find grid cells of winter stakes
      # run massbal_model (we cannot run it only on the stake grid cells due to avalanches!!)
  }
    
  # Here instead do the annual processing.
  # Find grid cells corresponding to the annual stakes.
  annual_stakes_cells <- fourCellsFromXY(data_dhms$elevation[[grid_id]], as.matrix(massbal_annual_cur[,4:5]))
  
  # Time specification of the annual run.
  model_annual_bounds <- model_time_bounds[1:2]
  

  

  
# }
