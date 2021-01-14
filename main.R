###################################################################################################
# Author:         Enrico Mattea (@unifr.ch), inspired by the IDL version by Matthias Huss.        #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the main loop and instructions.                              #
###################################################################################################

#### Include libraries ####
library(raster)
library(spatialEco)   # curvature()
library(scales)       # rescale()
library(Rcpp)         # avalanche function implemented in C++ for performance
library(gstat)

#### Load function definitions ####
invisible(sapply(paste("functions/", list.files("functions/", pattern = "\\.R$"), sep=""), source))
sourceCpp("functions/func_avalanche_gruber.cpp")

#### Set run parameters ####
run_params       <-   func_set_params()

#### Load input data ####
if (run_params$boot_file_read) {
  
  load(run_params$boot_file_name)
  
} else {
  
  data_weather               <-   func_load_weather(run_params)
  data_dems                  <-   func_load_elevation_grids(run_params, "dem")
  data_dhms                  <-   func_load_elevation_grids(run_params, "dhm")
  data_surftype              <-   func_load_surftype_grids(run_params)
  data_radiation             <-   func_load_radiation_grids(run_params)
  data_massbalance_annual    <-   func_load_massbalance_measurements(run_params, "annual")
  data_massbalance_winter    <-   func_load_massbalance_measurements(run_params, "winter")
  
  grids_avalanche            <-   func_compute_avalanche_fixed_grids(run_params, data_dhms)
  grids_snowdist_topographic <-   func_compute_snowdist_topographic(run_params, data_dhms, data_dems)
  
  if (run_params$boot_file_write) {
    save(list = c(apropos("^data_"), apropos("^grids_")), file = run_params$boot_file_name)
  }
}

# Cleanup memory (temporary variables during loading!)
gc()

#### Main loop ####
for (year_id in 1:run_params$n_years) {
  
  year_cur <- run_params$years[year_id]
  
  year_cur_params <- func_load_year_params(run_params, year_cur)
  
  massbal_annual_ids <- func_select_year_measurements(data_massbalance_annual, year_cur)
  massbal_winter_ids <- func_select_year_measurements(data_massbalance_winter, year_cur)
  
  
  # For the *computed* initial snow cover at the beginning of the modeling period:
    # parameter "initial_snowline_elevation" (in the IDL model it is a
      # constant, it could make sense to allow an annual override)
    # initial snow cover is 0 below the snowline elevation, above it is
      # given by a linear gradient with elevation (parameter!) multiplied by the relative distribution grid,
      # which is computed as normalize_over_glacier(reduce_effect(elevation_effect * curvature_effect * avalanche_effect) * smooth(IDW(winter_stakes_for_the_year)))
    # we can compute it as normalize_over_glacier(reduce_effect(avalanche(elevation_effect * curvature_effect)) * smooth(IDW(winter_stakes_for_the_year)))
    
  
}


