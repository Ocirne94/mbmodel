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

#### Load function definitions ####
invisible(sapply(paste("functions/", list.files("functions/", pattern = "\\.R$"), sep=""), source))
sourceCpp("functions/func_avalanche_gruber.cpp")

#### Set run parameters ####
run_params       <-   func_set_params()

#### Load input data ####
if (run_params$boot_file_read) {
  
  load(run_params$boot_file_name)
  
} else {
  
  data_weather            <-   func_load_weather(run_params)
  data_dems               <-   func_load_elevation_grids(run_params, "dem")
  data_dhms               <-   func_load_elevation_grids(run_params, "dhm")
  data_avalanche          <-   func_compute_avalanche_fixed_grids(run_params, data_dhms)
  data_surftype           <-   func_load_surftype_grids(run_params)
  data_radiation          <-   func_load_radiation_grids(run_params)
  data_massbalance_annual <-   func_load_massbalance_measurements(run_params, "annual")
  data_massbalance_winter <-   func_load_massbalance_measurements(run_params, "winter")
  
  if (run_params$boot_file_write) {
    save(list = apropos("^data_"), file = run_params$boot_file_name)
  }
}

# Cleanup memory (temporary variables during loading!)
gc()

#### Main loop ####
for (year_id in 1:run_params$n_years) {
  
  year_cur <- run_params$years[year_id]
  
  year_cur_params <- func_load_year_params(run_params, year_cur)
  
}

# TODO: continue main loop, with:
  # determination of year boundaries for modeling (should cover measurement period and hydrological year)
  # creation of initial snow cover (load from file, or use result from previous modeling, or estimate from parameters)
  # model run over the year with initial parameters
  # if told to optimize: repeat run, minimizing BIAS and RMS


