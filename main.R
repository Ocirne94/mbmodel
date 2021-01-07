###################################################################################################
# Author:         Enrico Mattea (@unifr.ch), inspired by the IDL version by Matthias Huss.        #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the main loop and instructions.                              #
# Latest update:  2021.1.6                                                                        #
###################################################################################################

## GLOBALS
# list run_params
# data.frame weather_series
# list DEMs (with raster objects, already processed: selected from closest year and possibly interpolated)
# list surface_type_maps (similar to DEMs)
# data.frame mass_balance_readings (directly read from file, optionally more than one file to separate summer and winter readings)

## UPDATED EACH YEAR
# class cur_year_parameters (read from file)
# raster cur_year_init_snowcover (either read from file, or estimated, or taken from previous year, depending on which year it is within the simulation and on the global options)

#### Include libraries ####
library(raster)

#### Load function definitions ####
invisible(sapply(paste("functions/", list.files("functions/"), sep=""), source))

#### Set run parameters ####
run_params       <-   func_set_params()

#### Load input data ####
data_weather     <-   func_load_weather(run_params)
data_dems        <-   func_load_elevation_grids(run_params, "dem")
data_dhms        <-   func_load_elevation_grids(run_params, "dhm")
data_surftype    <-   func_load_surftype_grids(run_params)
data_radiation   <-   func_load_radiation_grids(run_params)
data_massbalance <-   func_load_massbalance_measurements(run_params)

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


