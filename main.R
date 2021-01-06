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
source("functions/func_set_params.R")
source("functions/func_load_weather.R")
source("functions/func_load_elevation_grids.R")



run_params   <-   func_set_params()
data_weather <-   func_load_weather(run_params)
data_dems    <-   func_load_elevation_grids(run_params, "dem")
data_dhms    <-   func_load_elevation_grids(run_params, "dhm")
