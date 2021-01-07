###################################################################################################
# Author:         Enrico Mattea (@unifr.ch), inspired by the IDL version by Matthias Huss.        #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains only the commands to load all the other  function files.     #
# Latest update:  2021.1.6                                                                        #
################################################################################################### 


# When this file is read by main.R, we are in the main directory (not the functions/ one).
# So we have to specify the path here.
files.sources = list.files()
sapply(paste("functions/", list.files("functions/"), sep=""), source)
# source("functions/func_set_params.R")
# source("functions/func_load_weather.R")
# source("functions/func_load_elevation_grids.R")
# source("functions/func_load_surftype_grids.R")
# source("functions/func_load_radiation_grids.R")
# source("functions/func_load_massbalance_measurements.R")
