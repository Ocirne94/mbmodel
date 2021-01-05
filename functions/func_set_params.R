###################################################################################################
# Author:         Enrico Mattea (@unifr.ch), inspired by the IDL version by Matthias Huss.        #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the fixed parameter definitions for the model run.           #
# Latest update:  2021.1.5                                                                        #
###################################################################################################

run_params <- list(
  
  # Set paths.
  dir_data_dem =             "./data/dem/",
  dir_data_surf_type =       "./data/surf_type/",
  dir_data_radiation =       "./data/radiation/",
  dir_data_weather =         "./data/weather/",
  dir_data_mass_balance =    "./data/mass_balance/",
  
  
  
  first_year =               2017,                         # First modeled year (usually from October of the previous year to September of this year)
  last_year =                2020                          # Last modeled year (same as previous comment)
                            
)


