###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the routine which calls all the data loading sub-routines.   #
###################################################################################################  

# The <<- puts the variables in the global environment
# (else they would exist only within this function).
func_load_data_all <- function(run_params) {
  
  data_weather               <<-   func_load_weather(run_params)
  data_dems                  <<-   func_load_elevation_grids(run_params, "dem")
  data_dhms                  <<-   func_load_elevation_grids(run_params, "dhm")
  data_surftype              <<-   func_load_surftype_grids(run_params)
  data_outlines              <<-   func_load_outlines(run_params)
  data_radiation             <<-   func_load_radiation_grids(run_params)
  data_massbalance_annual    <<-   func_load_massbalance_measurements(run_params, "annual")
  data_massbalance_winter    <<-   func_load_massbalance_measurements(run_params, "winter")

}
