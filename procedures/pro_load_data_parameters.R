###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the code to load input data and parameters.                  #
###################################################################################################

# NOTE: this code is source()'d as part of main.R.
# We put code here just to make it more organized.

#### Load or set fixed run parameters ####
if (params_file_read) {
  load(params_file_name)
} else {
  run_params <- func_set_params()
}

#### Load input data from sources or reboot file ####
if (boot_file_read) {
  load(boot_file_name)
} else {
  func_load_data_all(run_params) # This calls all the data loading routines.
}
