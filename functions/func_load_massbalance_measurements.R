###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the loading routine for the point mass balance measurements. #
#                 As output we get a data.frame:                                                  #
#                   id start_date end_date x y z dh_cm density                                    #
#                 start_date = NA is interpreted as <end of previous ablation season>, useful     #
#                 for probe/snowpit measurements.                                                 #
# Latest update:  2021.1.6                                                                        #
################################################################################################### 

func_load_massbalance_measurements <- function(run_params) {
  
  massbalance_path <- paste(run_params$dir_data_massbalance,
                            run_params$filename_massbalance,
                            sep = "")
  
  data_massbalance <- read.table(massbalance_path, header = FALSE)
  
  names(data_massbalance) <- c("id", "start_date", "end_date", "x", "y", "z", "dh_cm", "density")
  
  # Convert timestamps to time objects.
  data_massbalance$start_date <- as.POSIXct(data_massbalance$start_date, format = "%d.%m.%Y", tz = "UTC")
  data_massbalance$end_date <- as.POSIXct(data_massbalance$end_date, format = "%d.%m.%Y", tz = "UTC")
  
  return(data_massbalance)
  
}