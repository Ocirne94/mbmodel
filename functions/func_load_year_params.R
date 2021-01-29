###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the annual parameter loading from file.                      #
#                 For the moment it works with the same file formatting as the IDL model.         #
###################################################################################################

func_load_year_params <- function(run_params, year_cur) {
  
  # Here goes the output.
  year_params <- list()
  
  filepath_params <- paste(run_params$dir_annual_params,
                           run_params$filename_params_prefix,
                           year_cur,
                           run_params$filename_params_suffix,
                           sep = "")
  
  params_raw <- read.delim(filepath_params,
                           header = FALSE,
                           sep = ";",
                           comment.char = "*",
                           stringsAsFactors = FALSE,
                           strip.white = TRUE
                           )
  
  # Assemble output, already converting to numeric types.
  for (param_id in 1:length(params_raw[,1])) {
    year_params[[params_raw[param_id,3]]] <- type.convert(params_raw[param_id,1], as.is = TRUE)
  }
  
  # If present, convert elevation bands for balance correction into a numeric vector.
  if (length(year_params$mb_corr_ele_bands) == 1) {
    year_params$mb_corr_ele_bands <- as.numeric(unlist(strsplit(year_params$mb_corr_ele_bands, ",")))
  }
  
  # Compute ratio of snow to ice radiation factors.
  # We will keep this ratio constant as we optimize
  # the radiation factors, then possibly use it
  # for subsequent optimization steps.
  year_params$rad_fact_ratio_snow_ice <- year_params$rad_fact_snow / year_params$rad_fact_ice
  
  
  # Compute start and end of current hydrological year.
  # The hydrological year starts on 1/10/<Y-1> at 00:00 and ends on 1/10/<Y> at 00:00.
  # Since we use Date objects which don't include the time of day, we can set the
  # hydro end to September 30.
  year_params$hydro_start <- as.Date(paste(year_cur-1, 10, 1), format="%Y %m %d")
  year_params$hydro_end   <- as.Date(paste(year_cur, 9, 30), format = "%Y %m %d")
  
  year_params$fixed_annual_start <- as.Date(paste(year_cur-1, run_params$massbal_fixed_annual_start), format = "%Y %m/%d")
  year_params$fixed_annual_end <- as.Date(paste(year_cur, run_params$massbal_fixed_annual_end), format = "%Y %m/%d")
  
  return(year_params)
  
}
