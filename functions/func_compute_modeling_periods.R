###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the routine to determine the modeling periods, both for      #
#                 the winter period and the annual model.                                         #
################################################################################################### 


# ANNUAL modeling period starts at the beginning of the observation
# period with the earliest start (i.e. the annual stake which was surveyed first
# on the previous year), but no later than Oct 1 (to include the whole hydrological year).
# Annual modeling period ends at the end of the observation period with the latest end
# (i.e. the annual stake which was surveyed last on the current year), but no earlier than Sep 30 (hydro year).
# WINTER modeling period starts at the beginning of the observation period with the earliest start
# (among the winter ones) and ends at the end of the period with the latest end (among the winter ones).
func_compute_modeling_periods <- function(run_params, massbal_annual, massbal_winter, year_cur) {
  
  hydro_start <- as.POSIXct(paste(year_cur-1, 10, 1), format="%Y %m %d")
  hydro_end   <- as.POSIXct(paste(year_cur, 9, 30), format = "%Y %m %d")
  
  # na.rm because we support NA as start date, meaning
  # "end of previous ablation season" i.e. mass balance minimum.
  # NOTE: we don't change the modeling period to include stakes
  # which start at NA. So if the mass balance minimum happens
  # before Sep 30 of the previous year and there isn't a measurement
  # before that date, the starting date of those stakes will be
  # set at Sep 30.
  # In the future we may support having a custom starting date.
  annual_start <- min(c(hydro_start, massbal_annual$start_date), na.rm = T)
  annual_end   <- max(c(hydro_end, massbal_annual$end_date))
  
  winter_start <- NA
  winter_end   <- NA
  if (length(massbal_winter[,1]) > 0) {
    winter_start <- min(massbal_winter$start_date)
    winter_end   <- max(massbal_winter$end_date)
  }
  
  return(c(annual_start, annual_end, winter_start, winter_end))
  
}
