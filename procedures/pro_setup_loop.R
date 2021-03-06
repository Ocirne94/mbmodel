###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the code to prepare variables before the main loop.          #
###################################################################################################

# NOTE: this code is source()'d as part of main.R.
# We put code here just to make it more organized.


# Cleanup memory (temporary variables during loading!)
invisible(gc())


# In this data frame we put all the annual results
# for the overview document: glacier-wide mass balance
# (annual and winter, all versions), ELA, AAR, RMSE,
# optimized parameters, cumulative mass balance.
df_overview <- data.frame(year                = run_params$years,
                          mb_annual_meas_corr = NA,
                          mb_annual_meas      = NA,
                          mb_annual_hydro     = NA,
                          mb_annual_fixed     = NA,
                          mb_winter_meas      = NA, # This stays NA unless winter measurements are available.
                          mb_winter_fixed     = NA,
                          ela                 = NA,
                          aar                 = NA,
                          rmse                = NA,
                          melt_factor         = NA,
                          rad_fact_ice        = NA,
                          rad_fact_snow       = NA,
                          prec_corr           = NA,
                          mb_cumul            = NA)

# This will be set to TRUE after the first modeled year,
# to enable re-using of the modeled SWE as starting condition.
swe_prev_available <- FALSE

# Here we will put just the final mass balance for each
# year, to produce the overview_areaplot multi-page PDF file.
overview_areaplots <- list()

# Create output directory.
dir.create(file.path(run_params$output_dirname, "annual_results"), recursive = TRUE)
