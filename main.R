###################################################################################################
# Author:         Enrico Mattea (@unifr.ch), inspired by the IDL version by Matthias Huss.        #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the main loop and instructions.                              #
# Latest update:  2021.1.5                                                                        #
###################################################################################################

## GLOBALS
# class generalParameters
# data.frame weatherSeries
# list DEMs (with raster objects, already processed: selected from closest year and possibly interpolated)
# list surfaceTypeMaps (similar to DEMs)
# data.frame massBalanceReadings (directly read from file, optionally more than one file to separate summer and winter readings)

## UPDATED EACH YEAR
# class curYearParameters (read from file)
# raster curYearInitSnowCover (either read from file, or estimated, or taken from previous year, depending on which year it is within the simulation and on the global options)
