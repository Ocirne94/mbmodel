###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the code to load all required modules.                       #
###################################################################################################

# NOTE: this code is source()'d as part of main.R.
# We put code here just to make it more organized.

library(raster)
library(spatialEco)   # curvature()
library(scales)       # rescale()
library(Rcpp)         # avalanche function implemented in C++ for performance
library(gstat)        # IDW of snow probing data
library(Rfast)        # rowSort() of the stake cells indices
library(stats)        # uniroot()
library(sf)           # geom_sf(), to plot the glacier outline in the maps
library(metR)         # geom_text_contour()
library(ggplot2)      # Base plotting library
library(ggpubr)       # Additional plotting functions
library(grid)         # Additional plotting functions
library(cowplot)      # Additional plotting functions
library(ggpattern)    # Additional plotting functions (install with remotes::install_github("coolbutuseless/ggpattern"))
library(shadowtext)   # Additional plotting functions
library(reshape2)     # melt()
library(RStoolbox)    # For the surface type basemap under the SWE plots
