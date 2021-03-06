###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the routine to load all required modules.                    #
###################################################################################################


func_load_libraries <- function() {
  require(raster)
  require(spatialEco)   # curvature()
  require(scales)       # rescale()
  require(Rcpp)         # avalanche function implemented in C++ for performance
  require(gstat)        # IDW of snow probing data
  require(Rfast)        # rowSort() of the stake cells indices
  require(stats)        # uniroot()
  require(sf)           # geom_sf(), to plot the glacier outline in the maps
  require(metR)         # geom_text_contour()
  require(ggplot2)      # Base plotting require
  require(ggpubr)       # Additional plotting functions
  require(grid)         # Additional plotting functions
  require(cowplot)      # Additional plotting functions
  require(ggpattern)    # Additional plotting functions (install with remotes::install_github("coolbutuseless/ggpattern"))
  require(shadowtext)   # Additional plotting functions
  require(reshape2)     # melt()
  require(RStoolbox)    # For the surface type basemap under the SWE plots
}
