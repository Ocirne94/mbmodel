###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the fixed parameter definitions for the model run.           #
# Latest update:  2021.1.6                                                                        #
###################################################################################################

func_set_params <- function() {

  run_params <- list(
    
    # name_glacier              =    "barkrak",                    # Glacier name, currently unused
    
    #### INPUT-related parameters ####
    # Set paths (don't forget the final / or \).
    dir_data_weather          =    "./data/weather/",            # The weather series goes here
    dir_data_dem              =    "./data/dem/",                # Path to the DEM(s) = elevation grid(s) (masked to the glacier surface, nodata outside)
    dir_data_dhm              =    "./data/dhm/",                # Path to the DHM(s) = elevation grids(s) (rectangular, to compute slopes and curvatures)
    dir_data_surftype         =    "./data/surftype/",           # Grids of surface type (snow/ice/firn/rock/debris) go here
    dir_data_radiation        =    "./data/radiation/",          # Grids of potential direct radiation sum go here
    dir_data_massbalance      =    "./data/massbalance/",        # The mass balance observations go here
    
    # Set filenames and input file properties.
    filename_weather          =    "barkrak_barkrak_d.dat",      # File name of the weather series
    file_weather_nskip        =    4,                            # Number of lines to skip in the weather file
    
    filename_dem_prefix       =    "dem_barkrak_",
    filename_dem_suffix       =    ".tif",                       # DEM name is <prefix><year><suffix>
    dem_years                 =    c(2015),                      # Years for which a DEM is available. These whould be sorted in ascending order. Syntax is like this: c(2010, 2015, 2020)
    dem_interpolate           =    FALSE,                        # Should we use linear interpolation to compute each year's DEM? (WARNING: if the glacier area changes, the interpolation between a NA cell and a non-NA one gives NA!)
    
    filename_dhm_prefix       =    "dhm_barkrak_",
    filename_dhm_suffix       =    ".tif",                       # DHM name is <prefix><year><suffix>
    dhm_years                 =    c(2015),                      # Years for which a DHM is available. These whould be sorted in ascending order.
    dhm_interpolate           =    FALSE,                        # Should we use linear interpolation to compute each year's DHM?
    
    filename_surftype_prefix  =    "surf_type_barkrak",
    filename_surftype_suffix  =    ".grid",                      # Surface type filename is <prefix><year><suffix>
    surftype_years            =    c(2015),                      # Years for which a surface type file is available. These whould be sorted in ascending order.
    
    filename_radiation_prefix =    "dir",
    filename_radiation_suffix =    "24.grid",                    # Radiation files are called <prefix><doy><suffix> where <doy> is the day of year, zero-padded to length 3 (e.g. 001).

    filename_massbalance      =    "peg_barkrak.dat",            # File name of the mass balance observations
    
    
    #### TIME-related parameters ####
    first_year                =    2017,                         # First modeled year (usually from October of the previous year to September of this year)
    last_year                 =    2020                          # Last modeled year (same as previous comment)
    
  )
  
  
  #### DERIVED parameters, automatically computed ####
  run_params$years <- run_params$first_year:run_params$last_year
  run_params$n_years <- length(run_params$years)
  
  
  return(run_params)

}
