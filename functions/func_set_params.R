###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the fixed parameter definitions for the model run.           #
###################################################################################################

func_set_params <- function() {

  run_params <- list(
    
    # name_glacier                =    "barkrak",                    # Glacier name, currently unused
    
    #### INPUT-related parameters ####
    # Set paths (don't forget the final / or \).
    dir_data_weather             =   "./input/data/weather/",      # The weather series goes here
    dir_data_dem                 =   "./input/data/dem/",          # Path to the DEM(s)   = elevation grid(s) (masked to the glacier surface, nodata outside)
    dir_data_dhm                 =   "./input/data/dhm/",          # Path to the DHM(s)   = elevation grids(s) (rectangular, to compute slopes and curvatures)
    dir_data_surftype            =   "./input/data/surftype/",     # Grids of surface type (snow/ice/firn/rock/debris) go here
    dir_data_radiation           =   "./input/data/radiation/",    # Grids of potential direct radiation sum go here
    dir_data_massbalance         =   "./input/data/massbalance/",  # The mass balance observations go here
    dir_annual_params            =   "./input/params/",            # The annual model parameter files go here
    
    # Set filenames and input file properties.
    filename_weather             =   "barkrak_barkrak_d.dat",      # File name of the weather series
    file_weather_nskip           =   4,                            # [-]: number of lines to skip in the weather file
    weather_aws_elevation        =   3463,                         # [m a.s.l.]: AWS elevation
    weather_snowfall_temp        =   1.5,                          # [°C]: at this temperature precipitation is half rain, half snow. One degree above it is all rain, one degree below it is all snow (snow fraction is linearly interpolated).
    weather_max_precip_ele       =   3700,                         # [m a.s.l.]: above this altitude, precipitation does not increase any more but becomes constant (cutoff).
    
    grids_crs                    =   CRS(SRS_string  ="EPSG:32642"), # Reference system of the grids, used in slope/aspect computations. Overrides any CRS info reported from the grid files.
    
    filename_dem_prefix          =   "dem_barkrak_",
    filename_dem_suffix          =   ".tif",                       # DEM name is <prefix><year><suffix>
    dem_years                    =   c(2015),                      # Years for which a DEM is available. These whould be sorted in ascending order. Syntax is like this: c(2010, 2015, 2020)
    dem_interpolate              =   FALSE,                        # Should we use linear interpolation to compute each year's DEM? (WARNING: if the glacier area changes, the interpolation between a NA cell and a non-NA one gives NA!)
    
    filename_dhm_prefix          =   "dhm_barkrak_",
    filename_dhm_suffix          =   ".tif",                       # DHM name is <prefix><year><suffix>
    dhm_years                    =   c(2015),                      # Years for which a DHM is available. These whould be sorted in ascending order.
    dhm_interpolate              =   FALSE,                        # Should we use linear interpolation to compute each year's DHM?
    
    filename_surftype_prefix     =   "surf_type_barkrak",
    filename_surftype_suffix     =   ".grid",                      # Surface type filename is <prefix><year><suffix>
    surftype_years               =   c(2015),                      # Years for which a surface type file is available. These whould be sorted in ascending order.
    
    filename_radiation_prefix    =   "dir",
    filename_radiation_suffix    =   "24.grid",                    # Radiation files are called <prefix><doy><suffix> where <doy> is the day of year, zero-padded to length 3 (e.g. 001).

    filename_massbalance_annual  =   "peg_barkrak_orig.dat",       # File name of the annual mass balance observations
    filename_massbalance_winter  =   "peg_barkrak_w.dat",          # File name of the winter mass balance observations
    
    filename_params_prefix       =   "param_",
    filename_params_suffix       =   ".dat",                       # Annual parameters filename is <prefix><year><suffix>
    
    
    #### TOPOGRAPHICAL SNOW DISTRIBUTION-related parameters ####
    curvature_dhm_smooth         =   1.0,                          # [cells]: amount of gaussian smoothing applied before computing curvature (which is very sensitive to DEM noise, unlike slope). Can be non-integer. 1.0 is good for a normal 20 m DEM.
    curvature_cutoff_fact        =   1.2,                          # [-]: multiplier for the curvature cutoff threshold at which the snow distribution is not further changed. The threshold is given by the smaller of the two curvature extremes (positive and negative) divided by this factor. Only values >  = 1 make sense.
    curvature_effect_limit       =   0.5,                          # [-]: maximum effect of curvature, i.e. the curvature multiplier will be within [1 ± curvature_effect_limit]. Only values between 0 and 1 make sense.
    
    elevation_effect_threshold   =   3700,                         # [m]: elevation above which snow accumulation decreases (wind effect)
    elevation_effect_fact        =   1.0,                          # [-]: strength of snow accumulation decrease at very high altitude. Only values between 0 and 1 make sense. At 0 accumulation does not decrease, at 1 accumulation decreases to 0 at the highest point in the DHM.
    
    
    #### AVALANCHE model parameters ####
    elevation_equal_threshold    =   1e-3,                         # [m]: threshold for considering two elevation values equal when we look for problematic flat patches
    deposition_slope_lim         =   40,                           # [°]: at or above this slope value snow will not be deposited during an avalanche. A lower value makes avalanches travel farther. Called beta_lim in Gruber (2007).
    deposition_mass_lim          =   4000,                         # [kg m-2]: maximum deposition during an avalanche. A lower value makes avalanches travel farther. Called D_lim in Gruber (2007).
    movable_slope_lim_lower      =   30,                           # [°]: above this slope value, there is a linearly increasing movable fraction in the initial mass distribution, for avalanches. A lower value makes avalanches start also on more gentle slopes.
    movable_slope_lim_upper      =   60,                           # [°]: above this slope value, all input snow is movable in the avalanche routine.
    deposition_max_ratio_init    =   12,                           # [-]: ONLY for the initial snow distribution grid, how much accumulation can locally result from an avalanche relative to the mean snow distribution before the avalanche? This controls how far avalanches travel, it should be set to a value low enough that avalanches don't bring snow below the marked snow line elevation, and high enough that avalanche deposits look plausible. An exploratory value of 10 can make sense.
    
    
    #### WINTER SNOW PROBES interpolation parameters ####
    snow_probes_idw_exp          =   0.75,                         # [-]: exponent for the IDW interpolation of winter snow measurements
    
    
    #### INITIAL SNOW COVER parameters ####
    initial_snowline_elevation   =   3700,                         # [m]: initial snow line elevation, at the beginning of each simulated year. In the future it will be customizable for each year, and it will be possible to use as initial snow cover the result of the previous year (so that this elevation is only used for the first modeled year).
    initial_snow_gradient        =   200,                          # [mm w.e. (100 m)-1]: increase of the initial snow amount for every 100 m elevation above the snow line
    initial_snow_dist_red_fac    =   0.5,                          # [-]: reduction factor to decrease the importance of the snow distribution variability (all components except winter snow probes), for the computed initial snow cover (of each year). 0 means uniform snow distribution, 1 means no reduction.
    
    
    #### MELT MODEL fixed parameters ####
    debris_red_fac               =   0.6,                          # [-]: reduction factor of melt over debris-covered ice.
    
    
    #### TIME-related parameters ####
    first_year                   =   2017,                         # First modeled year (usually from October of the previous year to September of this year)
    last_year                    =   2020                          # Last modeled year (same as previous comment)
    
  )
  
  
  #### DERIVED parameters, automatically computed ####
  run_params$years <- run_params$first_year:run_params$last_year
  run_params$n_years <- length(run_params$years)
  
  run_params$curvature_dhm_smooth <- max(1e-9,run_params$curvature_dhm_smooth) # The gaussian smoothing fails if sigma   = 0 (but 1e-9 still corresponds to no smoothing!)
  run_params$dhm_smooth_windowsize <- max(5, 2 * run_params$curvature_dhm_smooth + 1)
  
  return(run_params)

}
