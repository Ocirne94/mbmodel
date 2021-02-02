###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the fixed parameter definitions for the model run.           #
###################################################################################################

dir_data_base                <-   "./input/data_gries/"                   # The base directory for all the data

func_set_params <- function() {

  run_params <- list(
    
    name_glacier                =    "gries",                    # Glacier name, used for output directory naming.
    
    #### INPUT-related parameters ####
    # Set paths (don't forget the final / or \).
    dir_data_weather             =   paste0(dir_data_base, "weather/"),      # The weather series goes here
    dir_data_dem                 =   paste0(dir_data_base, "dem/"),          # Path to the DEM(s) = elevation grid(s) (masked to the glacier surface, nodata outside)
    dir_data_dhm                 =   paste0(dir_data_base, "dhm/"),          # Path to the DHM(s) = elevation grids(s) (rectangular, to compute slopes and curvatures)
    dir_data_surftype            =   paste0(dir_data_base, "surftype/"),     # Path to the grids of surface type (snow/ice/firn/rock/debris) go here
    dir_data_outline             =   paste0(dir_data_base, "outline/"),      # Path to the outlines
    dir_data_radiation           =   paste0(dir_data_base, "radiation/"),    # Path to the grids of potential direct radiation (daily sums)
    dir_data_massbalance         =   paste0(dir_data_base, "massbalance/"),  # The mass balance observations go here
    dir_annual_params            =   "./input/params_gries/",            # The annual model parameter files go here
    
    # Set filenames and input file properties.
    filename_weather             =   "gries_climap_d.dat",      # File name of the weather series
    file_weather_nskip           =   2,                            # [-]: number of lines to skip in the weather file
    weather_aws_elevation        =   2800,                         # [m a.s.l.]: AWS elevation
    weather_snowfall_temp        =   1.5,                          # [°C]: at this temperature precipitation is half rain, half snow. One degree above it is all rain, one degree below it is all snow (snow fraction is linearly interpolated).
    weather_max_precip_ele       =   4000,                         # [m a.s.l.]: above this altitude, precipitation does not increase any more but becomes constant (cutoff).
    
    grids_crs                    =   CRS(SRS_string  ="EPSG:21781"), # Reference system of the grids, used in slope/aspect computations. Overrides any CRS info reported from the grid files.
    
    filename_dem_prefix          =   "gl_gries",
    filename_dem_suffix          =   ".grid",                       # DEM name is <prefix><year><suffix>
    dem_years                    =   c(2012,2015),                      # Years for which a DEM is available. These should be sorted in ascending order. Syntax is like this: c(2010, 2015, 2020)
    dem_interpolate              =   FALSE,                        # [TRUE/FALSE]: should we use linear interpolation to compute each year's DEM? (WARNING: if the glacier area changes, the interpolation between a NA cell and a non-NA one gives NA!)
    
    filename_dhm_prefix          =   "dhm_gries",
    filename_dhm_suffix          =   ".grid",                       # DHM name is <prefix><year><suffix>
    dhm_years                    =   c(2012,2015),                      # Years for which a DHM is available. These should be sorted in ascending order.
    dhm_interpolate              =   FALSE,                        # [TRUE/FALSE]: should we use linear interpolation to compute each year's DHM?
    
    filename_surftype_prefix     =   "firn_gries",
    filename_surftype_suffix     =   ".grid",                      # Surface type filename is <prefix><year><suffix>
    surftype_years               =   c(2012,2015),                      # Years for which a surface type file is available. These whould be sorted in ascending order.
    
    filename_outline_prefix      =   "gries",
    filename_outline_suffix      =   "_gltot.xyzn",                # Outline name is <prefix><year><suffix>
    outline_years                =   c(2012,2015),                      # Years for which an outline is available.
    
    filename_radiation_prefix    =   "dir",
    filename_radiation_suffix    =   "24.grid",                    # Radiation files are called <prefix><doy><suffix> where <doy> is the day of year, zero-padded to length 3 (e.g. 001).

    filename_massbalance_annual  =   "peg_gries.dat",       # File name of the annual mass balance observations
    filename_massbalance_winter  =   "peg_gries_w.dat",          # File name of the winter mass balance observations
    
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
    deposition_mass_lim          =   2000,                         # [kg m-2]: maximum deposition during an avalanche. A lower value makes avalanches travel farther. Called D_lim in Gruber (2007).
    movable_slope_lim_lower      =   30,                           # [°]: above this slope value, there is a linearly increasing movable fraction in the initial mass distribution, for avalanches. A lower value makes avalanches start also on more gentle slopes.
    movable_slope_lim_upper      =   60,                           # [°]: above this slope value, all input snow is movable in the avalanche routine.
    deposition_max_ratio_init    =   12,                           # [-]: ONLY for the initial snow distribution grid, how much accumulation can locally result from an avalanche relative to the mean snow distribution before the avalanche? This controls how far avalanches travel, it should be set to a value low enough that avalanches don't bring snow below the marked snow line elevation, and high enough that avalanche deposits look plausible. An exploratory value of 10 can make sense.
    
    
    #### WINTER SNOW PROBES interpolation parameters ####
    snow_probes_idw_exp          =   0.75,                         # [-]: exponent for the IDW interpolation of winter snow measurements
    
    
    #### INITIAL SNOW COVER parameters ####
    initial_snowline_elevation   =   2800,                         # [m]: initial snow line elevation, at the beginning of each simulated year. In the future it will be customizable for each year, and it will be possible to use as initial snow cover the result of the previous year (so that this elevation is only used for the first modeled year).
    initial_snow_gradient        =   50,                           # [mm w.e. (100 m)-1]: increase of the initial snow amount for every 100 m elevation above the snow line
    initial_snow_dist_red_fac    =   0.5,                          # [-]: reduction factor to decrease the importance of the snow distribution variability (all components except winter snow probes), for the computed initial snow cover (of each year). 0 means uniform snow distribution, 1 means no reduction.
    initial_snow_dist_from_model =   TRUE,                         # [TRUE/FALSE]: if TRUE, use the simulated SWE of the previous year as starting condition for the simulation. If FALSE, compute initial SWE from topography and given parameters. The first simulated year always uses a computed initial SWE since there is no previous modeled year.
    
    
    #### ACCUMULATION and MELT MODEL fixed parameters ####
    debris_red_fac               =   0.6,                          # [-]: reduction factor of melt over debris-covered ice.
    accum_probes_red_fac         =   0.5,                          # [-]: reduction factor to decrease the importance of the snow probes distribution when distributing snowfall over the grid, in case those are measured also over avalanche deposits (else we would be accounting twice for avalanche redistribution, since we run a process-based avalanche model). 0 means uniform distribution, 1 means no redution in variability.
    accum_snow_dist_red_fac      =   0.5,                          # [-]: reduction factor to decrease the importance of the topographic snow distribution variability (curvature and elevation cutoff) when distributing snowfall over the grid. 0 means uniform snow distribution, 1 means no reduction.
    model_avalanche_dates        =   c("1/31", "4/30", "6/30", "7/31", "8/31"),  # [month/day]: dates at which an avalanche is simulated. Usually one at the end of winter (but before winter stakes are measured), and one or more in summer to avoid overloading the slopes with summer snowfall.
    
    albedo_ice_decrease_elev     =   2900.,                           # [m]: below this altitude, the ice albedo decreases linearly with altitude (darker ice).
    albedo_ice_decrease_fact     =   0.0114,                       # [fraction m-1]: rate of increase above 1 (with decreasing altitude) of the ice albedo factor (multiplying ice melt).

    
    #### STAKES COMPARISON parameters ####
    stakes_unknown_latest_start  =   "2/28",                       # [month/day]: in the automatic search of the start date for snow pits and depth probings without a measured start date, we search no later than this day of year. The starting date will be set to the day of the minimum cumulative mass balance between the start of the simulation and the date set here. Something like end of February should be safe for all stakes. 
    
    
    #### MODEL OPTIMIZATION parameters ####
    optim_max_corr_fact          =   1,                            # [-]: maximum allowable positive correction to the melt factor and the radiation factor during optimization, in units of the factors themselves (i.e. by how many times these can be increased). Only positive values make sense. A larger value is safer if a reasonable value for the melt factors is not known, but the optimization will be a bit slower. There is no parameter for the negative correction: it is automatically set to maximum 0.
    optim_bias_threshold         =   1,                            # [mm w.e.]: if abs(bias) is below this threshold then we stop the optimization. This saves us a couple iterations since the optim() function will stop when the value *change* is less than a threshold, not the value itself.
    
    
    #### FIXED MASS BALANCE PERIODS choice ####
    massbal_fixed_annual_start   =   "10/1",                       # [month/day]: start of the user-defined fixed period for annual mass balance evaluation. This is referred to (<year_cur> - 1).
    massbal_fixed_annual_end     =   "8/12",                       # [month/day]: end of the user-defined fixed period for annual mass balance evaluation. This is referred to <year_cur>.
    massbal_fixed_winter_start   =   "10/1",                       # [month/day]: start of the user-defined fixed period for winter mass balance evaluation. This is referred to (<year_cur> - 1).
    massbal_fixed_winter_end     =   "4/30",                       # [month/day]: end of the user-defined fixed period for winter mass balance evaluation. This is referred to <year_cur>.
    
    
    #### MASS BALANCE PROCESSING choices ####
    mb_optimization_skip         =   TRUE,                         # [TRUE/FALSE]: CURRENTLY NOT IMPLEMENTED, SHOULD WE SKIP THE OPTIMIZATION OF THE MASS BALANCE MODEL?
    mb_corr_bands_skip           =   FALSE,                        # [TRUE/FALSE]: CURRENTLY NOT IMPLEMENTED, SHOULD WE SKIP THE CORRECTION OF MASS BALANCE BASED ON ELEVATION BANDS?
    
    
    #### MODELED YEARS choice ####
    first_year                   =   2019,                         # First modeled year (usually from October of the previous year to September of this year)
    last_year                    =   2020                          # Last modeled year (same as previous comment)
    
  )
  
  
  #### DERIVED parameters, automatically computed ####
  run_params$years                       <- run_params$first_year:run_params$last_year
  run_params$n_years                     <- length(run_params$years)
  
  run_params$curvature_dhm_smooth        <- max(1e-9,run_params$curvature_dhm_smooth) # The gaussian smoothing fails if sigma   = 0 (but 1e-9 still corresponds to no smoothing!)
  run_params$dhm_smooth_windowsize       <- max(5, 2 * run_params$curvature_dhm_smooth + 1)
  
  run_params$model_avalanche_dates       <- format(as.Date(run_params$model_avalanche_dates, format = "%m/%d"), format = "%m/%d") # Add leading zeroes to single-digit values if needed.
  
  run_params$stakes_unknown_latest_start <- format(as.Date(run_params$stakes_unknown_latest_start, format = "%m/%d"), format = "%m/%d") # Same.
  
  run_params$massbal_fixed_annual_start <- format(as.Date(run_params$massbal_fixed_annual_start, format = "%m/%d"), format = "%m/%d")
  run_params$massbal_fixed_annual_end <- format(as.Date(run_params$massbal_fixed_annual_end, format = "%m/%d"), format = "%m/%d")
  run_params$massbal_fixed_winter_start <- format(as.Date(run_params$massbal_fixed_winter_start, format = "%m/%d"), format = "%m/%d")
  run_params$massbal_fixed_winter_end <- format(as.Date(run_params$massbal_fixed_winter_end, format = "%m/%d"), format = "%m/%d")
  
  return(run_params)

}
