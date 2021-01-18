###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the routine to run the mass balance model over a period      #
#                 (winter or year).                                                               #
###################################################################################################  


# year_cur_params are the melt/accumulation model parameters which we will optimize!
# (For the first iteration we take the ones loaded from file).
func_massbal_model <- function(run_params,
                               year_cur_params,
                               model_time_bounds,
                               dhm_values,
                               surftype_init_values,
                               snowdist_init_values,
                               radiation_grids,
                               weather_series) {
  
  year_cur_params$rad_fact_firn <- (year_cur_params$rad_fact_ice + year_cur_params$rad_fact_snow) / 2
  
  # TESTING #
  dhm_values <- getValues(data_dhms$elevation[[1]])
  surftype_init_values <- getValues(data_surftype[[1]])
  snowdist_init_values <- getValues(snowdist_init)
  ###########
  
  
  
  #### PROCESS AND EXTRAPOLATE WEATHER SERIES ####
  # Select weather series period.
  # If the weather series time specification is wrong
  # this step is where it all falls apart.
  weather_series_cur <- data_weather[which(data_weather$timestamp == model_time_bounds[1]):(which(data_weather$timestamp == model_time_bounds[2]) - 1),]
  
  # Correct precipitation undercatch with given parameters
  # (lower correction in summer).
  weather_series_cur$precip_corr <- weather_series_cur$precip * (1 + (year_cur_params$prec_corr / 100.))
  ids_summer <- (as.integer(format(weather_series_cur$timestamp, "%m")) %in% 5:9) # Logical indices: TRUE in May to September, FALSE elsewhere.
  weather_series_cur$precip_corr[ids_summer] <- weather_series_cur$precip[ids_summer] * (1 + (year_cur_params$prec_summer_fact * year_cur_params$prec_corr / 100.))
  
  # These two vectors hold the whole gridded temperature and precipitation series.
  # The computation uses the automatic repetition of a vector when it is
  # multiplied element-wise by a longer vector.
  # Temperature in °C, snowfall in mm w.e.
  vec_temperature     <- rep(weather_series_cur$t2m_mean, each = run_params$grid_ncells) + year_cur_params$temp_elegrad * (dhm_values - run_params$weather_aws_elevation) / 100
  vec_solid_prec_frac <- pmax(0, pmin(1, ((1 + run_params$weather_snowfall_temp) - vec_temperature) / 2))
  vec_snowfall        <- vec_solid_prec_frac * rep(weather_series_cur$precip_corr, each = run_params$grid_ncells) * (1 + (pmin(run_params$weather_max_precip_ele, dhm_values) - run_params$weather_aws_elevation) * year_cur_params$prec_elegrad / 1e4 ) # 1e4: gradient is in [% / 100 m], we want [fraction / m].
  
  
  #### CREATE OUTPUT VECTORS ####
  # This should be equal to the length of the selected weather series!
  model_days_n <- as.integer(difftime(model_time_bounds[2], model_time_bounds[1], "days"))
  
  # Length of each modeled vector.
  # NOTE: modeled vectors have one timestep more than weather vectors because
  # they also store the initial timestep.
  # We store each daily map in
  # (nrow * ncol) components of the vector, going sequentially.
  # The first stored map holds the initial conditions (before the first day),
  # the last holds the final conditions (after the last day).
  # The vector element at index (2*nrow*ncol) + (ncol + 1) is then
  # the first (leftmost) cell of the second row, after the second day
  # of melt and accumulation.
  # Remember that C++ indices start at 0 instead of 1!
  vec_items_n <- run_params$grid_ncol * run_params$grid_nrow * (model_days_n + 1)
  
  # These three vectors will hold all the output of the mass balance model.
  vec_snow_swe      <- rep(NA_real_, vec_items_n)
  vec_surf_type     <- rep(NA_real_, vec_items_n) # We use 0 for ice, 1 for firn, 2 for snow, 4 for rock, 5 for debris.
  vec_massbal_cumul <- rep(NA_real_, vec_items_n)
  melt_cur          <- rep(NA_real_, run_params$grid_ncells) # This instead holds only a single timestep. Here goes the daily melt amount.

  # Fill vectors with initial conditions.
  vec_snow_swe[1:run_params$grid_ncells]  <- snowdist_init_values
  vec_surf_type[1:run_params$grid_ncells] <- surftype_init_values
  vec_surf_type[which(vec_snow_swe[1:run_params$grid_ncells] > 0)] <- 2  # Add computed snow to the initial ice/firn/debris map.
  vec_massbal_cumul[1:run_params$grid_ncells] <- 0
  
  
  #### MAIN SIMULATION LOOP ####
  # day_id <- 1
  for (day_id in 1:model_days_n) { # DISABLED FOR TESTING
    cat(day_id, "/", model_days_n, "\n")
    
    doy <- weather_series_cur$doy[day_id]
    
    radiation_cur <- getValues(data_radiation[[doy]])
    
    offset_cur <- day_id * run_params$grid_ncells
    offset_prev <- (day_id - 1) * run_params$grid_ncells
    cells_cur  <- offset_cur + 1:run_params$grid_ncells # Indices of all the grid cells with values at the end of the current day.
    cells_prev <- cells_cur - run_params$grid_ncells    # Indices of all the grid cells with values at the beginning of the current day.
    
    # Set the entire melt_cur to NA before computing,
    # to avoid any possible problems from values of the
    # previous iteration.
    melt_cur[1:run_params$grid_ncells] <- NA_real_
    
    # TODO
    # Compute condition for avalanche
    # Check condition for avalanche, if true run avalanche
    # END TODO
    
    # Melt ice, firn and debris-covered cells.
    # We take the surf type from the previous
    # timestep to find out which cells to consider.
    # The result is a vector of cell indices
    # directly applicable to the melt_cur vector
    # (i.e. indices starting at 1).
    cells_ice    <- which(vec_surf_type[cells_prev] == 0)
    cells_firn   <- which(vec_surf_type[cells_prev] == 1)
    cells_snow   <- which(vec_surf_type[cells_prev] == 2)
    cells_debris <- which(vec_surf_type[cells_prev] == 5)
    
    # Compute melt amounts.
    melt_cur[cells_ice]    <- year_cur_params$melt_factor + 24 * year_cur_params$rad_fact_ice / 1000. * radiation_cur[cells_ice] * vec_temperature[offset_prev + cells_ice] # We use offset_prev in the temperature vector because it has one timestep less than the modeled grids (which also have the initial conditions as first timestep).
    melt_cur[cells_firn]   <- year_cur_params$melt_factor + 24 * year_cur_params$rad_fact_firn / 1000. * radiation_cur[cells_firn] * vec_temperature[offset_prev + cells_firn]
    melt_cur[cells_snow]   <- year_cur_params$melt_factor + 24 * year_cur_params$rad_fact_firn / 1000. * radiation_cur[cells_snow] * vec_temperature[offset_prev + cells_snow]
    melt_cur[cells_debris] <- year_cur_params$melt_factor + 24 * run_params$debris_red_fac * year_cur_params$rad_fact_ice / 1000. * radiation_cur[cells_debris] * vec_temperature[offset_prev + cells_debris]
    melt_cur[is.na(melt_cur)] <- 0.0 # Don't melt rock, but never go into the NAs (we care about the SWE over rock, for avalanches!)
    melt_cur <- pmax(0.0, melt_cur)  # Clamp to positive values: negative PDDs do not add mass.
    
    # We check which cells have had their snow cover depleted
    # at the current time step, to change their surface type.
    # We ignore the "mixed" melting regime arising from a day
    # where snow cover is depleted (radiation factor should in
    # principle be partly snow, partly ice or firn or debris).
    ids_swe_depleted <- which(melt_cur[cells_snow] > vec_snow_swe[cells_prev][cells_snow])
    vec_snow_swe[cells_cur] <- pmax(0, vec_snow_swe[cells_prev] - melt_cur)
    vec_surf_type[cells_cur] <- vec_surf_type[cells_prev]
    vec_surf_type[cells_cur][cells_snow][ids_swe_depleted] <- surftype_init_values[cells_snow][ids_swe_depleted]
    
    # Add accumulation and update cumulative mass balance.
    vec_snow_swe[cells_cur] <- vec_snow_swe[cells_cur] + vec_snowfall[cells_prev]
    vec_massbal_cumul[cells_cur] <- vec_massbal_cumul[cells_prev] - melt_cur + vec_snowfall[cells_prev]
    vec_surf_type[cells_cur][which(vec_snowfall[cells_prev] > 0.0)] <- 2 # Mark surface as snow after snowfall.
    
    
    # Experimental plot of daily evolution.
    ras <- setValues(data_dhms$elevation[[1]], vec_surf_type[cells_cur])
    plot_df <- data.frame(coordinates(ras))
    max_swe <- 1500
    plot_df$swe <- clamp(vec_snow_swe[cells_cur], -Inf, max_swe)
    plot_df$snow <- as.integer(plot_df$swe > 0)
    plot_df$surf <- vec_surf_type[cells_cur]
    date_text <- format(weather_series_cur$timestamp[day_id], "%Y/%m/%d")
    ggplot(plot_df) +
      surf_base +
      geom_raster(aes(x = x, y = y, fill = swe, alpha = as.character(snow))) +
      scale_alpha_manual(values = c("0" = 0, "1" = 1)) +
      annotate("label", x = Inf, y = Inf, hjust = 1.3, vjust = 1.5, label = date_text) +
      scale_fill_distiller(palette = "RdPu", direction = 1, limits = c(0,max_swe)) +
      guides(alpha = "none") +
      theme_void()
    ggsave(paste("output/surftype/", sprintf("%03d", day_id), ".png", sep=""), width = 5, height = 3)

  }
  
}
