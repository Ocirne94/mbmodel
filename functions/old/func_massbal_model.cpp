//###################################################################################################
//# Author:         Enrico Mattea (@unifr.ch)                                                       #
//# Description:    this program models the distributed mass balance of a glacier at daily          #
//#                 resolution, optimizing model parameters towards the best fit with point         #
//#                 mass balance measurements.                                                      #
//#                 This file contains the routine to run the mass balance model over a period      #
//#                 (winter or year).                                                               #
//###################################################################################################   

#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List func_massbal_model(List run_params,
                        List year_cur_params,
                        NumericVector dhm_values,
                        IntegerVector glacier_cell_ids,
                        NumericVector surftype_init_values,
                        NumericVector snowdist_init_values,
                        List radiation_values_list,
                        DataFrame weather_series_cur,
                        NumericVector snowdist_topographic_values_red,
                        NumericVector snowdist_probes_norm_values_red,
                        List grids_avalanche_cur) {
  
  
  
  double year_cur_params_rad_fact_ice = year_cur_params["rad_fact_ice"];
  double year_cur_params_rad_fact_snow = year_cur_params["rad_fact_snow"];
  double year_cur_params_rad_fact_firn = (year_cur_params_rad_fact_ice + year_cur_params_rad_fact_snow) / 2.;
  double year_cur_params_prec_corr = year_cur_params["prec_corr"];
  double year_cur_params_prec_summer_fact = year_cur_params["prec_summer_fact"];
  double year_cur_params_temp_elegrad = year_cur_params["temp_elegrad"];
  double year_cur_params_prec_elegrad = year_cur_params["prec_elegrad"];
  double year_cur_params_melt_factor = year_cur_params["melt_factor"];

  unsigned int run_params_grid_ncells = run_params["grid_ncells"];
  double run_params_weather_aws_elevation = run_params["weather_aws_elevation"];
  double run_params_weather_snowfall_temp = run_params["weather_snowfall_temp"];
  double run_params_weather_max_precip_ele = run_params["weather_max_precip_ele"];
  unsigned int run_params_grid_ncol = run_params["grid_ncol"];
  unsigned int run_params_grid_nrow = run_params["grid_nrow"];
  double run_params_debris_red_fac = run_params["debris_red_fac"];
  

  
  // --- PROCESS AND EXTRAPOLATE WEATHER SERIES --- //
  NumericVector  weather_series_cur_t2m_mean  = weather_series_cur["t2m_mean"];
  NumericVector  weather_series_cur_precip    = weather_series_cur["precip"];
  DatetimeVector weather_series_cur_timestamp = weather_series_cur["timestamp"];
  IntegerVector  weather_series_cur_doy       = weather_series_cur["doy"];

  // Correct precipitation undercatch with given parameters
  // (lower correction in summer).
  NumericVector weather_series_cur_precip_corr = weather_series_cur_precip * (1. + (year_cur_params_prec_corr / 100.));
    
  int model_days_n = weather_series_cur_timestamp.size();
  
  IntegerVector weather_series_cur_months(model_days_n);
  for (int day_id = 0; day_id < model_days_n; day_id ++) {
    Datetime timestamp_cur = weather_series_cur_timestamp[day_id];
    weather_series_cur_months[day_id] = timestamp_cur.getMonth();
  }
  
  LogicalVector ids_summer_logi = (weather_series_cur_months >= 5) & (weather_series_cur_months <= 9);
  NumericVector weather_series_cur_precip_summer = weather_series_cur_precip[ids_summer_logi];
  NumericVector weather_series_cur_precip_summer_corr = weather_series_cur_precip_summer * (1 + (year_cur_params_prec_summer_fact * year_cur_params_prec_corr / 100.));
  weather_series_cur_precip_corr[ids_summer_logi] = weather_series_cur_precip_summer_corr;

  
  // These two vectors (temperature and accumulation) hold the whole
  // gridded temperature and solid precipitation series.
  // The computation uses the automatic repetition of a vector when it is
  // multiplied element-wise by a longer vector.
  // Temperature in °C, snowfall in mm w.e.
  NumericVector vec_temperature     = rep_each(weather_series_cur_t2m_mean, run_params_grid_ncells) + rep(year_cur_params_temp_elegrad * (dhm_values - run_params_weather_aws_elevation) / 100, model_days_n);
  NumericVector vec_solid_prec_frac = pmax(0, pmin(1, ((1 + run_params_weather_snowfall_temp) - vec_temperature) / 2));
  // Unlike in R, here we have to multiply equal-length vectors, so we use rep() as needed.
  NumericVector vec_accumulation    = rep(snowdist_probes_norm_values_red, model_days_n) * rep(snowdist_topographic_values_red, model_days_n) * vec_solid_prec_frac * rep_each(weather_series_cur_precip_corr, run_params_grid_ncells) * rep(1 + (pmin(run_params_weather_max_precip_ele, dhm_values) - run_params_weather_aws_elevation) * year_cur_params_prec_elegrad / 1e4, model_days_n); // 1e4: gradient is in [% / 100 m], we want [fraction / m].
  
  
  // --- CREATE OUTPUT VECTORS --- //
  // Length of each modeled vector.
  // NOTE: modeled vectors have one timestep more than weather vectors because
  // they also store the initial timestep.
  // We store each daily map in
  // (nrow * ncol) components of the vector, going sequentially.
  // The first stored map holds the initial conditions (before the first day),
  // the last holds the final conditions (after the last day).
  // The vector element at index (2*nrow*ncol) + (ncol + 1) is then
  // the first (leftmost) cell of the second row, after the second day
  // of melt and accumulation.
  // Remember that C++ indices start at 0 instead of 1!
  int vec_items_n = run_params_grid_ncol * run_params_grid_nrow * (model_days_n + 1);
  IntegerVector grid_indices = seq(0,run_params_grid_ncells-1);
    
  // These three vectors will hold all the output of the mass balance model.
  NumericVector vec_snow_swe(vec_items_n);
  NumericVector vec_surf_type(vec_items_n); // We use 0 for ice, 1 for firn, 2 for snow, 4 for rock, 5 for debris.
  NumericVector vec_massbal_cumul(vec_items_n);
  NumericVector melt_cur(run_params_grid_ncells); // This instead holds only a single timestep. Here goes the daily melt amount.
  NumericVector gl_massbal_cumul(model_days_n + 1);  // This holds the cumulative glacier-wide mass balance.
  gl_massbal_cumul[0] = 0.0; // Initial cumulative mass balance is 0.0.
  

  // Fill vectors with initial conditions.
  vec_snow_swe[seq(0,run_params_grid_ncells-1)]  = snowdist_init_values;
  vec_surf_type[seq(0,run_params_grid_ncells-1)] = surftype_init_values;
  // Add computed snow to the initial ice/firn/debris map.
  LogicalVector swe_positive = (as<NumericVector>(vec_snow_swe[grid_indices]) > 0.0); // as<NumericVector> because the [] operator does NOT return a NumericVector already (see https://teuder.github.io/rcpp4everyone_en/201_caution_vector.html#return-type-of-operator).
  NumericVector vec_surf_type_cur_new = vec_surf_type[grid_indices];
  vec_surf_type_cur_new[swe_positive] = 2;
  vec_surf_type[grid_indices] = vec_surf_type_cur_new;
  IntegerVector vec_massbal_cumul_indices = seq(0,run_params_grid_ncells-1);
  vec_massbal_cumul[vec_massbal_cumul_indices] = 0.0; // Unlike in R, we can only assign a vector subset if we have already generated the indices (i.e. "a[seq(1,3)] = 2" doesn't work; a working alternative is a[seq(1,3)] = rep(2,3)).
    

  // This below is a vector[ncells] to keep track of the swe
  // after the last avalanche: we assume that snow can be
  // avalanched only once, so the input grid to the avalanche
  // model will be the difference between the current swe
  // and the latest swe_post_previous_avalanche (clamped to
  // positive values in case of significant melt).
  NumericVector swe_post_previous_avalanche = clone(snowdist_init_values);
  
  
  // --- MAIN SIMULATION LOOP --- //
  // NOTE: we start day_id at 1! So we correct below as needed.
  for (int day_id = 1; day_id <= model_days_n; day_id ++) {
    
    // Rcout << day_id << "/" << model_days_n << std::endl;
    
    // --- SETUP INDICES --- //
    int doy = weather_series_cur_doy[day_id - 1] - 1; // doy is compute in R within [1,366], we want it in [0,365].
    
    NumericVector radiation_cur = radiation_values_list[doy];
    
    // offset_cur is the index (in the output vectors)
    // of the last cell of the previous day.
    // offset_prev is the same but one more day before.
    int offset_cur  = (day_id * run_params_grid_ncells) - 1; // The R version does NOT have -1 (C++ indices start at 0!)
    int offset_prev = ((day_id - 1) * run_params_grid_ncells) - 1; // We cannot use unsigned since at the first iteration this will be -1.
    IntegerVector cells_cur  = offset_cur + seq(1,run_params_grid_ncells); // Indices of all the grid cells with values at the end of the current day.
    IntegerVector cells_prev = cells_cur - run_params_grid_ncells;    // Indices of all the grid cells with values at the beginning of the current day.
    

    // TODO: ADD AVALANCHE ROUTINE!!
    
    
    
    // --- MELT MODEL --- //
    // Set the entire melt_cur to NaN before computing,
    // to avoid any possible problems from values of the
    // previous iteration.
    melt_cur.fill(R_NaN);
    // Melt ice, firn and debris-covered cells.
    // We take the surf type from the previous
    // timestep to find out which cells to consider.
    // The result is a logical vector as long as one grid,
    // so directly applicable to the melt_cur vector
    // and to any single-day subset of our vectors.
    NumericVector surf_type_prev = vec_surf_type[cells_prev];
    LogicalVector cells_ice      = surf_type_prev == 0;
    LogicalVector cells_firn     = surf_type_prev == 1;
    LogicalVector cells_snow     = surf_type_prev == 2;
    LogicalVector cells_debris   = surf_type_prev == 5;
    
    NumericVector temperature_prev = vec_temperature[cells_prev];
    
    
    // Compute melt amounts.
    NumericVector melt_ice     = year_cur_params_melt_factor + 24 * year_cur_params_rad_fact_ice / 1000. * as<NumericVector>(radiation_cur[cells_ice]) * as<NumericVector>(temperature_prev[cells_ice]); // We use offset_prev in the temperature vector because it has one timestep less than the modeled grids (which also have the initial conditions as first timestep).
    melt_cur[cells_ice]        = melt_ice;
    NumericVector melt_firn    = year_cur_params_melt_factor + 24 * year_cur_params_rad_fact_firn / 1000. * as<NumericVector>(radiation_cur[cells_firn]) * as<NumericVector>(temperature_prev[cells_firn]);
    melt_cur[cells_firn]       = melt_firn;
    NumericVector melt_snow    = year_cur_params_melt_factor + 24 * year_cur_params_rad_fact_snow / 1000. * as<NumericVector>(radiation_cur[cells_snow]) * as<NumericVector>(temperature_prev[cells_snow]);
    melt_cur[cells_snow]       = melt_snow;
    NumericVector melt_debris  = year_cur_params_melt_factor + 24 * run_params_debris_red_fac * year_cur_params_rad_fact_ice / 1000. * as<NumericVector>(radiation_cur[cells_debris]) * as<NumericVector>(temperature_prev[cells_debris]);
    melt_cur[cells_debris]     = melt_debris;
    melt_cur[is_nan(melt_cur)] = 0.0; // Don't melt rock, but never go into the NAs (we care about the SWE over rock, for avalanches!)
    melt_cur                   = pmax(0.0, melt_cur);  // Clamp to positive values: negative PDDs do not add mass.
  
    
    // We check which cells have had their snow cover depleted
    // at the current time step, to change their surface type.
    // We ignore the "mixed" melting regime arising from a day
    // where snow cover is depleted (radiation factor should in
    // principle be partly snow, partly ice or firn or debris).
    // Unlike the R version, here we use only indices referred
    // to the whole grid (since multiple inline sub-setting is
    // not available in Rcpp).
    // Like this we set as free from snow also all the cells which
    // were already free, but it is not an issue.
    NumericVector snow_swe_prev = vec_snow_swe[cells_prev];
    LogicalVector swe_depleted_logi = (melt_cur >= snow_swe_prev);
    NumericVector snow_swe_new = pmax(0, as<NumericVector>(vec_snow_swe[cells_prev]) - melt_cur);
    vec_snow_swe[cells_cur] = snow_swe_new;
    vec_surf_type[cells_cur] = vec_surf_type[cells_prev];
    IntegerVector cells_cur_swe_depleted = cells_cur[swe_depleted_logi]; // Indices relative to the whole vec_surf_type vector!
    vec_surf_type[cells_cur_swe_depleted] = surftype_init_values[swe_depleted_logi];
    

    // --- ACCUMULATION and MASS BALANCE --- //
    // Add accumulation and update cumulative mass balance.
    vec_snow_swe[cells_cur] = vec_snow_swe[cells_cur] + vec_accumulation[cells_prev];
    NumericVector massbal_cumul_new = as<NumericVector>(vec_massbal_cumul[cells_prev]) - melt_cur + as<NumericVector>(vec_accumulation[cells_prev]);
    vec_massbal_cumul[cells_cur] = massbal_cumul_new;
    IntegerVector cells_cur_accumulation = cells_cur[as<NumericVector>(vec_accumulation[cells_prev]) > 0.0];
    vec_surf_type[cells_cur_accumulation] = 2; // Mark surface as snow after snowfall.
    gl_massbal_cumul[day_id + 1] = mean(as<NumericVector>(vec_massbal_cumul[offset_cur + glacier_cell_ids]));
    
    
  }
    
  List mb_model_output = List::create(Named("vec_massbal_cumul") = vec_massbal_cumul,
                                      Named("gl_massbal_cumul")  = gl_massbal_cumul);
  
  return mb_model_output;
    
}
