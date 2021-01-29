###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the routine to produce daily plots of SWE and mass balance.  #
#                 These can be directly turned into a nice animation.                             #
#                 The routine has no return value, it is used as a procedure.                     #
###################################################################################################

func_plot_daily_maps <- function(run_params,
                                 weather_series_cur,
                                 data_dems,
                                 data_outlines,
                                 mb_model_output,
                                 surf_base,
                                 elevation_grid_id,
                                 outline_id) {
  
  n_days <- length(weather_series_cur[,1])
  
  plot_df <- data.frame(coordinates(data_dems$elevation[[elevation_grid_id]]))
  
  # elevation_df is to plot the contours.
  elevation_df <- data.frame(plot_df, z = getValues(data_dems$elevation[[elevation_grid_id]]))
  
  # Plot of daily SWE evolution.
  for (day_id in 1:(n_days + 1)) {
    cat("\rGenerating daily SWE plots...", day_id, "/", n_days+1)
    cells_cur <- (day_id-1) * run_params$grid_ncells + 1:(run_params$grid_ncells)
    max_swe <- 3500
    plot_df$swe <- clamp(mb_model_output$vec_swe_all[cells_cur], -Inf, max_swe)
    plot_df$snow <- as.integer(plot_df$swe > 0)
    plot_df$surf <- mb_model_output$vec_surftype_all[cells_cur]
    date_text <- format(c(weather_series_cur$timestamp, weather_series_cur$timestamp[length(weather_series_cur[,1])] + 1)[day_id], "%Y/%m/%d")
    ggplot(plot_df) +
      surf_base +
      geom_raster(aes(x = x, y = y, fill = swe, alpha = as.character(snow))) +
      scale_alpha_manual(values = c("0" = 0, "1" = 1)) +
      geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = 0.2) +
      geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = 0.15) +
      geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.2, stroke.color = "#FFFFFF", size = 1.6, min.size = 10) +
      annotate("label", x = Inf, y = Inf, hjust = 1.3, vjust = 1.5, label = date_text) +
      scale_fill_fermenter(name = "SWE [mm]", palette = "RdPu",
                           direction = 1, limits = c(0,max_swe),
                           breaks = c(100,200,400,800,1400,2000,2800,3600)) +
      guides(alpha = "none") +
      theme_void()
    ggsave(paste("output/surftype/", sprintf("%03d", day_id), ".png", sep=""), width = 5, height = 3)
  }
  
  cat("\n")
  
  
  # Plot of daily cumulative SMB.
  for (day_id in 1:(n_days+1)) {
    cat("\rGenerating daily SMB plots...", day_id, "/", n_days+1)
    cells_cur <- (day_id-1) * run_params$grid_ncells + 1:(run_params$grid_ncells)
    max_mb <- 3999
    plot_df$massbal <- mb_model_output$vec_massbal_cumul[cells_cur]
    date_text <- format(c(weather_series_cur$timestamp, weather_series_cur$timestamp[length(weather_series_cur[,1])] + 86400)[day_id], "%Y/%m/%d")
    ggplot(plot_df) +
      surf_base +
      geom_raster(aes(x = x, y = y, fill = massbal)) +
      geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = 0.2) +
      geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = 0.15) +
      geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.2, stroke.color = "#FFFFFF", size = 1.6, min.size = 10) +      annotate("label", x = Inf, y = Inf, hjust = 1.3, vjust = 1.5, label = date_text) +
      scale_fill_fermenter(name = "Cumulative\nSMB [mm w.e.]", palette = "RdBu",
                           direction = 1, limits = c(-max_mb,max_mb),
                           breaks = c(-3000,-1600,-800,-300,0,300,800,1600,3000)) +
      theme_void()
    ggsave(paste("output/massbal/", sprintf("%03d", day_id), ".png", sep=""), width = 5, height = 3)
  }
  
  cat("\n")
  
}
