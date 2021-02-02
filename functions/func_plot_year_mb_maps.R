###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the routine which plots modeled mass balance maps.           #
################################################################################################### 

func_plot_year_mb_maps <- function(run_params,
                                   year_cur,
                                   data_dems,
                                   data_outlines,
                                   elevation_grid_id,
                                   outline_id,
                                   massbal_annual_maps,
                                   massbal_winter_maps,
                                   massbal_annual_values,
                                   massbal_winter_values,
                                   massbal_annual_meas_period,
                                   massbal_winter_meas_period,
                                   process_winter) {
  
  dir.create(file.path("output", run_params$name_glacier, "massbal"), recursive = TRUE)
  
  base_size <- 6
  
  plot_df_base <- data.frame(coordinates(data_dems$elevation[[elevation_grid_id]]))
  elevation_df <- data.frame(plot_df_base, z = getValues(data_dems$elevation[[elevation_grid_id]]))
  
  palette_RdBu_ext <- c("#33000F", RColorBrewer::brewer.pal(11, "RdBu")[c(1:4,6,8:11)], "#011830")
  
  
  #### HYDROLOGICAL YEAR ####
  mb_hydro_lab <- sprintf("%.3f",massbal_annual_values[["hydro"]] / 1000.)
  plot_df <- plot_df_base
  max_mb <- 6.0
  plot_df$massbal <- getValues(massbal_annual_maps$hydro)
  ggplot(plot_df[data_dems$glacier_cell_ids[[elevation_grid_id]],]) +
    geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
    geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = 0.2) +
    geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = 0.15) +
    geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.2, stroke.color = "#FFFFFF", size = 1.6, min.size = 10) +
    labs(title    = paste0("Hydrological year ", year_cur, ": 10/01 - 09/30"),
         subtitle = bquote(b[n]*" = "*.(mb_hydro_lab)*" m w.e.")) +
    scale_fill_stepsn(name = "SMB [m w.e.]", colors = palette_RdBu_ext,
                      limits = max_mb*c(-1,1),
                      breaks = c(-5000,-3000,-2000,-1000,-400,0,400,1000,2000,3000,5000)/1000) +
    theme_void(base_size = base_size) +
    theme(legend.position = "bottom",
          legend.key.width = unit(30, "pt"),
          legend.key.height = unit(2, "pt"),
          legend.box.margin = margin(0,0,5,0),
          legend.title = element_text(vjust = 1),
          plot.margin = margin(0,0,0,0))
  ggsave(paste0("output/", run_params$name_glacier, "/massbal/", year_cur, "_hydro.png"), width = 3, height = 3)
  
  
  
  #### MEASUREMENT PERIOD - ANNUAL ####
  mb_meas_period_annual_lab <- paste(format(massbal_annual_meas_period, "%m/%d"), collapse = " - ")
  mb_meas_annual_lab <- sprintf("%.3f",massbal_annual_values[["meas_period"]] / 1000.)
  plot_df <- plot_df_base
  max_mb <- 6.0
  plot_df$massbal <- getValues(massbal_annual_maps$meas_period)
  ggplot(plot_df[data_dems$glacier_cell_ids[[elevation_grid_id]],]) +
    geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
    geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = 0.2) +
    geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = 0.15) +
    geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.2, stroke.color = "#FFFFFF", size = 1.6, min.size = 10) +
    labs(title    = paste0("Measurement period (annual) ", year_cur, ": ", mb_meas_period_annual_lab),
         subtitle = bquote(b[n]*" = "*.(mb_meas_annual_lab)*" m w.e.")) +
    scale_fill_stepsn(name = "SMB [m w.e.]", colors = palette_RdBu_ext,
                      limits = max_mb*c(-1,1),
                      breaks = c(-5000,-3000,-2000,-1000,-400,0,400,1000,2000,3000,5000)/1000) +
    theme_void(base_size = base_size) +
    theme(legend.position = "bottom",
          legend.key.width = unit(30, "pt"),
          legend.key.height = unit(2, "pt"),
          legend.box.margin = margin(0,0,5,0),
          legend.title = element_text(vjust = 1),
          plot.margin = margin(0,0,0,0))
  ggsave(paste0("output/", run_params$name_glacier, "/massbal/", year_cur, "_meas_annual.png"), width = 3, height = 3)
  
  
  #### MEASUREMENT PERIOD CORRECTED - ANNUAL ####
  mb_meas_corr_annual_lab <- sprintf("%.3f",massbal_annual_values[["meas_period_corr"]] / 1000.)
  plot_df <- plot_df_base
  max_mb <- 6.0
  plot_df$massbal <- getValues(massbal_annual_maps$meas_period_corr)
  ggplot(plot_df[data_dems$glacier_cell_ids[[elevation_grid_id]],]) +
    geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
    geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = 0.2) +
    geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = 0.15) +
    geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.2, stroke.color = "#FFFFFF", size = 1.6, min.size = 10) +
    labs(title    = paste0("Measurement period (annual, corrected) ", year_cur, ": ", mb_meas_period_annual_lab),
         subtitle = bquote(b[n]*" = "*.(mb_meas_corr_annual_lab)*" m w.e.")) +
    scale_fill_stepsn(name = "SMB [m w.e.]", colors = palette_RdBu_ext,
                      limits = max_mb*c(-1,1),
                      breaks = c(-5000,-3000,-2000,-1000,-400,0,400,1000,2000,3000,5000)/1000) +
    theme_void(base_size = base_size) +
    theme(legend.position = "bottom",
          legend.key.width = unit(30, "pt"),
          legend.key.height = unit(2, "pt"),
          legend.box.margin = margin(0,0,5,0),
          legend.title = element_text(vjust = 1),
          plot.margin = margin(0,0,0,0))
  ggsave(paste0("output/", run_params$name_glacier, "/massbal/", year_cur, "_meas_annual_corr.png"), width = 3, height = 3)
  
  
  
  #### USER-DEFINED FIXED PERIOD - ANNUAL ####
  mb_fixed_period_annual_lab <- paste(run_params$massbal_fixed_annual_start, run_params$massbal_fixed_annual_end, sep = " - ")
  mb_fixed_annual_lab <- sprintf("%.3f",massbal_annual_values[["fixed"]] / 1000.)
  plot_df <- plot_df_base
  max_mb <- 6.0
  plot_df$massbal <- getValues(massbal_annual_maps$fixed)
  ggplot(plot_df[data_dems$glacier_cell_ids[[elevation_grid_id]],]) +
    geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
    geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = 0.2) +
    geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = 0.15) +
    geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.2, stroke.color = "#FFFFFF", size = 1.6, min.size = 10) +
    labs(title    = paste0("Fixed annual period ", year_cur, ": ", mb_fixed_period_annual_lab),
         subtitle = bquote(b[n]*" = "*.(mb_fixed_annual_lab)*" m w.e.")) +
    scale_fill_stepsn(name = "SMB [m w.e.]", colors = palette_RdBu_ext,
                      limits = max_mb*c(-1,1),
                      breaks = c(-5000,-3000,-2000,-1000,-400,0,400,1000,2000,3000,5000)/1000) +
    theme_void(base_size = base_size) +
    theme(legend.position = "bottom",
          legend.key.width = unit(30, "pt"),
          legend.key.height = unit(2, "pt"),
          legend.box.margin = margin(0,0,5,0),
          legend.title = element_text(vjust = 1),
          plot.margin = margin(0,0,0,0))
  ggsave(paste0("output/", run_params$name_glacier, "/massbal/", year_cur, "_fixed_annual.png"), width = 3, height = 3)
  
  
  #### USER-DEFINED FIXED PERIOD - WINTER ####
  mb_fixed_period_winter_lab <- paste(run_params$massbal_fixed_winter_start, run_params$massbal_fixed_winter_end, sep = " - ")
  mb_fixed_winter_lab <- sprintf("%.3f",massbal_winter_values[["fixed"]] / 1000.)
  plot_df <- plot_df_base
  max_mb <- 6.0
  plot_df$massbal <- getValues(massbal_winter_maps$fixed)
  ggplot(plot_df[data_dems$glacier_cell_ids[[elevation_grid_id]],]) +
    geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
    geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = 0.2) +
    geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = 0.15) +
    geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.2, stroke.color = "#FFFFFF", size = 1.6, min.size = 10) +
    labs(title    = paste0("Fixed winter period ", year_cur, ": ", mb_fixed_period_winter_lab),
         subtitle = bquote(b[n]*" = "*.(mb_fixed_winter_lab)*" m w.e.")) +
    scale_fill_stepsn(name = "SMB [m w.e.]", colors = palette_RdBu_ext,
                      limits = max_mb*c(-1,1),
                      breaks = c(-5000,-3000,-2000,-1000,-400,0,400,1000,2000,3000,5000)/1000) +
    theme_void(base_size = base_size) +
    theme(legend.position = "bottom",
          legend.key.width = unit(30, "pt"),
          legend.key.height = unit(2, "pt"),
          legend.box.margin = margin(0,0,5,0),
          legend.title = element_text(vjust = 1),
          plot.margin = margin(0,0,0,0))
  ggsave(paste0("output/", run_params$name_glacier, "/massbal/", year_cur, "_fixed_winter.png"), width = 3, height = 3)
  
  
  if (process_winter) {
    #### MEASUREMENT PERIOD - WINTER ####
    mb_meas_period_winter_lab <- paste(format(massbal_winter_meas_period, "%m/%d"), collapse = " - ")
    mb_meas_winter_lab <- sprintf("%.3f",massbal_winter_values[["meas_period"]] / 1000.)
    plot_df <- plot_df_base
    max_mb <- 6.0
    plot_df$massbal <- getValues(massbal_winter_maps$meas_period)
    ggplot(plot_df[data_dems$glacier_cell_ids[[elevation_grid_id]],]) +
      geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
      geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = 0.2) +
      geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = 0.15) +
      geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.2, stroke.color = "#FFFFFF", size = 1.6, min.size = 10) +
      labs(title    = paste0("Measurement period (winter) ", year_cur, ": ", mb_meas_period_winter_lab),
           subtitle = bquote(b[n]*" = "*.(mb_meas_winter_lab)*" m w.e.")) +
      scale_fill_stepsn(name = "SMB [m w.e.]", colors = palette_RdBu_ext,
                        limits = max_mb*c(-1,1),
                        breaks = c(-5000,-3000,-2000,-1000,-400,0,400,1000,2000,3000,5000)/1000) +
      theme_void(base_size = base_size) +
      theme(legend.position = "bottom",
            legend.key.width = unit(30, "pt"),
            legend.key.height = unit(2, "pt"),
            legend.box.margin = margin(0,0,5,0),
            legend.title = element_text(vjust = 1),
            plot.margin = margin(0,0,0,0))
    ggsave(paste0("output/", run_params$name_glacier, "/massbal/", year_cur, "_meas_winter.png"), width = 3, height = 3)
  }
  
}
