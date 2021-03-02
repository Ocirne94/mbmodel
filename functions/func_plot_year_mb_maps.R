###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the routine which plots modeled mass balance maps.           #
################################################################################################### 


# NOTE: in ggplot2, the geom_sf() command
# which plots the glacier outline is forcing
# the glacier image proportions so that the
# glacier is not distorted.
# This means that the output images can get white
# margins (either above/below or left/right,
# depending on whether the glacier is larger in the
# X or in the Y coordinate).
# Without geom_sf(), the glacier is distorted
# until the image is filled.
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
  
  base_size <- 16 # For the plots.
  theme_map_massbal <- theme_void(base_size = base_size) +
                       theme(legend.position = "bottom",
                             legend.key.width = unit(3, "cm"),
                             legend.key.height = unit(0.25, "cm"),
                             legend.box.margin = margin(0,0,5,0),
                             legend.title = element_text(vjust = 1, face = "bold", size = 16),
                             legend.text = element_text(face = "bold", size = 12),
                             plot.margin = margin(0,0,0,0))
  
  contour_label_textsize <- 4
  contour_linesize <- 0.4
  outline_linesize <- 0.7
  
  palette_RdBu_ext <- c("#33000F", RColorBrewer::brewer.pal(11, "RdBu")[c(1:4,6,8:11)], "#011830")
  
  
  
  plot_df_base <- data.frame(coordinates(data_dems$elevation[[elevation_grid_id]]))
  elevation_df <- data.frame(plot_df_base, z = getValues(data_dems$elevation[[elevation_grid_id]]))
  
  plots <- list()
  
  #### HYDROLOGICAL YEAR ####
  mb_hydro_lab <- sprintf("%.3f",massbal_annual_values[["hydro"]] / 1000.)
  plot_df <- plot_df_base
  max_mb <- 6.0
  plot_df$massbal <- getValues(massbal_annual_maps$hydro)
  plots[[1]] <- ggplot(plot_df[data_dems$glacier_cell_ids[[elevation_grid_id]],]) +
    geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
    geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = outline_linesize) +
    coord_sf(clip = "off") +
    geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = contour_linesize) +
    geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.1, stroke.color = "#FFFFFF", size = contour_label_textsize, min.size = 10, fontface = "bold") +
    annotation_custom(grobTree(textGrob(paste0(year_cur-1, "/", year_cur),
                                        x=0.05,  y=1.15, hjust=0, gp = gpar(fontsize = 2 * base_size, fontface = "bold")))) +
    annotation_custom(grobTree(textGrob("Hydrological year: 10/01 - 09/30",
                                        x=0.05,  y=1.06, hjust=0, gp = gpar(fontsize = 1 * base_size, fontface = "bold")))) +
    annotation_custom(grobTree(textGrob(bquote(bold(b[n]*" = "*.(mb_hydro_lab)*" m w.e.")),
                                        x = 0.05, y = 1.0, hjust = 0, gp = gpar(fontsize = 1 * base_size)))) +
    labs(title    = " ", # Empty title to preserve spacing. We add the real title just above, with annotation_custom().
         subtitle = " ") +
    scale_fill_stepsn(name = "SMB [m w.e.]", colors = palette_RdBu_ext,
                      limits = max_mb*c(-1,1),
                      breaks = c(-5000,-3000,-2000,-1000,-400,0,400,1000,2000,3000,5000)/1000) +
    theme_map_massbal
  # ggsave(file.path(run_params$output_dirname, paste0(year_cur, "_hydro.pdf")), width = 3, height = 3)

  
  
  #### MEASUREMENT PERIOD - ANNUAL ####
  mb_meas_period_annual_lab <- paste(format(massbal_annual_meas_period, "%m/%d"), collapse = " - ")
  mb_meas_annual_lab <- sprintf("%.3f",massbal_annual_values[["meas_period"]] / 1000.)
  plot_df <- plot_df_base
  max_mb <- 6.0
  plot_df$massbal <- getValues(massbal_annual_maps$meas_period)
  plots[[2]] <- ggplot(plot_df[data_dems$glacier_cell_ids[[elevation_grid_id]],]) +
    geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
    geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = outline_linesize) +
    coord_sf(clip = "off") +
    geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = contour_linesize) +
    geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.1, stroke.color = "#FFFFFF", size = contour_label_textsize, min.size = 10, fontface = "bold") +
    annotation_custom(grobTree(textGrob(paste0(year_cur-1, "/", year_cur),
                                        x=0.05,  y=1.15, hjust=0, gp = gpar(fontsize = 2 * base_size, fontface = "bold")))) +
    annotation_custom(grobTree(textGrob(paste0("Measurement period (annual): ", mb_meas_period_annual_lab),
                                        x=0.05,  y=1.06, hjust=0, gp = gpar(fontsize = 1 * base_size, fontface = "bold")))) +
    annotation_custom(grobTree(textGrob(bquote(bold(b[n]*" = "*.(mb_meas_annual_lab)*" m w.e.")),
                                        x = 0.05, y = 1.0, hjust = 0, gp = gpar(fontsize = 1 * base_size)))) +
    labs(title    = " ", # Empty title to preserve spacing. We add the real title just above, with annotation_custom().
         subtitle = " ") +
    scale_fill_stepsn(name = "SMB [m w.e.]", colors = palette_RdBu_ext,
                      limits = max_mb*c(-1,1),
                      breaks = c(-5000,-3000,-2000,-1000,-400,0,400,1000,2000,3000,5000)/1000) +
    theme_map_massbal
  # ggsave(file.path(run_params$output_dirname, paste0(year_cur, "_meas_annual.pdf")), width = 3, height = 3)

  
  
  #### MEASUREMENT PERIOD CORRECTED - ANNUAL ####
  mb_meas_corr_annual_lab <- sprintf("%.3f",massbal_annual_values[["meas_period_corr"]] / 1000.)
  plot_df <- plot_df_base
  max_mb <- 6.0
  plot_df$massbal <- getValues(massbal_annual_maps$meas_period_corr)
  plots[[3]] <- ggplot(plot_df[data_dems$glacier_cell_ids[[elevation_grid_id]],]) +
    geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
    geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = outline_linesize) +
    coord_sf(clip = "off") +
    geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = contour_linesize) +
    geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.1, stroke.color = "#FFFFFF", size = contour_label_textsize, min.size = 10, fontface = "bold") +
    annotation_custom(grobTree(textGrob(paste0(year_cur-1, "/", year_cur),
                                        x=0.05,  y=1.15, hjust=0, gp = gpar(fontsize = 2 * base_size, fontface = "bold")))) +
    annotation_custom(grobTree(textGrob(paste0("Measurement period (annual, corrected): ", mb_meas_period_annual_lab),
                                        x=0.05,  y=1.06, hjust=0, gp = gpar(fontsize = 1 * base_size, fontface = "bold")))) +
    annotation_custom(grobTree(textGrob(bquote(bold(b[n]*" = "*.(mb_meas_corr_annual_lab)*" m w.e.")),
                                        x = 0.05, y = 1.0, hjust = 0, gp = gpar(fontsize = 1 * base_size)))) +
    labs(title    = " ", # Empty title to preserve spacing. We add the real title just above, with annotation_custom().
         subtitle = " ") +
    scale_fill_stepsn(name = "SMB [m w.e.]", colors = palette_RdBu_ext,
                      limits = max_mb*c(-1,1),
                      breaks = c(-5000,-3000,-2000,-1000,-400,0,400,1000,2000,3000,5000)/1000) +
    theme_map_massbal
  # ggsave(file.path(run_params$output_dirname, paste0(year_cur, "_meas_annual_corr.pdf")), width = 3, height = 3)

  
  
  #### USER-DEFINED FIXED PERIOD - ANNUAL ####
  mb_fixed_period_annual_lab <- paste(run_params$massbal_fixed_annual_start, run_params$massbal_fixed_annual_end, sep = " - ")
  mb_fixed_annual_lab <- sprintf("%.3f",massbal_annual_values[["fixed"]] / 1000.)
  plot_df <- plot_df_base
  max_mb <- 6.0
  plot_df$massbal <- getValues(massbal_annual_maps$fixed)
  plots[[4]] <- ggplot(plot_df[data_dems$glacier_cell_ids[[elevation_grid_id]],]) +
    geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
    geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = outline_linesize) +
    coord_sf(clip = "off") +
    geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = contour_linesize) +
    geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.1, stroke.color = "#FFFFFF", size = contour_label_textsize, min.size = 10, fontface = "bold") +
    annotation_custom(grobTree(textGrob(paste0(year_cur-1, "/", year_cur),
                                        x=0.05,  y=1.15, hjust=0, gp = gpar(fontsize = 2 * base_size, fontface = "bold")))) +
    annotation_custom(grobTree(textGrob(paste0("Fixed period (annual): ", mb_fixed_period_annual_lab),
                                        x=0.05,  y=1.06, hjust=0, gp = gpar(fontsize = 1 * base_size, fontface = "bold")))) +
    annotation_custom(grobTree(textGrob(bquote(bold(b[n]*" = "*.(mb_fixed_annual_lab)*" m w.e.")),
                                        x = 0.05, y = 1.0, hjust = 0, gp = gpar(fontsize = 1 * base_size)))) +
    labs(title    = " ", # Empty title to preserve spacing. We add the real title just above, with annotation_custom().
         subtitle = " ") +
    scale_fill_stepsn(name = "SMB [m w.e.]", colors = palette_RdBu_ext,
                      limits = max_mb*c(-1,1),
                      breaks = c(-5000,-3000,-2000,-1000,-400,0,400,1000,2000,3000,5000)/1000) +
    theme_map_massbal
  # ggsave(file.path(run_params$output_dirname, paste0(year_cur, "_fixed_annual.pdf")), width = 3, height = 3)

  
  
  #### USER-DEFINED FIXED PERIOD - WINTER ####
  mb_fixed_period_winter_lab <- paste(run_params$massbal_fixed_winter_start, run_params$massbal_fixed_winter_end, sep = " - ")
  mb_fixed_winter_lab <- sprintf("%.3f",massbal_winter_values[["fixed"]] / 1000.)
  plot_df <- plot_df_base
  max_mb <- 6.0
  plot_df$massbal <- getValues(massbal_winter_maps$fixed)
  plots[[5]] <- ggplot(plot_df[data_dems$glacier_cell_ids[[elevation_grid_id]],]) +
    geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
    geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = outline_linesize) +
    coord_sf(clip = "off") +
    geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = contour_linesize) +
    geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.1, stroke.color = "#FFFFFF", size = contour_label_textsize, min.size = 10, fontface = "bold") +
    annotation_custom(grobTree(textGrob(paste0(year_cur-1, "/", year_cur),
                                        x=0.05,  y=1.15, hjust=0, gp = gpar(fontsize = 2 * base_size, fontface = "bold")))) +
    annotation_custom(grobTree(textGrob(paste0("Fixed period (winter): ", mb_fixed_period_winter_lab),
                                        x=0.05,  y=1.06, hjust=0, gp = gpar(fontsize = 1 * base_size, fontface = "bold")))) +
    annotation_custom(grobTree(textGrob(bquote(bold(b[w]*" = "*.(mb_fixed_winter_lab)*" m w.e.")),
                                        x = 0.05, y = 1.0, hjust = 0, gp = gpar(fontsize = 1 * base_size)))) +
    labs(title    = " ", # Empty title to preserve spacing. We add the real title just above, with annotation_custom().
         subtitle = " ") +
    scale_fill_stepsn(name = "SMB [m w.e.]", colors = palette_RdBu_ext,
                      limits = max_mb*c(-1,1),
                      breaks = c(-5000,-3000,-2000,-1000,-400,0,400,1000,2000,3000,5000)/1000) +
    theme_map_massbal
  # ggsave(file.path(run_params$output_dirname, paste0(year_cur, "_fixed_winter.pdf")), width = 3, height = 3)
  
  
  if (process_winter) {
    #### MEASUREMENT PERIOD - WINTER ####
    mb_meas_period_winter_lab <- paste(format(massbal_winter_meas_period, "%m/%d"), collapse = " - ")
    mb_meas_winter_lab <- sprintf("%.3f",massbal_winter_values[["meas_period"]] / 1000.)
    plot_df <- plot_df_base
    max_mb <- 6.0
    plot_df$massbal <- getValues(massbal_winter_maps$meas_period)
    plots[[6]] <- ggplot(plot_df[data_dems$glacier_cell_ids[[elevation_grid_id]],]) +
      geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
      geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = outline_linesize) +
      coord_sf(clip = "off") +
      geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = contour_linesize) +
      geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.1, stroke.color = "#FFFFFF", size = contour_label_textsize, min.size = 10, fontface = "bold") +
      annotation_custom(grobTree(textGrob(paste0(year_cur-1, "/", year_cur),
                                          x=0.05,  y=1.15, hjust=0, gp = gpar(fontsize = 2 * base_size, fontface = "bold")))) +
      annotation_custom(grobTree(textGrob(paste0("Measurement period (winter): ", mb_meas_period_winter_lab),
                                          x=0.05,  y=1.06, hjust=0, gp = gpar(fontsize = 1 * base_size, fontface = "bold")))) +
      annotation_custom(grobTree(textGrob(bquote(bold(b[w]*" = "*.(mb_meas_winter_lab)*" m w.e.")),
                                          x = 0.05, y = 1.0, hjust = 0, gp = gpar(fontsize = 1 * base_size)))) +
      labs(title    = " ", # Empty title to preserve spacing. We add the real title just above, with annotation_custom().
           subtitle = " ") +
      scale_fill_stepsn(name = "SMB [m w.e.]", colors = palette_RdBu_ext,
                        limits = max_mb*c(-1,1),
                        breaks = c(-5000,-3000,-2000,-1000,-400,0,400,1000,2000,3000,5000)/1000) +
      theme_map_massbal
    # ggsave(file.path(run_params$output_dirname, paste0(year_cur, "_meas_winter.pdf")), width = 3, height = 3)
  }
  
  return(plots)
  
}
