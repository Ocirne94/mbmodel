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
                                   dem_grid_id,
                                   outline_id,
                                   massbal_annual_maps,
                                   massbal_winter_maps,
                                   massbal_annual_values,
                                   massbal_winter_values,
                                   massbal_annual_meas_period,
                                   massbal_winter_meas_period,
                                   massbal_annual_meas_cur,
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
  max_mb <- abs(2*run_params$mb_colorscale_breaks[1] - run_params$mb_colorscale_breaks[2]) # Values exceeding +/- max_mb will be clamped.
  
  plot_df_base <- data.frame(coordinates(data_dems$elevation[[dem_grid_id]]))
  elevation_df <- data.frame(plot_df_base, z = getValues(data_dems$elevation[[dem_grid_id]]))
  
  plots <- list()
  
  #### HYDROLOGICAL YEAR ####
  mb_hydro_lab <- sprintf("%.3f",massbal_annual_values[["hydro"]] / 1000.)
  plot_df <- plot_df_base
  plot_df$massbal <- getValues(massbal_annual_maps$hydro)
  plots[[1]] <- ggplot(plot_df[data_dems$glacier_cell_ids[[dem_grid_id]],]) +
    geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
    geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = outline_linesize) +
    coord_sf(clip = "off") +
    geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = contour_linesize) +
    geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.1, stroke.color = "#FFFFFF", size = contour_label_textsize, min.size = 15, fontface = "bold") +
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
                      breaks = run_params$mb_colorscale_breaks) +
    theme_map_massbal

  
  
  #### MEASUREMENT PERIOD - ANNUAL ####
  mb_meas_period_annual_lab <- paste(format(massbal_annual_meas_period, "%m/%d"), collapse = " - ")
  mb_meas_annual_lab <- sprintf("%.3f",massbal_annual_values[["meas_period"]] / 1000.)
  plot_df <- plot_df_base
  plot_df$massbal <- getValues(massbal_annual_maps$meas_period)
  plots[[2]] <- ggplot(plot_df[data_dems$glacier_cell_ids[[dem_grid_id]],]) +
    geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
    geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = outline_linesize) +
    coord_sf(clip = "off") +
    geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = contour_linesize) +
    geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.1, stroke.color = "#FFFFFF", size = contour_label_textsize, min.size = 15, fontface = "bold") +
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
                      breaks = run_params$mb_colorscale_breaks) +
    theme_map_massbal

  
  
  #### MEASUREMENT PERIOD - ANNUAL, WITH STAKES ####
  # Also RMS.
  plots[[3]] <- ggplot(plot_df[data_dems$glacier_cell_ids[[dem_grid_id]],]) +
    geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
    geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = outline_linesize) +
    coord_sf(clip = "off") +
    geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = contour_linesize) +
    geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.1, stroke.color = "#FFFFFF", size = contour_label_textsize, min.size = 15, fontface = "bold") +
    geom_point(data = massbal_annual_meas_cur, aes(x = x, y = y), shape = 3, stroke = 1.5, size = 0) +
    geom_shadowtext(data = massbal_annual_meas_cur, aes(x = x, y = y, label = sprintf("%.2f", massbal_standardized/1e3)), size = 3, fontface = "bold", color = "#000000", hjust = -0.12, vjust = -0.12, bg.color = "#FFFFFF") +
    annotation_custom(grobTree(textGrob(paste0(year_cur-1, "/", year_cur),
                                        x=0.05,  y=1.15, hjust=0, gp = gpar(fontsize = 2 * base_size, fontface = "bold")))) +
    annotation_custom(grobTree(textGrob(paste0("Measurement period (annual): ", mb_meas_period_annual_lab),
                                        x=0.05,  y=1.06, hjust=0, gp = gpar(fontsize = 1 * base_size, fontface = "bold")))) +
    annotation_custom(grobTree(textGrob(bquote(bold(b[n]*" = "*.(mb_meas_annual_lab)*" m w.e.")),
                                        x = 0.05, y = 1.0, hjust = 0, gp = gpar(fontsize = 1 * base_size)))) +
    annotation_custom(grobTree(textGrob(paste0("RMS: ", sprintf("%.3f", mod_output_annual_cur$global_rms/1e3), " m w.e."),
                                        x = 0.05, y = 0.94, hjust = 0, gp = gpar(fontsize = 1 * base_size, fontface = "bold")))) +
    labs(title    = " ", # Empty title to preserve spacing. We add the real title just above, with annotation_custom().
         subtitle = " ") +
    scale_fill_stepsn(name = "SMB [m w.e.]", colors = palette_RdBu_ext,
                      limits = max_mb*c(-1,1),
                      breaks = run_params$mb_colorscale_breaks) +
    theme_map_massbal
    
  
  
  
  #### MEASUREMENT PERIOD - ANNUAL CORRECTED ####
  mb_meas_corr_annual_lab <- sprintf("%.3f",massbal_annual_values[["meas_period_corr"]] / 1000.)
  plot_df <- plot_df_base
  plot_df$massbal <- getValues(massbal_annual_maps$meas_period_corr)
  plots[[4]] <- ggplot(plot_df[data_dems$glacier_cell_ids[[dem_grid_id]],]) +
    geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
    geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = outline_linesize) +
    coord_sf(clip = "off") +
    geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = contour_linesize) +
    geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.1, stroke.color = "#FFFFFF", size = contour_label_textsize, min.size = 15, fontface = "bold") +
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
                      breaks = run_params$mb_colorscale_breaks) +
    theme_map_massbal

  
  
  #### MEASUREMENT PERIOD - ANNUAL CORRECTED, WITH STAKES ####
  rmse_bandcorr <- sqrt(mean((massbal_annual_meas_cur$massbal_standardized - extract(massbal_annual_maps$meas_period_corr, cbind(massbal_annual_meas_cur$x, massbal_annual_meas_cur$y), method = "bilinear"))^2))
  plots[[5]] <- ggplot(plot_df[data_dems$glacier_cell_ids[[dem_grid_id]],]) +
    geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
    geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = outline_linesize) +
    coord_sf(clip = "off") +
    geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = contour_linesize) +
    geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.1, stroke.color = "#FFFFFF", size = contour_label_textsize, min.size = 15, fontface = "bold") +
    geom_point(data = massbal_annual_meas_cur, aes(x = x, y = y), shape = 3, stroke = 1.5, size = 0) +
    geom_shadowtext(data = massbal_annual_meas_cur, aes(x = x, y = y, label = sprintf("%.2f", massbal_standardized/1e3)), size = 3, fontface = "bold", color = "#000000", hjust = -0.12, vjust = -0.12, bg.color = "#FFFFFF") +
    annotation_custom(grobTree(textGrob(paste0(year_cur-1, "/", year_cur),
                                        x=0.05,  y=1.15, hjust=0, gp = gpar(fontsize = 2 * base_size, fontface = "bold")))) +
    annotation_custom(grobTree(textGrob(paste0("Measurement period (annual, corrected): ", mb_meas_period_annual_lab),
                                        x=0.05,  y=1.06, hjust=0, gp = gpar(fontsize = 1 * base_size, fontface = "bold")))) +
    annotation_custom(grobTree(textGrob(bquote(bold(b[n]*" = "*.(mb_meas_corr_annual_lab)*" m w.e.")),
                                        x = 0.05, y = 1.0, hjust = 0, gp = gpar(fontsize = 1 * base_size)))) +
    annotation_custom(grobTree(textGrob(paste0("RMS: ", sprintf("%.3f", rmse_bandcorr/1e3), " m w.e."),
                                        x = 0.05, y = 0.94, hjust = 0, gp = gpar(fontsize = 1 * base_size, fontface = "bold")))) +
    labs(title    = " ", # Empty title to preserve spacing. We add the real title just above, with annotation_custom().
         subtitle = " ") +
    scale_fill_stepsn(name = "SMB [m w.e.]", colors = palette_RdBu_ext,
                      limits = max_mb*c(-1,1),
                      breaks = run_params$mb_colorscale_breaks) +
    theme_map_massbal
  
  
  
  #### USER-DEFINED FIXED PERIOD - ANNUAL ####
  mb_fixed_period_annual_lab <- paste(run_params$massbal_fixed_annual_start, run_params$massbal_fixed_annual_end, sep = " - ")
  mb_fixed_annual_lab <- sprintf("%.3f",massbal_annual_values[["fixed"]] / 1000.)
  plot_df <- plot_df_base
  plot_df$massbal <- getValues(massbal_annual_maps$fixed)
  plots[[6]] <- ggplot(plot_df[data_dems$glacier_cell_ids[[dem_grid_id]],]) +
    geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
    geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = outline_linesize) +
    coord_sf(clip = "off") +
    geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = contour_linesize) +
    geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.1, stroke.color = "#FFFFFF", size = contour_label_textsize, min.size = 15, fontface = "bold") +
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
                      breaks = run_params$mb_colorscale_breaks) +
    theme_map_massbal

  
  
  #### USER-DEFINED FIXED PERIOD - WINTER ####
  mb_fixed_period_winter_lab <- paste(run_params$massbal_fixed_winter_start, run_params$massbal_fixed_winter_end, sep = " - ")
  mb_fixed_winter_lab <- sprintf("%.3f",massbal_winter_values[["fixed"]] / 1000.)
  plot_df <- plot_df_base
  plot_df$massbal <- getValues(massbal_winter_maps$fixed)
  plots[[7]] <- ggplot(plot_df[data_dems$glacier_cell_ids[[dem_grid_id]],]) +
    geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
    geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = outline_linesize) +
    coord_sf(clip = "off") +
    geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = contour_linesize) +
    geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.1, stroke.color = "#FFFFFF", size = contour_label_textsize, min.size = 15, fontface = "bold") +
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
                      breaks = run_params$mb_colorscale_breaks) +
    theme_map_massbal

  
  if (process_winter) {
    #### MEASUREMENT PERIOD - WINTER ####
    mb_meas_period_winter_lab <- paste(format(massbal_winter_meas_period, "%m/%d"), collapse = " - ")
    mb_meas_winter_lab <- sprintf("%.3f",massbal_winter_values[["meas_period"]] / 1000.)
    plot_df <- plot_df_base
    plot_df$massbal <- getValues(massbal_winter_maps$meas_period)
    plots[[8]] <- ggplot(plot_df[data_dems$glacier_cell_ids[[dem_grid_id]],]) +
      geom_raster(aes(x = x, y = y, fill = massbal/1000)) +
      geom_sf(data = as(data_outlines$outlines[[outline_id]], "sf"), fill = NA, color = "#202020", size = outline_linesize) +
      coord_sf(clip = "off") +
      geom_contour(data = elevation_df, aes(x = x, y = y, z = z), color = "#202020", size = contour_linesize) +
      geom_text_contour(data = elevation_df, aes(x = x, y = y, z = z), check_overlap = TRUE, stroke = 0.1, stroke.color = "#FFFFFF", size = contour_label_textsize, min.size = 15, fontface = "bold") +
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
                        breaks = run_params$mb_colorscale_breaks) +
      theme_map_massbal
  }
  
  return(plots)
  
}
