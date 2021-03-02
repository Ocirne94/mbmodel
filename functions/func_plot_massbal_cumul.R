###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the generation of the two plots of mass balance evolution    #
#                 over a year: overall mass balance alone, and also with accumulation and         #
#                 ablation.                                                                       #
###################################################################################################  

func_plot_massbal_cumul <- function(run_params,
                                    mod_output_annual_cur,
                                    model_annual_bounds) {
  
  # Here will go the two plots of annual mass balance
  # ((1) mass balance only, and (2) also with accumulation
  # and ablation components).
  plots_mb <- list()
  
  
  # Prepare the data for plotting.
  massbal_cumul_df <- data.frame(date = seq.Date(model_annual_bounds[1]-1, model_annual_bounds[2], by = "1 day"),
                                 mb = mod_output_annual_cur$gl_massbal_cumul)
  massbal_cumul_df$day_id <- seq_along(massbal_cumul_df$date) - (length(massbal_cumul_df$date) - as.integer(format(massbal_cumul_df$date[length(massbal_cumul_df$date)], "%j"))) + 1
  
  
  # Setup month labels.
  months_labels_all <- format(massbal_cumul_df$date, "%b")
  months_doy <- c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)
  
  months_labels_ids <- which(as.integer(format(massbal_cumul_df$date, "%j")) %in% months_doy) # Select the day at the middle of each month.
  months_labels_df <- data.frame(day_id = massbal_cumul_df$day_id[months_labels_ids],
                                 label = months_labels_all[months_labels_ids])
  # Don't add label for first month unless it is
  # represented by at least 28 days, and same for last month.
  months_cur_rle <- rle(as.integer(format(massbal_cumul_df$date, "%m")))
  if (months_cur_rle$lengths[1] < 28) { 
    months_labels_df <- months_labels_df[-1,]
  }
  if (months_cur_rle$lengths[length(months_cur_rle$lengths)] < 28) { # Same, for last month.
    months_labels_df <- months_labels_df[-length(months_labels_df[,1]),]
  }
  
  
  base_size <- 16 # For the plots
  theme_mbcumul_plots <- theme_bw(base_size = base_size) +
    theme(axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          text = element_text(face = "bold"))
  
  # Generate plot of mass balance alone.
  plots_mb[[1]] <- ggplot(massbal_cumul_df) +
    annotate("text", x = months_labels_df$day_id, y = -Inf, label = months_labels_df$label, vjust = -1, fontface = "bold") +
    geom_line(aes(x = day_id, y = mb / 1e3)) +
    scale_y_continuous(breaks = pretty(massbal_cumul_df$mb/1e3)) +
    ylab("Mass balance [m w.e.]") +
    theme_mbcumul_plots
  
  
  # Generate plot of mass balance with accumulation and ablation.
  
  
  
  return(plots_mb)
}