###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the generation of the final overview plots.                  #
###################################################################################################  

func_plot_overview <- function(df_overview) {
  
  base_size <- 16 # For the plots
  
  theme_overview_plots <- theme_bw(base_size = base_size) +
                          theme(axis.title.x = element_blank(),
                                plot.title = element_text(hjust = 0.5),
                                text = element_text(face = "bold"))
  
  plots <- list()
  
  # Show 4-5 years on the x axis, using only integer values.
  x_breaks <- seq(df_overview$year[1], df_overview$year[length(df_overview$year)], by = max(1, floor(length(df_overview$year) / 4)))
  
  # Do we have a single year to show?
  single_year <- nrow(df_overview) == 1
  
  # Time series of annual mass balance over the measurement
  # period, corrected within elevation bands.
  # Also horizontal line with mean over the period.
  # If we model just one year, add point plot so that something is visible.
  single_year_point <- NULL
  if (single_year) {single_year_point <- geom_point(aes(x = year, y = mb_annual_meas_corr))}
  plots[[1]] <- ggplot(df_overview) +
    geom_line(aes(x = year, y = mb_annual_meas_corr), size = 1) +
    geom_segment(x = df_overview$year[1], xend = df_overview$year[length(df_overview$year)],
                 y = mean(df_overview$mb_annual_meas_corr), yend = mean(df_overview$mb_annual_meas_corr),
                 linetype = "dashed", size = 1) +
    single_year_point +
    scale_x_continuous(breaks = x_breaks) +
    ylab("Mass balance [m w.e.]") +
    ggtitle("Final mass balance") +
    theme_overview_plots

  
  # Time series of other annual mass balances:
  # over the measurement period with no correction,
  # over the hydrological year,
  # over a fixed (user-defined) period.
  # If we model just one year, add point plot so that something is visible.
  single_year_point1 <- NULL
  single_year_point2 <- NULL
  single_year_point3 <- NULL
  if (single_year) {
    single_year_point1 <- geom_point(aes(x = year, y = mb_annual_meas), color = "#FF00FF")
    single_year_point2 <- geom_point(aes(x = year, y = mb_annual_hydro), color = "#0000FF")
    single_year_point3 <- geom_point(aes(x = year, y = mb_annual_fixed), color = "#00FFFF")
  }
  plots[[2]] <- ggplot(df_overview) +
    geom_line(aes(x = year, y = mb_annual_meas), color = "#FF00FF", size = 1) +
    geom_line(aes(x = year, y = mb_annual_hydro), color = "#0000FF", size = 1) +
    geom_line(aes(x = year, y = mb_annual_fixed), color = "#00FFFF", size = 1) +
    single_year_point1 +
    single_year_point2 +
    single_year_point3 +
    ylab("Mass balance [m w.e.]") +
    scale_y_continuous(expand = expansion(0.3, 0)) +
    scale_x_continuous(breaks = x_breaks) +
    ggtitle("Annual mass balance") +
    annotation_custom(grobTree(textGrob("Measurement period", x=0.05, y = 0.19, hjust = 0,
                                        gp=gpar(col="#FF00FF", fontsize = base_size * 1., fontface="bold")))) +
    annotation_custom(grobTree(textGrob("Hydrological year", x=0.05, y = 0.12, hjust = 0,
                                        gp=gpar(col="#0000FF", fontsize = base_size * 1., fontface="bold")))) +
    annotation_custom(grobTree(textGrob("Fixed period", x=0.05, y = 0.05, hjust = 0,
                                        gp=gpar(col="#00FFFF", fontsize = base_size * 1., fontface="bold")))) +
    theme_overview_plots
  
  
  # Time series of winter mass balances:
  # over a fixed (user-defined) winter period,
  # over the measurement period (only if winter measurements are available).
  # If we model just one year, add point plot so that something is visible.
  single_year_point1 <- NULL
  single_year_point2 <- NULL
  if (single_year) {
    single_year_point1 <- geom_point(aes(x = year, y = mb_winter_fixed), color = "#00FFFF")
    if(any(!is.na(df_overview$mb_winter_meas))) {
      single_year_point2 <- geom_point(aes(x = year, y = mb_winter_meas), color = "#FF00FF")
    }
  }
  plots[[3]] <- ggplot(df_overview) +
    {if(any(!is.na(df_overview$mb_winter_meas))) geom_line(aes(x = year, y = mb_winter_meas), color = "#FF00FF", size = 1)} +
    geom_line(aes(x = year, y = mb_winter_fixed), color = "#00FFFF", size = 1) +
    single_year_point1 +
    single_year_point2 +
    ylab("Mass balance [m w.e.]") +
    scale_y_continuous(expand = expansion(0.3, 0)) +
    scale_x_continuous(breaks = x_breaks) +
    ggtitle("Winter mass balance") +
    {if(any(!is.na(df_overview$mb_winter_meas))) annotation_custom(grobTree(textGrob("Measurement period", x=0.05, y = 0.12, hjust = 0,
                                        gp=gpar(col="#FF00FF", fontsize = base_size * 1., fontface="bold"))))} +
    annotation_custom(grobTree(textGrob("Fixed period", x=0.05, y = 0.05, hjust = 0,
                                        gp=gpar(col="#00FFFF", fontsize = base_size * 1., fontface="bold")))) +
    theme_overview_plots
  
  
  # Time series of ELA.
  # If we model just one year, add point plot so that something is visible.
  single_year_point <- NULL
  if (single_year) {single_year_point <- geom_point(aes(x = year, y = ela))}
  plots[[4]] <- ggplot(df_overview) +
    geom_line(aes(x = year, y = ela), size = 1) +
    single_year_point +
    ylab("Equilibrium Line Altitude [m a.s.l.]") +
    scale_y_continuous(expand = expansion(0.5, 0)) +
    scale_x_continuous(breaks = x_breaks) +
    ggtitle("Equilibrium Line Altitude") +
    theme_overview_plots
  
  
  # Time series of AAR.
  # If we model just one year, add point plot so that something is visible.
  single_year_point <- NULL
  if (single_year) {single_year_point <- geom_point(aes(x = year, y = aar))}
  plots[[5]] <- ggplot(df_overview) +
    geom_line(aes(x = year, y = aar), size = 1) +
    single_year_point +
    ylab("Accumulation-Area Ratio [%]") +
    scale_y_continuous(expand = expansion(0.5, 0)) +
    scale_x_continuous(breaks = x_breaks) +
    ggtitle("Accumulation-Area Ratio") +
    theme_overview_plots
  
  
  # Time series of RMSE.
  # If we model just one year, add point plot so that something is visible.
  single_year_point <- NULL
  if (single_year) {single_year_point <- geom_point(aes(x = year, y = rmse))}
  plots[[6]] <- ggplot(df_overview) +
    geom_line(aes(x = year, y = rmse), size = 1) +
    single_year_point +
    ylab("RMSE [m w.e.]") +
    scale_y_continuous(limits = c(0,NA), expand = expansion(mult = c(0, 0.2))) +
    scale_x_continuous(breaks = x_breaks) +
    ggtitle("Root-Mean-Square Error") +
    theme_overview_plots
  
  
  # Time series of the melt parameters.
  # If we model just one year, add point plot so that something is visible.
  single_year_point1 <- NULL
  single_year_point2 <- NULL
  single_year_point3 <- NULL
  if (single_year) {
    single_year_point1 <- geom_point(aes(x = year, y = rad_fact_snow), color = "#0000FF")
    single_year_point2 <- geom_point(aes(x = year, y = melt_factor), color = "#FF00FF")
    single_year_point3 <- geom_point(aes(x = year, y = rad_fact_ice), color = "#00FFFF")
  }
  plots[[7]] <- ggplot(df_overview) +
    geom_line(aes(x = year, y = rad_fact_snow), color = "#0000FF", size = 1) +
    geom_line(aes(x = year, y = melt_factor), color = "#FF00FF", size = 1.25) + # Different size since the melt factor is sometimes the same as the rad_fact_ice.
    geom_line(aes(x = year, y = rad_fact_ice), color = "#00FFFF", size = 0.5) +
    single_year_point1 +
    single_year_point2 +
    single_year_point3 +
    ylab("Parameter value [different units]") +
    scale_y_continuous(expand = expansion(0.5, 0)) +
    scale_x_continuous(breaks = x_breaks) +
    ggtitle("Melt parameters") +
    annotation_custom(grobTree(textGrob("Degree-day melt factor", x=0.05, y = 0.19, hjust = 0,
                                        gp=gpar(col="#FF00FF", fontsize = base_size * 1., fontface="bold")))) +
    annotation_custom(grobTree(textGrob("Radiation factor (snow)", x=0.05, y = 0.12, hjust = 0,
                                        gp=gpar(col="#0000FF", fontsize = base_size * 1., fontface="bold")))) +
    annotation_custom(grobTree(textGrob("Radiation factor (ice)", x=0.05, y = 0.05, hjust = 0,
                                        gp=gpar(col="#00FFFF", fontsize = base_size * 1., fontface="bold")))) +
    theme_overview_plots
  
  
  # Time series of the precipitation correction.
  # We use a slightly more complex formula for the
  # y-axis limits so that when we have a single value
  # the limits still make sense.
  # If we model just one year, add point plot so that something is visible.
  single_year_point <- NULL
  if (single_year) {single_year_point <- geom_point(aes(x = year, y = prec_corr), color = "#0000FF")}
  plots[[8]] <- ggplot(df_overview) +
    geom_line(aes(x = year, y = prec_corr), color = "#0000FF", size = 1) +
    single_year_point +
    ylab("Precipitation correction [%]") +
    scale_y_continuous(limits = c(min(min(df_overview$prec_corr), mean(df_overview$prec_corr) - 0.1 * mean(df_overview$prec_corr)), max(max(df_overview$prec_corr), mean(df_overview$prec_corr + 0.1 * df_overview$prec_corr)))) +
    scale_x_continuous(breaks = x_breaks) +
    ggtitle("Precipitation correction factor") +
    theme_overview_plots
  
  

  x_breaks_cumul <- seq(df_overview$year[1]-1, df_overview$year[length(df_overview$year)], by = max(1, floor((length(df_overview$year)+1) / 4)))
  # Time series of cumulative mass balance.
  plots[[9]] <- ggplot(data.frame(year = c(df_overview$year[1]-1, df_overview$year),
                          mb_cumul = c(0, df_overview$mb_cumul))) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
    geom_line(aes(x = year, y = mb_cumul), color = "#FF0000", size = 1) +
    geom_point(aes(x = year, y = mb_cumul), color = "#FF0000", shape = 2, size = 3, stroke = 1.2) +
    scale_y_continuous(breaks = pretty(c(0, df_overview$mb_cumul))) +
    scale_x_continuous(breaks = x_breaks_cumul) +
    ylab("Cumulative mass balance [m w.e.]") +
    ggtitle("Cumulative mass balance") +
    theme_overview_plots
  
  overview <- ggarrange(plotlist = plots, ncol = 1, nrow = 3, align = "hv")

  suppressMessages(ggexport(overview,
                            filename = file.path(run_params$output_dirname, "overview.pdf"),
                            width = 21 * run_params$size_mult,
                            height = 29.7 * run_params$size_mult))
  
}
