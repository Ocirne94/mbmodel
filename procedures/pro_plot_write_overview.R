###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file contains the code to plot and write the model output overview.        #
###################################################################################################

# NOTE: this code is source()'d as part of main.R.
# We put code here just to make it more organized.

cat("** Drawing overview plots... **\n")

df_overview$mb_cumul <- cumsum(df_overview$mb_annual_meas_corr)
func_plot_overview(df_overview)
df_overview_out <- data.frame(year = df_overview$year,
                              apply(df_overview[,2:7], 2, sprintf, fmt="%.3f"),
                              df_overview[,8],
                              sprintf("%.1f", df_overview[,9]),
                              apply(df_overview[,10:13], 2, sprintf, fmt="%.3f"),
                              df_overview[,14],
                              sprintf("%.1f", df_overview[,15]))
names(df_overview_out) <- names(df_overview)
write.csv(df_overview_out,
          file.path(run_params$output_dirname, "overview.csv"),
          quote = FALSE,
          row.names = FALSE)

# Save to a separate file the annual maps of final annual mass balance
# (same as already saved within each annual PDF).
overview_areaplot <- ggarrange(plotlist = overview_areaplots, ncol = 1, nrow = 1, align = "hv")
suppressMessages(ggexport(overview_areaplot,
                          filename = file.path(run_params$output_dirname, "overview_areaplot.pdf"),
                          width = 21 * run_params$size_mult,
                          height = 29.7 * run_params$size_mult)) 
