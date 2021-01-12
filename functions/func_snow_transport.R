###################################################################################################
# Author:         Enrico Mattea (@unifr.ch)                                                       #
# Description:    this program models the distributed mass balance of a glacier at daily          #
#                 resolution, optimizing model parameters towards the best fit with point         #
#                 mass balance measurements.                                                      #
#                 This file implements the mass transport algorithm by Gruber (2007),             #
#                 to compute avalanche redistribution of snow.                                    #
#                 Input: raster grids of elevation, slope, aspect and initial mass.               #
#                 NOTE: we compute here a (heavily parametrized/simplified) *movable* mass        #
#                 which for each cell is a fraction of the initial mass.                          #
#                 NOTE-2: the 4 directions from a cell are indexed as 1 = up, 2 = left,           #
#                 3 = right, 4 = bottom. We use lists of rasters for the processing,              #
#                 list[[1]] is the raster related to the up direction.                            #                                                                     #
# Latest update:  2021.1.11                                                                       #
################################################################################################### 



func_snow_transport <- function(elevation_raw, slope, aspect, mass_initial) {
  
  # Original IDL implementation: snow accumulation factor
  # decreased linearly from 1 to 0 between slopes of 40° and 60°.
  # Snow was uniformly redistributed over the entire grid (unphysical!).
  # Here we take a similar approach to compute the movable mass fraction:
  # it is 0 below 35° (will be customizable), and increases linearly to 1 at 60°.
  
  elevation_raw <- raster("functions/dhm_barkrak2015.grid")
  elevation <- func_preprocess_elevation(elevation_raw, run_params)
  
  grid_ncol <- ncol(elevation)
  grid_nrow <- nrow(elevation)
  
  crs(elevation) <- CRS(SRS_string="EPSG:32642")
  slope <- terrain(elevation, "slope", "degrees")
  aspect <- terrain(elevation, "aspect", "degrees")
  
  deposition_slope_lim <- 40 # [°], beta_lim in the paper: maximum angle at which snow is deposited. Will be customizable.
  deposition_mass_lim <- 4000 # [kg m-2], D_lim in the paper: maximum deposition (if the avalanche routine
    # is fast enough to be used also for each (snowy) daily timestep of the MB model, then this D_lim should
    # become a rate: [kg m-2 d-1], with a low-ish value so that avalanches during winter can be realistic).
    # ALSO, TO RUN THIS REPEATEDLY: WE SHOULD COMPUTE THE FLOW GRIDS ONLY ONCE, PUT THEM IN A STRUCTURE, AND
    # THEN PASS THEM TO A FUNCTION ACTUALLY DOING THE COMPUTATION ON THE AVAILABLE SNOW.
  movable_slope_lim_lower <- 30
  movable_slope_lim_upper <- 60
  
  movable_frac <- setValues(slope,
                            pmax(0, pmin(1, rescale(getValues(slope),
                                                    to = c(0, 1),
                                                    from = c(movable_slope_lim_lower, movable_slope_lim_upper)))))

  # Compute flow widths to 4-neighbors.
  # Indices: 1 top, 2 left, 3 right, 4 bottom (as in Gruber, 2007, Fig. 1).
  conv <- pi / 180
  grid_cell_size <- xres(aspect)
  flow_width <- list() # This is called L_NB by Gruber (2007).
  flow_width[[1]] <- cos(conv * aspect) * grid_cell_size
  flow_width[[2]] <- -sin(conv * aspect) * grid_cell_size
  flow_width[[3]] <- -flow_width[[2]]
  flow_width[[4]] <- -flow_width[[1]]
  
  # Compute elevation difference to 4-neighbors.
  # Same indexing as above.
  dz <- list()
  dz[[1]] <- setValues(elevation, c(rep(NA, grid_ncol), elevation[2:grid_nrow,] - elevation[1:(grid_nrow - 1),]))
  dz[[2]] <- dz[[1]]
  dz[[2]][,1:grid_ncol] <- c(rep(NA, grid_nrow), elevation[,2:grid_ncol] - elevation[,1:(grid_ncol - 1)]) # We cannot use setValues() here because we compute by column and setValues sets by row.
  dz[[3]] <- dz[[1]]
  dz[[3]][,1:grid_ncol] <- c(elevation[,1:(grid_ncol - 1)] - elevation[,2:grid_ncol], rep(NA, grid_nrow)) # We cannot use setValues() here because we compute by column and setValues sets by row.
  dz[[4]] <- setValues(elevation, c(elevation[1:(grid_nrow - 1),] - elevation[2:grid_nrow,], rep(NA, grid_ncol)))
  
  # Compute draining fractions (first the raw coefficient,
  # then normalized) to 4-neighbors (Eqs. 7 and 9 in Gruber, 2007).
  draining_coeff <- list() # This is called C_NB by Gruber (2007).
  for (dir_id in 1:4) {
    draining_coeff[[dir_id]] <- flow_width[[dir_id]] * (dz[[dir_id]] > 0) * (flow_width[[dir_id]] > 0)
  }
  draining_coeff_sum <- draining_coeff[[1]] + draining_coeff[[2]] + draining_coeff[[3]] + draining_coeff[[4]]
  
  # Despite all the DEM pre-processing, there can be cells
  # with Csum = 0, i.e. no drainage possible.
  # These arise when the computed flow width grids do not match
  # the computed dz grids (flow width is computed from aspect, which is
  # determined via curve fitting since a plane is determined
  # by 3 points, so on a regular grid with 4 neighbors aspect is not
  # well defined!). In those cases, the direction for which dz<i> > 0
  # does not match the direction for which L<i> > 0, so that no drainage
  # would be possible from those cells. These would behave as infinite
  # sinks, stealing mass as it enters but cannot exit.
  # To solve this, we take those (hopefully rare) cells and
  # we force drainage towards the most likely direction
  # (i.e., of all the directions with dz<i> > 0, the direction
  # for which flow width is least negative).
  residual_sink_cell_ids <- which(getValues(draining_coeff_sum) == 0)
  residual_sinks_n <- length(residual_sink_cell_ids)
  cat("Residual sinks detected:", residual_sinks_n)
  
  if (residual_sinks_n) {
    for (residual_sink_id in 1:residual_sinks_n) {
      residual_sink_cell_id <- residual_sink_cell_ids[residual_sink_id]
      
      # Get dz and flow widths of problem cell (two 4-vectors).
      # There won't be corresponding elements which are both positive
      # (else the cell would not be problematic).
      # We want the element with least negative flow width,
      # and positive dz.
      cell_orig_ids <- 1:4
      cell_dzs <- sapply(dz, `[`, residual_sink_cell_id)
      cell_flow_widths <- sapply(flow_width, `[`, residual_sink_cell_id)
      
      cell_df <- data.frame(cell_orig_ids, cell_dzs, cell_flow_widths)
      cell_df_downslope <- cell_df[which(cell_df$cell_dzs >= 0),]
      cell_df_id_sel <- which.max(cell_df_downslope$cell_flow_widths)
      
      dir_id <- cell_df_downslope$cell_orig_ids[cell_df_id_sel]
      
      # Send all snow to the selected cell.
      draining_coeff[[dir_id]][residual_sink_cell_id] <- 1
      draining_coeff_sum[residual_sink_cell_id] <- 1
    }
  }
  
  # Compute normalized draining fractions for the 4 directions.
  draining_fraction <- list()
  for (dir_id in 1:4) {
    draining_fraction[[dir_id]] <- draining_coeff[[dir_id]] / draining_coeff_sum
  }

  
  # Compute max deposition (Eq. 10 in Gruber, 2007); [kg m-2].
  deposition_max <- (1 - slope / deposition_slope_lim) * deposition_mass_lim * (slope < deposition_slope_lim)

  # Compute indices for loop from highest to lowest domain cell.
  # Cells are indexed by row (i.e. first all the first row, then
  # all the second row; the last cell is the bottom-right corner).
  # We exclude all cells along the border since their drainage is problematic.
  # Of course the glacier should not reach the grid border.
  elevation_sorted_ids_raw <- sort(elevation[], decreasing = TRUE, index.return = TRUE)[[2]]
  elevation_ids_border <- unique(c(cellFromRow(elevation, c(1, grid_nrow)), cellFromCol(elevation, c(1, grid_ncol))))
  elevation_sorted_ids <- setdiff(elevation_sorted_ids_raw, elevation_ids_border)
  
  # ALL THE CODE ABOVE THIS LINE SHOULD GO INTO A PREPROCESSING FUNCTION CALLED ONLY ONCE AT THE BEGINNING OF THE WHOLE RUN
  # (or at the beginning of each year if we change DEMs!!)
  # Then the computed grids can be used several times to compute avalanches.
  mass_initial <- setValues(elevation, 1000) # [kg m-2] FOR TESTING, will be replaced by the actual snow cover pre-avalanches.
  
  deposition <- setValues(elevation, 0.0)
  mass_movable <- mass_initial * movable_frac
  mass_fixed <- mass_initial * (1 - movable_frac) # Mass which stays in place no matter what.
  
  # The snow transport loop is implemented in C++ for performance (about 5000 times faster than pure R).
  # The legacy R implementation is in file "func_snow_transport_gruber_legacy.R.old",
  # to use it (you better don't do that) you should replace the lines below with all
  # the lines from that file, and be prepared to wait a long time before the loop is complete.
  deposition <- setValues(deposition, transport_deposit_mass(elevation_sorted_ids,
                                                             grid_ncol,
                                                             getValues(deposition),
                                                             getValues(mass_movable),
                                                             getValues(deposition_max),
                                                             getValues(draining_fraction[[1]]),
                                                             getValues(draining_fraction[[2]]),
                                                             getValues(draining_fraction[[3]]),
                                                             getValues(draining_fraction[[4]])))
  
  mass_final <- mass_fixed + deposition

  return(mass_final)
  
}