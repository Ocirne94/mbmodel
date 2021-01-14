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
#                 list[[1]] is the raster related to the up direction.                            #
################################################################################################### 


# This function shall be called as follows:
  # avalanche_grids <- sapply(avalanche_grids, `[[`, grid_id)  # To select all the appropriate grids.
  # mass_final <- func_avalanche(avalanche_grids, mass_initial, mult)
# The deposition_max_multiplier can be used to enable consistent modelling on different input grids
# (e.g., the normalized snow distribution grid, which has cell values close to 1; an actual snow cover grid,
# with cell values of maybe 1000 (kg m-2); and a seasonal sum (if we choose to have a single avalanche over
# a whole year), with cell values of maybe 5000 (kg m-2)).
# This addresses the discussion of Section 4.4 in Gruber (2007).
func_avalanche <- function(avalanche_grids, mass_initial, deposition_max_multiplier = 1.0) {
  
  deposition <- setValues(avalanche_grids$elevation_proc, 0.0)
  mass_movable <- mass_initial * avalanche_grids$movable_frac
  mass_fixed <- mass_initial - mass_movable  # Mass which stays in place no matter what.
  
  # The snow transport loop is implemented in C++ for performance (about 5000 times faster than pure R).
  # The legacy R implementation is in file "func_avalanche_gruber_legacy.R.old" (not updated after move to list() structures).
  deposition <- setValues(deposition, transport_deposit_mass(avalanche_grids$elevation_sorted_ids,
                                                             run_params$grid_ncol,
                                                             getValues(deposition),
                                                             getValues(mass_movable),
                                                             getValues(avalanche_grids$deposition_max) * deposition_max_multiplier,
                                                             getValues(avalanche_grids$draining_fraction[[1]]),
                                                             getValues(avalanche_grids$draining_fraction[[2]]),
                                                             getValues(avalanche_grids$draining_fraction[[3]]),
                                                             getValues(avalanche_grids$draining_fraction[[4]])))
  
  mass_final <- mass_fixed + deposition

  return(mass_final)
  
}