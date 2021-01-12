//###################################################################################################
//# Author:         Enrico Mattea (@unifr.ch)                                                       #
//# Description:    this program models the distributed mass balance of a glacier at daily          #
//#                 resolution, optimizing model parameters towards the best fit with point         #
//#                 mass balance measurements.                                                      #
//#                 This file contains the computation of a grid of relative (normalized)           #
//#                 snow distribution from elevation, slope and curvature.                          #
//# Latest update:  2021.1.10                                                                       #
//###################################################################################################

#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector transport_deposit_mass(IntegerVector elevation_sorted_ids,
                                     unsigned int grid_ncol,
                                     NumericVector deposition,
                                     NumericVector mass_movable,
                                     NumericVector deposition_max,
                                     NumericVector draining_fraction1,
                                     NumericVector draining_fraction2,
                                     NumericVector draining_fraction3,
                                     NumericVector draining_fraction4
                                    ) {
    
    elevation_sorted_ids = elevation_sorted_ids - 1; // C++ indexing starts from 0, R indexing starts from 1!
    unsigned int n_cells = elevation_sorted_ids.size();
    
    // Iterate over the cells from the highest to the lowest,
    // moving mass downstream according to the geometries and thresholds.
    for (unsigned int cell_cur_id_id = 0; cell_cur_id_id < n_cells; cell_cur_id_id++) {
        
        // Make loop interruptible.
        if (cell_cur_id_id % 1000u == 0u) Rcpp::checkUserInterrupt();
    
        unsigned int cell_cur_id = elevation_sorted_ids[cell_cur_id_id];
        
        // 1 = top, 2 = left, 3 = right, 4 = bottom.
        // We are iterating only on the non-border cells,
        // so the adjacent ones are always easy to compute.
        unsigned int cell_adjacent_id1 = cell_cur_id - grid_ncol;
        unsigned int cell_adjacent_id2 = cell_cur_id - 1;
        unsigned int cell_adjacent_id3 = cell_cur_id + 1;
        unsigned int cell_adjacent_id4 = cell_cur_id + grid_ncol;
        
        deposition[cell_cur_id] = std::min(deposition_max[cell_cur_id], mass_movable[cell_cur_id]); // Eq. 1 in Gruber (2007).
        double flow1 = (mass_movable[cell_cur_id] - deposition[cell_cur_id]) * draining_fraction1[cell_cur_id]; // Eq. 2 in Gruber (2007).
        double flow2 = (mass_movable[cell_cur_id] - deposition[cell_cur_id]) * draining_fraction2[cell_cur_id];
        double flow3 = (mass_movable[cell_cur_id] - deposition[cell_cur_id]) * draining_fraction3[cell_cur_id];
        double flow4 = (mass_movable[cell_cur_id] - deposition[cell_cur_id]) * draining_fraction4[cell_cur_id];
        
        mass_movable[cell_adjacent_id1] = mass_movable[cell_adjacent_id1] + flow1;
        mass_movable[cell_adjacent_id2] = mass_movable[cell_adjacent_id2] + flow2;
        mass_movable[cell_adjacent_id3] = mass_movable[cell_adjacent_id3] + flow3;
        mass_movable[cell_adjacent_id4] = mass_movable[cell_adjacent_id4] + flow4;
        
    }
    
    return deposition;
}
