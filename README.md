# Glacier mass balance model for interpolation of point measurements

The idea is to run a mass balance model (accumulation and radiation-corrected degree-day melt), optimizing the model parameters towards the best fit with a set of point mass balance measurements.


## Model description

We model at daily scale, one year at a time (typically from September Y-1 to September Y).

First we compute an initial snow distribution for the year, from snowline elevation and snow elevation gradient, and a distribution grid (function of curvature and elevation, as well as a process-based avalanche model and a correction based on winter snow probes if available).

Specifically, the initial snow distribution can be expressed as:

D_probes_norm * reduce(avalanche(D_curv_ele * D_tsl_snowgrad))          (1)

where\
**D_probes_norm** is the normalized distribution from IDW interpolation of the winter snow probes (only if available),\
**reduce()** is a function which optionally brings the distribution closer to the mean, to reduce variability,\
**avalanche()** is a function to simulate an avalanche (Gruber, 2007), with actual cell-wise downstream mass transport,\
**D_curv_ele** is a normalized grid computed from surface curvature plus an elevation cutoff for accumulation,\
**D_tsl_snowgrad** is the distribution from snow line elevation and a snow height gradient with altitude.

By contrast, snow distribution following each precipitation event (during the actual model run) can be expressed as:

reduce(D_probes_norm) * D_curv_ele * solid_fraction * distribute_precgrad_cutoff(correct(precipitation)),            (2)\
while the avalanche is simulated at user-defined fixed dates of the year (typically end-of-winter and sometime in summer).

In (2),\
**reduce()** is as in (1), and can be used (with a custom reduction factor) in case snow probes are acquired on avalanche deposits, to avoid accounting twice for avalanche redistribution,\
**D_probes_norm** is as in (1),\
**D_curv_ele** is as in (1),\
**solid_fraction** is the solid fraction of each precipitation event, computed from temperature for each grid cell and each time-step,\
**distribute_precgrad_cutoff()** is a function which increases precipitation amount with altitude, up to a cutoff altitude,\
**correct()** is a function to correct for precipitation under-catch at the rain gauge, with the option to reduce the correction in summer (differential under-catch of snow and rain),\
**precipitation** is the rain gauge reading.
