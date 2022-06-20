# Optimal_Stein_GoF
Optimal design of Stein-type GoF-tests relying on asymptotic computations.

Provides the R codes used for the simulations and data analyses in the article "Optimal Stein-type Goodness-of-Fit Tests for Count Data", which is a joint work by Christian H. Weiß, Pedro Puig, and Boris Aleksandrov.

The subsequent codes have been written by Christian H. Weiß.

E-mail: weissc@hsu-hh.de

The codes have been evaluated and tested under

R version 4.1.0 (2021-05-18)

Platform: x86_64-w64-mingw32/x64 (64-bit)

Running under: Windows 10 x64 (build 19044)

The R-package "dgof" is of version dgof_1.2.



## Files in this repository:
Note that all R codes can be run independently of each other.


### AsymptPowerSteinChen.r:
Refers to Stein-Chen statistic for Poi-null.
Code for computing the asymptotic power curves of Section 3 by using Theorem 3.1, and code for simulation results of Section 5 and Supplement S1.

### AsymptPowerBinStein.r:
Refers to Stein statistic for Bin-null.
Code for simulation results of Section 5 and Supplement S1, which are computed by using Theorem 3.1.

### AsymptPowerNBStein.r:
Refers to Stein statistic for NB-null.
Code for computing asymptotic power curves according to Appendix A.5.

### Example_Dicentrics_OptimalPoiStein.r:
Code for analyzing the unbounded dicentrics counts of Section 6.1.

### Example_Foci_OptimalPoiStein.r:
Code for analyzing the unbounded foci counts of Section 6.1.

### Example_DMFT_OptimalBinStein.r:
Code for analyzing the bounded DMFT counts of Section 6.2.
Requires to load data from file "dmft_bohning99.txt".

### dmft_bohning99.txt:
DMFT data of Böhning et al. (1999). These data are also offered by the R-package "flexmix: Flexible Mixture Modeling" on https://CRAN.R-project.org/package=flexmix, use data("dmft").

Böhning D, Dietz E, Schlattmann P, Mendonca L, Kirchner U. The zero-inflated Poisson model and the decayed, missing and filled teeth index in dental epidemiology. J Royal Stat Soc A. 1999;162(2):195–209.
