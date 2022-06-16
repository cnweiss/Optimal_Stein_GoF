# Optimal_Stein_GoF
Optimal design of Stein-type GoF-tests relying on asymptotic computations.

Provides the R codes used for the simulations and data analyses in the article "Optimal Stein-type Goodness-of-Fit Tests for Count Data", which is a joint work by Christian H. Weiß, Pedro Puig, and Boris Aleksandrov.

## Files in this repository:

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

### Example_DMFT_OptimalBinStein.r:
Code for analyzing the bounded DMFT counts of Section 6.2.
Requires to load data from file "dmft_bohning99.txt".

### dmft_bohning99.txt:
DMFT data of Böhning et al. (1999). Data also offered by R-package "flexmix: Flexible Mixture Modeling" on https://CRAN.R-project.org/package=flexmix, use data("dmft").
