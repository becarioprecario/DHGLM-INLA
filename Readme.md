# Introduction

This repository contains the `R` code for paper *Fitting Double Hierarchical Models with the Integrated Nested Laplace
 Approximation* by Morales-Otero et al. (2022) ([arXiv prep-rpint]()). The diffferent folders conatin the following examples:

* `Poisson`: Poisson model with random effects and a hierarchical structure on the precision of the random effects.

* `NegBin`: Negative binomial model with  a hierarchical structure on the precision parameter.

* `Gaussian5`: Multilevel Gaussian model with random effects for 5 groups and a  hierarchical structure (i.e., mixed-effects model) on the group precisions.

* `Colombia_sim`:

* `sleepstudy`:

* `slepstudy_rcoef`:


All folders contain the following `R` files:

* `01_data.R`: Script to collect the dataset and save it in `data.RData`.

* `02_AMIS.R`: Script to run the DHGLM model with AMIS with INLA and save the output.

* `03_JAGS.R`: Script to run the DHGLM model with JAGS and save the ou
tput.

* `04_summary.R`: Script to produce some summary results and plots.

* `05_output_paper.R`: Script to produce some summary results and the plots included in the paper.


