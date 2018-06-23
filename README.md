# check_mboost: Diagnostic function for mboost models

* check_mboost.R: script for the actual diagnostic function producing results
* plot.check_mboost.R: script for plotting function(s)
* diagnostic_functions.R: script containing all provided diagnostic functions

TODOs:

* [ ] write plot functions
* [ ] complete check_mboost function
* [ ] test functions
* [ ] hat values calculation more efficient (step 1: calculate first, then criteria, step 2: efficient getUps)
* [ ] simulation study -> KLD, special residuals, KLD Binomial / Poisson
* [ ] write documentary
* [ ] pull request to mboost
* [ ] check whether residual matrix / hat matrix stuff works, especially for losses other than L2-loss
