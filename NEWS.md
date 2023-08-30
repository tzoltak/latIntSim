# latIntSim 0.2.0 (30.08.2023)

## New features

-   It is possible to merge simulation results from different .RData files using the new `combine_results()` function.

## Bug fixes

-   Function `run_iteration()` works on even if model estimation throws an error.
-   Function `lasted()`, responsible for formating estimation time reported on the console computed number of seconds improperly.
-   Function `get_model_pars()` used name *std.all* instead of *lambda* in the (empty) data frame with measurement model parameters for Mplus models that did not converge.
-   Function `get_model_pars()` did not returned interaction effects of the structural model.

# latIntSim 0.1.0.9000 (30.06.2023)

* First usable version of the package.
