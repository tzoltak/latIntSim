# latIntSim 0.3.1 (13.04.2025)

## New features

-   It is possible to retrieve model specifications used within the simulation by running `run_simulation()` with the new argument `modelSpecsOnly=TRUE`.

## Documentation

-   Some information regarding the ways to use the package was included into README.md.
-   A vignette containing model syntax used in the simulation conditions included in the `conditions` built-in dataset was added.
-   LICENSE file updated inline with the current policy.

# latIntSim 0.3.0 (18.04.2024)

## New features

-   Changed SAM estimation from the approach using product indicators to the newly developed (and much better) Kronecker Product Method.
-   Changed PLS estimation to its consistent (under normality) version by switching to the *cSEM* package.
-   Implemented *double mean centering* approach.

## Bug fixes

-   Model specification in product indicators approach with mean (or *double mean*) centering and more than one latent interaction includes residual covariances of product indicators sharing one *first order* indicator.
-   Standard errors of standardized parameters are extracted (if available and otherwise computed) by `get_model_pars()` to get rid of the problem with division by 0 which (rarely) occurred while post-processing results from Mplus (which writes results in its output files using only three decimal places).
-   Function `get_model_summary()` checks *lavaan* results for negative variances (because of the model specification these will be always, if any, residual variances of observed indicators) and if it finds any, it writes a warning to the results.
-   Function `get_model_summary()` correctly deals with errors.
-   Function `check_conditions()` correctly specifies conditions for checking whether parameters are integers.

# latIntSim 0.2.0 (30.08.2023)

## New features

-   It is possible to merge simulation results from different .RData files using the new `combine_results()` function.

## Bug fixes

-   Function `run_iteration()` works on even if model estimation throws an error.
-   Function `lasted()`, responsible for formatting estimation time reported on the console computed number of seconds improperly.
-   Function `get_model_pars()` used name *std.all* instead of *lambda* in the (empty) data frame with measurement model parameters for Mplus models that did not converge.
-   Function `get_model_pars()` did not returned interaction effects of the structural model.

# latIntSim 0.1.0.9000 (30.06.2023)

* First usable version of the package.
