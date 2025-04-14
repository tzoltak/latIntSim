# latIntSim

Simulation studies on SEM with latent interactions (moderations)

## Installation/update

Currently package is not submitted to CRAN. Nevertheless you can easily install it using the remotes package:

1.  First you need to install the remotes package, if you don't have it already:
    `r install.packages('remotes')`
2.  Next you can install *logLime*: `remotes::install_github('tzoltak/logLime')`

## Usage

### The simplest usecase

The main function of the package is `run_simulation()`. You need to provide it with a data frame defining a set of simulation conditions (each condition described using a single row) and an integer specyfing the number of replications that sholud be repeted for each of the conditions.

F.g. code below uses first to rows from the built-in package dataset `conditions` as a simulation grid and run only one iteration for each of these two conditions:

```
library(latIntSim)
set.seed(12345) # always remember to set some seed before runing a simulation!
results <- run_simulation(conditions[1:2, ], nIterPerCond = 1L)
```

Function returns (invisibly) a list four of data frames:
1.  `conditions` - simulation grid that was passed to the function,
2.  `modelSummaries` - basic model summary statistics - see `?get_model_summary` for more information,
3.  `structPars` - estimates of the structural part parameters - see `?get_model_pars` for more information (the additional column `gen` stores the parameter value used in the data generating model),
4.  `measurePars` - estimates of the measurement part parameters - see `?get_model_pars` for more information (the additional column `gen` stores the parameter value used in the data generating model).

Moreover, after completing each simulation condition-iteration function will save the data frames listed above to the file "latentInteractions_results.RData".

### Defininig simulation conditions

A data frame containg a simulation grid (set of simulation conditions) can be prepared using a package built-in dataset `conditions` as a template. All the variables included in this data frame are required by `run_simulation()` and also it will be unable to make use of any other variables. Consequently, you need to specify:

-   `nObs` (an integer) - number of observations,
-   `nExogLVs` (an integer) - number of exogenous latent variables),
-   `corExogLVs` (a number) - value of pairwise correlations between exogenous latent variables),
-   `nLIs` (an integer) - number of latent (first-order) interactions,
-   `nIndic` (an integer) - number of observed indicators for each exogenous latent variable and for the dependent latent variable,
-   `lambda` (a number) - value of standardized factor loadings (the same for all observed indicators of all latent variables),
-   `partCorMain` (a number) - value of partial correlation between (each) exogenous latent variable and the dependent latent variable, given all the other exogenous latent variables (and latent interactions),
-   `partCorInt` (a number) - value of partial correlation between (each) latent interaction variable and the dependent latent variable, given all the other latent interaction variables (and exogenous latent variables).

If you want to set your simulation design as *fully crossed* given some sets of values for each of the aforementioned parameters, you may find useful an auxliary function `prepare_conditions()` that was included in the package (see `?prepare_conditions` for more information).

### Runing simulations in parallel

At the moment the only parallelisation strategy handled by the package is *parallelization by hand*, i.e. by manually deploying several processes (possibly on different machines) by the user (do not forget to set different seeds at the begining of each process!).

While I acknowledge that this may not be seen as a proper parallelisation (I won't argue), in my personal opinion it works good enough. With methods like QML and Bayesian estimation of interaction models often require large memory resources, so one will not be able to run many processes in parallel on a single machine anyway. Moreover, manually deploying several processes is not as big problem.

There are two things one need to remember while deploying parallel processess: 1) run each process in a separate folder, 2) set different random seed for each process.

To support this kind of *parallelization* there are two features included in the package:

1.  You can provide an additional, third argument `suffix` to `run_simulation()`, that will modify the name of the .RData file in which results will be stored.
2.  There is a function `combine_results()` enabling to gather results saved in several different .RData file into one set of data frames.
