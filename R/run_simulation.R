#' @title Simulation flow
#' @description
#' Runs the simulation.
#' @param conditions a data frame with simulation conditions, typically
#' constructed using [prepare_conditions]
#' @param nIterPerCond a positive integer - number of iterations to be run for
#' each condition
#' @param suffix optionally a string - suffix that will be added to the name
#' of a file storing simulation results (that will be saved to the disk)
#' @param modelSpecsOnly optionally a logical flag indicating that only list
#' of model specifications should be returned (for each condition), without
#' actually estimating the models
#' @details
#' For description of estimated models see [prepare_models]. For references to
#' additional details see [run_iteration].
#' @returns (invisibly) a list four of data frames:
#' \describe{
#'   \item{conditions}{a data frame that was passed to the function using the
#'                     `conditions` argument}
#'   \item{modelSummaries}{basic model summary statistics - see [get_model_summary]}
#'   \item{structPars}{estimates of the structural part parameters -
#'                     see [get_model_pars], the additional column `gen` stores
#'                     the parameter value used in the data generating model}
#'   \item{measurePars}{estimates of the measurement part parameters -
#'                      see [get_model_pars], the additional columns `gen` and
#'                      `genLambda` store values of the parameters used in the
#'                      data generating model}
#' }
#' Moreover, after completing each simulation condition-iteration function will
#' save the data frames listed above to the file
#' "latentInteractions_results\[suffix\].RData".
#' @examples
#' \dontrun{
#' str(conditions)
#' set.seed(12345)
#' run_simulation(conditions[1:2, ], nIterPerCond = 1L)
#' }
#' @export
run_simulation <- function(conditions, nIterPerCond, suffix = "",
                           modelSpecsOnly = FALSE) {
  check_conditions(conditions)
  stopifnot(is.numeric(nIterPerCond), length(nIterPerCond) == 1,
            as.integer(nIterPerCond) == nIterPerCond, nIterPerCond > 0,
            is.character(suffix), length(suffix) == 1L, !anyNA(suffix),
            is.logical(modelSpecsOnly), length(modelSpecsOnly) == 1L,
            modelSpecsOnly %in% c(FALSE, TRUE))
  defaultOptWarn <- options()$warn
  options(warn = 1)
  on.exit(options(warn = defaultOptWarn))
  modelSummaries <- structPars <- measurePars <- data.frame()
  if (modelSpecsOnly) {
    modelSpecifications <- vector(mode = "list",
                                  length = nrow(conditions))
    if (nIterPerCond != 1L) message("With `modelSpecsOnly = TRUE` value of the `nIterPerCond` was automatically set to 1.")
    nIterPerCond <- 1L
  }
  for (i in seq_len(nIterPerCond)) {
    for (j in seq_len(nrow(conditions))) {
      cat("\n###########################################################################\n Simulation iteration ",
          i, " (out of ", nIterPerCond,"),\n condition number ", j, " (out of ",
          nrow(conditions),")\n\n",
          sep = "")
      print(as.data.frame(conditions[j, ]), row.names = FALSE)
      cat("###########################################################################\n\n")
      startTime <- Sys.time()
      resultsIter <- run_iteration(conditions[j, ],
                                   modelSpecsOnly = modelSpecsOnly)
      if (modelSpecsOnly) {
        modelSpecifications[[j]] <- resultsIter
        next
      }
      cat("\n  ", lasted(startTime), "\n", sep = "")
      modelSummaries <- dplyr::bind_rows(modelSummaries,
                                         cbind(i = rep(i, nrow(resultsIter$modelSummaries)),
                                               cond = rep(j, nrow(resultsIter$modelSummaries)),
                                               resultsIter$modelSummaries))
      structPars <- dplyr::bind_rows(structPars,
                                     cbind(i = rep(i, nrow(resultsIter$structPars)),
                                           cond = rep(j, nrow(resultsIter$structPars)),
                                           resultsIter$structPars))
      measurePars <- dplyr::bind_rows(measurePars,
                                      cbind(i = rep(i, nrow(resultsIter$measurePars)),
                                            cond = rep(j, nrow(resultsIter$measurePars)),
                                            resultsIter$measurePars))
      save(conditions, modelSummaries, structPars, measurePars,
           file = paste0("latentInteractions", suffix, ".RData"))
    }
  }
  if (modelSpecsOnly) return(structure(modelSpecifications,
                                       conditions = conditions))
  invisible(list(conditions = conditions,
                 modelSummaries = modelSummaries,
                 structPars = structPars,
                 measurePars = measurePars))
}
