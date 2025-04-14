#' @title Simulation flow
#' @description
#' Performs a single iteration of the simulation.
#' @param condition a one-row data frame with additional parameters specifying
#' data generation process that will be passed to data generating functions
#' ([generate_latent], [generate_loadings] and [generate_observed])
#' @param modelSpecsOnly optionally a logical flag indicating that only list
#' of model specifications should be returned, without actually estimating the
#' models
#' @returns By default or with `modelSpecsOnly = FALSE` - a list of three data
#' frames:
#' \describe{
#'   \item{modelSummaries}{basic model summary statistics - see [get_model_summary]}
#'   \item{structPars}{estimates of the structural part parameters -
#'                     see [get_model_pars], the additional column `gen` stores
#'                     the parameter value used in the data generating model}
#'   \item{measurePars}{estimates of the measurement part parameters -
#'                      see [get_model_pars], the additional columns `gen` and
#'                      `genLambda` store values of the parameters used in the
#'                      data generating model}
#' }
#' With `modelSpecsOnly = TRUE` a list returned by [prepare_models].
#' @seealso [prepare_models], [get_model_summary], [get_model_pars]
run_iteration <- function(condition, modelSpecsOnly = FALSE) {
  stopifnot(is.data.frame(condition), nrow(condition) == 1L,
            is.logical(modelSpecsOnly), length(modelSpecsOnly) == 1L,
            modelSpecsOnly %in% c(FALSE, TRUE))
  condition <- as.list(condition)

  latent <- generate_latent(nObs = condition$nObs, nLV = condition$nExogLVs,
                            nLI = condition$nLIs, r = condition$corExogLVs,
                            partR = condition$partCorInt,
                            partRI = condition$partCorInt)
  nIndic <- rep(condition$nIndic, condition$nExogLVs + 1)
  names(nIndic) <- grep("^[xy][[:digit:]]*$", colnames(latent), value = TRUE)
  loadings <- generate_loadings(nIndic = nIndic, fun = condition$lambda)
  observed <- generate_observed(latent[, grep("^[xy][l[:digit:]]*$",
                                              colnames(latent))],
                                loadings, loadingsStandardized = TRUE,
                                latentVariances =
                                  diag(attributes(latent)$cov)[grep("^[xy][l[:digit:]]*$",
                                                                    colnames(latent))])
  observedMeanC <- create_interaction_indicators(observed, loadings,
                                                 attributes(latent)$mapping,
                                                 center = "mean",
                                                 nonredundantOnly = TRUE)
  observedDoubleC <- create_interaction_indicators(observed, loadings,
                                                  attributes(latent)$mapping,
                                                  center = "double",
                                                  nonredundantOnly = TRUE)
  observedResidC <- create_interaction_indicators(observed, loadings,
                                                  attributes(latent)$mapping,
                                                  center = "residual",
                                                  nonredundantOnly = TRUE)
  models <- prepare_models(observed = observed,
                           observedMeanC = observedMeanC,
                           observedDoubleC = observedDoubleC,
                           observedResidC = observedResidC,
                           interactionsMapping = attributes(latent)$mapping)
  if (modelSpecsOnly) return(models)
  estimationTimes <- rep(NA_real_, length(models))
  for (i in seq_along(models)) {
    cat("Running model ", names(models)[i], sep = "")
    startTime <- Sys.time()
    models[[i]] <- suppressMessages(
      try(do.call(models[[i]]$fun, models[[i]][setdiff(names(models[[i]]),
                                                       "fun")])))
    estimationTimes[i] <- difftime(Sys.time(), startTime, units = "secs")
    cat(" ", lasted(startTime), "\n", sep = "")
    if (inherits(models[[i]], "mplusObject")) {
      unlink(models[[i]]$results$input$data$file)
      unlink(attributes(models[[i]]$results$summaries)$filename)
      unlink(sub("\\.out$", ".inp",
                 attributes(models[[i]]$results$summaries)$filename))
    }
  }
  modelSummaries <- dplyr::bind_rows(lapply(models, get_model_summary),
                                     .id = "model")
  modelSummaries$time <- estimationTimes
  modelPars <- lapply(models, get_model_pars)
  structPars = dplyr::bind_rows(lapply(modelPars,
                                       function(x) x$structural),
                                .id = "model")
  structPars <- merge(structPars,
                      data.frame(DV = "y",
                                 IV = names(attributes(latent)$slopes),
                                 gen = attributes(latent)$slopes))
  measurePars = dplyr::bind_rows(lapply(modelPars,
                                        function(x) x$measurement),
                                 .id = "model")
  measurePars <- merge(measurePars,
                       data.frame(LV = apply(loadings, 1,
                                             function(x, colnames) colnames[x != 0],
                                             colnames = colnames(loadings)),
                                  obsIndic = rownames(loadings),
                                  gen = rowSums(attributes(observed)$loadings),
                                  genLambda = rowSums(loadings),
                                  stringsAsFactors = FALSE))

  return(list(modelSummaries = modelSummaries,
              structPars = structPars,
              measurePars = measurePars))
}
