#' @title Model estimation
#' @description
#' Prepares a list of model estimating functions along with their arguments.
#' @param observed a matrix with values of observed indicators (excluding
#' product indicators for interaction terms)
#' @param observedMeanC a matrix with values of observed indicators including
#' product indicators for interaction terms that were created using *grand mean
#' centering* approach, typically returned by [create_interaction_indicators]
#' called with `center = "mean"`
#' @param observedResidC a matrix with values of observed indicators including
#' product indicators for interaction terms that were created using *residual
#' centering* approach, typically returned by [create_interaction_indicators]
#' called with `center = "residual"`
#' @param interactionsMapping a matrix describing which exogenous latent
#' variables construct interaction terms, typically the `mapping` attribute
#' of an object returned by the `generate_latent`
#' @details
#' At the moment function estimates models:
#' \describe{
#'   \item{unconstrMC.ML}{Product indicators, grand mean centered - unconstrained approach estimated by ML using *lavaan*}
#'   \item{unconstrMC.TSLS}{Product indicators, grand mean centered - unconstrained approach estimated by IVSEM using *MIIVsem*}
#'   \item{unconstrRC.ML}{Product indicators, residual centered - unconstrained approach estimated by ML using *lavaan*}
#'   \item{unconstrRC.TSLS}{Product indicators, residual centered - unconstrained approach estimated by IVSEM using *MIIVsem*}
#'   \item{unconstrRC.SAM}{Product indicators, residual centered - unconstrained approach estimated by Structural After Measurement (local) using *lavaan*}
#'   \item{LMS.ML}{Latent Moderated Structural implemented in *Mplus* - ML (EM) estimation}
#'   \item{LMS.Bayes}{Latent Moderated Structural implemented in *Mplus* - Bayesian estimation}
#'   \item{PLS}{Partial Least Squares-SEM implemented via *seminr* package}
#' }
#' @returns a list of lists: in each element the first element is a model
#' estimating function and other elements are arguments that should be passed
#' to it
#' @export
prepare_models <- function(observed, observedMeanC, observedResidC,
                           interactionsMapping) {
  syntaxUnconstrMeanC <-
    prepare_unconstrained_syntax(mapping = attributes(observedMeanC)$mapping,
                                 interactionsMapping = interactionsMapping)
  syntaxUnconstrResidC <-
    prepare_unconstrained_syntax(mapping = attributes(observedResidC)$mapping,
                                 interactionsMapping = NULL)
  syntaxLMS <- prepare_lms_syntax(mapping = attributes(observed)$mapping,
                                  interactionsMapping = interactionsMapping)
  # for some strange reasons lavaan::sem don't want to be called through do.call()
  # but lavaan::lavaan works well this way
  return(list(unconstrMC.ML = list(fun = lavaan::lavaan,
                                   model = syntaxUnconstrMeanC,
                                   data = observedMeanC,
                                   int.ov.free = TRUE, int.lv.free = FALSE,
                                   auto.fix.first = TRUE, auto.fix.single = TRUE,
                                   auto.var = TRUE, auto.cov.lv.x = TRUE,
                                   auto.efa = TRUE, auto.th = TRUE,
                                   auto.delta = TRUE, auto.cov.y = TRUE),
              unconstrMC.TSLS = list(fun = MIIVsem::miive,
                                     model = syntaxUnconstrMeanC,
                                     data = observedMeanC,
                                     estimator = "2SLS",
                                     se = "standard",
                                     var.cov = TRUE),
              # for some strange reasons lavaan::sem don't want to be called through do.call()
              # but lavaan::lavaan works well this way
              unconstrRC.ML = list(fun = lavaan::lavaan,
                                   model = syntaxUnconstrResidC,
                                   data = observedResidC,
                                   int.ov.free = TRUE, int.lv.free = FALSE,
                                   auto.fix.first = TRUE, auto.fix.single = TRUE,
                                   auto.var = TRUE, auto.cov.lv.x = TRUE,
                                   auto.efa = TRUE, auto.th = TRUE,
                                   auto.delta = TRUE, auto.cov.y = TRUE),
              unconstrRC.TSLS = list(fun = MIIVsem::miive,
                                     model = syntaxUnconstrResidC,
                                     data = observedResidC,
                                     estimator = "2SLS",
                                     se = "standard",
                                     var.cov = TRUE),
              unconstrRC.SAM = list(fun = lavaan::sam,
                                    model = syntaxUnconstrResidC,
                                    data = observedResidC,
                                    sam.method = "local",
                                    se = "twostep"),
              LMS.ML = list(fun = MplusAutomation::mplusModeler,
                            object = MplusAutomation::mplusObject(
                              ANALYSIS = "TYPE IS RANDOM;
ESTIMATOR IS MLR;
ALGORITHM IS INTEGRATION EM;
INTEGRATION IS STANDARD(15);
PROCESSORS ARE 4;",
                              OUTPUT = "STANDARDIZED (STDYX) TECH8;",
                              MODEL = syntaxLMS,
                              rdata = as.data.frame(observed),
                              autov = TRUE),
                            run = 1L,
                            modelout = "lms_em.inp"),
              LMS.Bayes = list(fun = MplusAutomation::mplusModeler,
                               object = MplusAutomation::mplusObject(
                                 ANALYSIS = "TYPE IS RANDOM;
ESTIMATOR IS BAYES;
CHAINS ARE 2;
PROCESSORS ARE 4;",
                                 OUTPUT = "STANDARDIZED (STDYX) TECH8;",
                                 MODEL = syntaxLMS,
                                 rdata = as.data.frame(observed),
                                 autov = TRUE),
                               run = 1L,
                               modelout = "lms_bayes.inp"),
              PLS = list(fun =
                           function(...) {
                             seminr::bootstrap_model(seminr::estimate_pls(...),
                                                     nboot = 500)
                           },
                         data = observed,
                         measurement_model =
                           prepare_pls_measurement_syntax(attributes(observed)$mapping,
                                                          interactionsMapping),
                         structural_model =
                           prepare_pls_structural_syntax(attributes(observed)$mapping,
                                                         interactionsMapping),
                         maxIt = 300)))
}
#' @title Model estimation
#' @description
#' Prepares *lavaan*/*miivsem* model syntax for the *unconstrained* approach
#' to latent interactions model specification.
#' @param mapping a matrix describing which observed indicators are loaded by
#' which latent variables, typically the `mapping` attribute of an object
#' returned by the [create_interaction_indicators]
#' @inheritParams prepare_models
#' @returns a string with model syntax
prepare_unconstrained_syntax <- function(mapping, interactionsMapping) {
  if (!is.null(interactionsMapping)) {
    stopifnot("Higher-order latent interactions are not supported" =
                all(rowSums(interactionsMapping) <= 2))
  }
  measurement <-
    paste0(colnames(mapping),
           " =~ ",
           # ifelse(grepl("^xi[[:digit:]]*$", colnames(mapping)),
           #        " =~ NA*", " =~ "),
           apply(mapping, 2,
                 function(x, rownames) paste(rownames[x != 0], collapse = " + "),
                 rownames = rownames(mapping)))
  structural <- paste0(grep("^y", colnames(mapping), value = TRUE), " ~ ",
                       paste(grep("^x", colnames(mapping), value = TRUE),
                             collapse = " + "))
  if (is.null(interactionsMapping)) {
    lV <- grep("^x[[:digit:]]*$", colnames(mapping), value = TRUE)
    lI <- grep("^xi[[:digit:]]*$", colnames(mapping), value = TRUE)
    covariances <- paste0(rep(lI, length(lV)),
                          " ~~ 0*", rep(lV, each = length(lI)))
    variances <- expectedValues <- constraints <-
      vector(mode = "character", length = 0L)
  } else {
    covariances <- apply(interactionsMapping, 1,
                         function(x, names) utils::combn(names[x != 0], 2),
                         names = colnames(interactionsMapping))
    variances <- constraints <- vector(mode = "character", length = 0L)
    # variances <-
    #   c(colnames(interactionsMapping), rownames(interactionsMapping))
    # variancesLabels <- paste0("v", seq_along(variances))
    # names(variancesLabels) <- variances
    # constraints <- paste0(variancesLabels[-seq_len(ncol(interactionsMapping))],
    #                       " == ", variancesLabels[covariances[1, ]], "*",
    #                       variancesLabels[covariances[2, ]],
    #                       "+c", seq_len(ncol(covariances)), "^2")
    covariances <- paste0(covariances[1, ], " ~~ c", seq_len(ncol(covariances)),
                          "*", covariances[2, ])
    # variances <- paste0(variances, " ~~ ", variancesLabels, "*", variances)
    expectedValues <- paste0(rownames(interactionsMapping), " ~ c",
                             seq_len(nrow(interactionsMapping)), "*1")
  }
  return(paste(c(measurement, variances, covariances, constraints,
                 expectedValues, structural),
               collapse = "\n"))
}
#' @title Model estimation
#' @description
#' Prepares *Mplus* model syntax for the *Latent Moderated Structural* (LMS)
#' approach to latent interactions model specification.
#' @inheritParams prepare_unconstrained_syntax
#' @returns a string with model syntax
prepare_lms_syntax <- function(mapping, interactionsMapping) {
  if (!is.null(interactionsMapping)) {
    stopifnot("Higher-order latent interactions are not supported" =
                all(rowSums(interactionsMapping) <= 2))
  }
  measurement <- paste0(colnames(mapping), " BY ",
                        apply(mapping, 2,
                              function(x, rownames) paste(rownames[x != 0],
                                                          collapse = " "),
                              rownames = rownames(mapping)), ";")
  interactions <- paste0(rownames(interactionsMapping), " | ",
                         apply(interactionsMapping, 1,
                               function(x, colnames) paste(colnames[x != 0],
                                                           collapse = " XWITH "),
                               colnames = colnames(interactionsMapping)), ";")
  structural <- paste0(grep("^y", colnames(mapping), value = TRUE), " ON ",
                       paste(c(grep("^[^y][[:digit:]]*$", colnames(mapping),
                                    value = TRUE),
                                    rownames(interactionsMapping)),
                             collapse = " "), ";")
  return(paste(c(measurement, interactions, structural), collapse = "\n"))
}
#' @title Model estimation
#' @description
#' Prepares *seminr* description of the measurement part of a PLS model.
#' @inheritParams prepare_unconstrained_syntax
#' @returns an object of classes *measurement_model* and *seminr_model* returned
#' by [seminr::constructs]
prepare_pls_measurement_syntax <- function(mapping, interactionsMapping = NULL) {
  mainEffects <- mapply(seminr::composite,
                        construct_name = colnames(mapping),
                        item_names = apply(mapping, 2,
                                           function(x, rownames) rownames[x != 0],
                                           rownames = rownames(mapping),
                                           simplify = FALSE),
                        MoreArgs = list(weights = seminr::correlation_weights),
                        SIMPLIFY = FALSE)
  if (!is.null(interactionsMapping)) {
    stopifnot("Higher-order latent interactions are not supported" =
                all(rowSums(interactionsMapping) <= 2))
    interactions <- apply(interactionsMapping, 1,
                          function(x, colnames) colnames[x != 0],
                          colnames = colnames(interactionsMapping))
    interactions <- mapply(seminr::interaction_term,
                           iv = interactions[1, ],
                           moderator = interactions[2, ],
                           MoreArgs = list(method = seminr::two_stage),
                           SIMPLIFY = FALSE)
  } else {
    interactions <- NULL
  }
  return(do.call(seminr::constructs, c(mainEffects, interactions)))
}
#' @title Model estimation
#' @description
#' Prepares *seminr* description of the structural part of a PLS model.
#' @inheritParams prepare_unconstrained_syntax
#' @returns an object of classes *structural_model* and *seminr_model* returned
#' by [seminr::relationships]
prepare_pls_structural_syntax <- function(mapping, interactionsMapping = NULL) {
  return(seminr::relationships(seminr::paths(
    from = c(grep("^x[[:digit:]]*$", colnames(mapping), value = TRUE),
             apply(interactionsMapping, 1,
                   function(x, colnames) paste(colnames[x != 0], collapse = "*"),
                   colnames = colnames(interactionsMapping))),
    to = grep("^y[[:digit:]]*$", colnames(mapping), value = TRUE))))
}
