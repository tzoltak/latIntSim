#' @title Model estimation
#' @description
#' Prepares a list of model estimating functions along with their arguments.
#' @param observed a matrix with values of observed indicators (excluding
#' product indicators for interaction terms)
#' @param observedMeanC a matrix with values of observed indicators including
#' product indicators for interaction terms that were created using *grand mean
#' centering* approach, typically returned by [create_interaction_indicators]
#' called with `center = "mean"`
#' @param observedDoubleC a matrix with values of observed indicators including
#' product indicators for interaction terms that were created using *double mean
#' centering* approach, typically returned by [create_interaction_indicators]
#' called with `center = "double"`
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
#'   \item{PLS}{Partial Least Squares-SEM implemented via *cSEM* package}
#' }
#' @returns a list of lists: in each element the first element is a model
#' estimating function and other elements are arguments that should be passed
#' to it
#' @export
prepare_models <- function(observed, observedMeanC, observedDoubleC,
                           observedResidC, interactionsMapping) {
  syntaxUnconstrMeanC <-
    prepare_unconstrained_syntax(mapping = attributes(observedMeanC)$mapping,
                                 interactionsMapping = interactionsMapping,
                                 centering = "mean")
  syntaxUnconstrDoubleC <-
    prepare_unconstrained_syntax(mapping = attributes(observedMeanC)$mapping,
                                 interactionsMapping = interactionsMapping,
                                 centering = "double")
  syntaxUnconstrResidC <-
    prepare_unconstrained_syntax(mapping = attributes(observedResidC)$mapping,
                                 interactionsMapping = NULL,
                                 centering = "resid")
  syntaxLav <- prepare_lav_syntax(mapping = attributes(observed)$mapping,
                                  interactionsMapping = interactionsMapping)
  syntaxLMS <- prepare_lms_syntax(mapping = attributes(observed)$mapping,
                                  interactionsMapping = interactionsMapping)
  # for some strange reasons lavaan::sem don't want to be called through do.call()
  # but lavaan::lavaan works well this way
  return(list(unconstrMC.ML = list(fun = lavaan::lavaan,
                                   model = syntaxUnconstrMeanC,
                                   data = observedMeanC,
                                   se = "robust.sem",
                                   int.ov.free = TRUE, int.lv.free = FALSE,
                                   auto.fix.first = TRUE, auto.fix.single = TRUE,
                                   auto.var = TRUE, auto.cov.lv.x = TRUE,
                                   auto.cov.y = TRUE),
              unconstrMC.TSLS = list(fun = MIIVsem::miive,
                                     model = syntaxUnconstrMeanC,
                                     data = observedMeanC,
                                     estimator = "2SLS",
                                     se = "standard",
                                     var.cov = TRUE),
              # for some strange reasons lavaan::sem don't want to be called through do.call()
              # but lavaan::lavaan works well this way
              unconstrDC.ML = list(fun = lavaan::lavaan,
                                   model = syntaxUnconstrDoubleC,
                                   data = observedDoubleC,
                                   se = "robust.sem",
                                   int.ov.free = TRUE, int.lv.free = FALSE,
                                   auto.fix.first = TRUE, auto.fix.single = TRUE,
                                   auto.var = TRUE, auto.cov.lv.x = TRUE,
                                   auto.cov.y = TRUE),
              unconstrDC.TSLS = list(fun = MIIVsem::miive,
                                     model = syntaxUnconstrDoubleC,
                                     data = observedDoubleC,
                                     estimator = "2SLS",
                                     se = "standard",
                                     var.cov = TRUE),
              # for some strange reasons lavaan::sem don't want to be called through do.call()
              # but lavaan::lavaan works well this way
              unconstrRC.ML = list(fun = lavaan::lavaan,
                                   model = syntaxUnconstrResidC,
                                   data = observedResidC,
                                   se = "robust.sem",
                                   int.ov.free = TRUE, int.lv.free = FALSE,
                                   auto.fix.first = TRUE, auto.fix.single = TRUE,
                                   auto.var = TRUE, auto.cov.lv.x = TRUE,
                                   auto.cov.y = TRUE),
              unconstrRC.TSLS = list(fun = MIIVsem::miive,
                                     model = syntaxUnconstrResidC,
                                     data = observedResidC,
                                     estimator = "2SLS",
                                     se = "standard",
                                     var.cov = TRUE),
              SAM = list(fun = lavaan::sam,
                         model = syntaxLav,
                         data = observed,
                         sam.method = "local",
                         se = "twostep"),
              LMS.Mplus = list(fun = MplusAutomation::mplusModeler,
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
              Bayes.Mplus = list(fun = MplusAutomation::mplusModeler,
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
              # PLS.seminr = list(fun =
              #                   function(...) {
              #                     seminr::bootstrap_model(seminr::estimate_pls(...),
              #                                             nboot = 500)
              #                   },
              #                   data = observed,
              #                   measurement_model =
              #                     prepare_pls_measurement_syntax(attributes(observed)$mapping,
              #                                                    interactionsMapping),
              #                   structural_model =
              #                     prepare_pls_structural_syntax(attributes(observed)$mapping,
              #                                                   interactionsMapping),
              #                   maxIt = 300),
              PLS.cSEM = list(fun =
                                function(...) {
                                  cSEM::resamplecSEMResults(cSEM::csem(...),
                                                            R = 499,
                                                            .resample_method = "bootstrap")
                                },
                              .data = observed,
                              .model = gsub(":", ".", syntaxLav),
                              .disattenuate = TRUE,
                              .iter_max = 100)))
}
#' @title Model estimation
#' @description
#' Prepares *lavaan*/*miivsem* model syntax for the *unconstrained* approach
#' to latent interactions model specification.
#' @param mapping a matrix describing which observed indicators are loaded by
#' which latent variables, typically the `mapping` attribute of an object
#' returned by the [create_interaction_indicators]
#' @inheritParams prepare_models
#' @param centering character vector describing method of centering the observed
#' inidicators: "mean" - *first-order* indicators are mean centered (but
#' product indicators are not), "double" - both *first-order* and product
#' indicators are centered, "resid" - product indicators are residuals from the
#' regression of products of *first-order* indicators on these *first-order*
#' indicators
#' @returns a string with model syntax
prepare_unconstrained_syntax <- function(mapping, interactionsMapping,
                                         centering = c("mean", "double", "resid")) {
  centering <- match.arg(centering)
  if (!is.null(interactionsMapping)) {
    stopifnot("Higher-order latent interactions are not supported" =
                all(rowSums(interactionsMapping) <= 2))
  }
  measurement <-
    paste0(colnames(mapping),
           " =~ ",
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
    expectedValues <- vector(mode = "character", length = 0L)
  } else {
    if (centering == "mean") {
      covariances <- apply(interactionsMapping, 1,
                           function(x, names) utils::combn(names[x != 0], 2),
                           names = colnames(interactionsMapping))
      covariances <- paste0(covariances[1, ], " ~~ c", seq_len(ncol(covariances)),
                            "*", covariances[2, ])
      expectedValues <- c(paste0(rownames(interactionsMapping), " ~ c",
                                 seq_len(nrow(interactionsMapping)), "*1"),
                          "y ~ NA*1")
    } else {
      covariances <- expectedValues <- vector(mode = "character", length = 0L)
    }
    if (any(colSums(interactionsMapping) > 1)) {
      for (i in which(colSums(interactionsMapping) == 2)) {
        interactIndic <- apply(mapping[, rownames(interactionsMapping)[interactionsMapping[, i] > 0]],
                               2,
                               function(x, names) {names(x)[x > 0]},
                               names = rownames(mapping),
                               simplify = FALSE)
        interactIndic <- lapply(interactIndic,
                                function(x) {
                                  data.frame(indic = x,
                                             order = sub("^.*_", "", x))})
        for (j in seq_along(interactIndic[-1])) {
          interactIndic[[1]] <- merge(interactIndic[[1]],
                                      interactIndic[[1 + j]],
                                      by = "order", all = TRUE)
        }
        interactIndic <- as.matrix(interactIndic[[1]][, -1])
        interactIndic <- interactIndic[rowSums(!is.na(interactIndic)) > 1, ]

        covariances <- paste0(
          c(covariances,
            unlist(apply(interactIndic, 1,
                         function(x) {
                           apply(utils::combn(x, m = 2), 2,
                                 paste, collapse = " ~~ ", simplify = FALSE)
                         }, simplify = FALSE))),
          collapse = "\n")
      }
    }
  }
  return(paste(c(measurement, covariances, expectedValues, structural),
               collapse = "\n"))
}
#' @title Model estimation
#' @description
#' Prepares *lavaan* model syntax for the *Kronecker/Tensor Product Method*
#' approach to latent interactions model specification.
#' @inheritParams prepare_unconstrained_syntax
#' @returns a string with model syntax
prepare_lav_syntax <- function(mapping, interactionsMapping) {
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
                       paste(c(grep("^x", colnames(mapping), value = TRUE),
                               apply(interactionsMapping, 1,
                                     function(x, names) {
                                       paste(names[x != 0], collapse = ":")
                                     }, names = colnames(interactionsMapping))),
                             collapse = " + "))
  return(paste(c(measurement, structural), collapse = "\n"))
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
