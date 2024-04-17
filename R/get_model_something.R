#' @title Extracting results from the estimated models
#' @description
#' Gets summaries from estimated models.
#' @param x an object of one of the following classes: *lavaan*, *mplusObject*,
#' *miive* or *seminr_model*
#' @returns a data frame with the following columns:
#' \describe{
#'   \item{time}{estimation time in seconds (only for *lavaan* and *mplusObject*)}
#'   \item{converged}{logical value indicating whether model converged}
#'   \describe{
#'     \item{*lavaan*}{taken from information included explicitly in the object}
#'     \item{*mplusObject*}{determined by checking whether standard errors (or
#'                          posterior standard deviations) of model parameters
#'                          are reported in the results - if yes, it is assumed
#'                          that model converged}
#'     \item{*miive*}{determined by checking whether parameters covariance
#'                    matrix is reported for every equation estimated in the
#'                    model - if yes, it is assumed that model converged}
#'     \item{*seminr_model*}{determined by checking whether reported number of
#'                           iterations of model estimation is no greater than
#'                           the maximal number of estimations (both are
#'                           reported explicitly in the object)}
#'   }
#'   \item{niter}{the number of iterations spent to estimated the model (not
#'                reported for *miive*, because its estimation procedure is not
#'                iterative)}
#'   \item{warnings}{warning messages, if any (these are not captured through
#'                   the estimation, but retrieved from the object storing the
#'                   estimated model - for *miive* the only possible warning is
#'                   that about negative variance estimates, for *seminr* no
#'                   warnings are reported)}
#' }
get_model_summary <- function(x) {
  if (inherits(x, "try-error")) {
    return(data.frame(time = NA_real_,
                      converged = FALSE,
                      niter = NA_integer_,
                      warnings = as.character(x),
                      stringsAsFactors = FALSE))
  } else if (inherits(x, "lavaan")) {
    variances <- lavaan::summary(x)$pe
    variances <- variances[variances$op == "~~" & variances$lhs == variances$rhs, ]
    return(data.frame(time = x@timing$total,
                      converged = x@optim$converged,
                      niter = x@optim$iterations,
                      warnings =
                        paste0(x@optim$warn.txt,
                               ifelse(x@optim$warn.txt == "" |
                                        all(variances$est >= 0), "", "; "),
                               ifelse(all(variances$est >= 0), "",
                                      "Some estimated variances are negative.")),
                      stringsAsFactors = FALSE))
  } else if (inherits(x, "mplusObject")) {
    return(data.frame(time = get_mplus_time(x$results$output),
                      converged = any(c("se", "posterior_sd") %in%
                                        names(x$results$parameters$unstandardized)),
                      niter = get_mplus_niter(x$results$output),
                      warnings = sub("USE THE FBITERATIONS OPTION TO INCREASE THE NUMBER OF ITERATIONS BY A FACTOR\nOF AT LEAST TWO TO CHECK CONVERGENCE AND THAT THE PSR VALUE DOES NOT INCREASE.",
                                     "",
                                     paste(c(unlist(x$results$errors),
                                             unlist(x$results$warnings)),
                                           collapse = "\n")),
                      stringsAsFactors = FALSE))
  } else if (inherits(x, "miive")) {
    return(data.frame(time = NA_real_,
                      converged = all(sapply(x$eqn, function(x) !anyNA(x$XX1))),
                      niter = NA_integer_,
                      warnings = ifelse(!all(x$v$coefficients >= 0),
                                        "Some estimated variances are negative.",
                                        ""),
                      stringsAsFactors = FALSE))
  } else if (inherits(x, "seminr_model")) {
    return(data.frame(time = NA_real_,
                      converged = x$iterations <= x$settings$maxIt,
                      niter = x$iterations,
                      warnings = NA_character_,
                      stringsAsFactors = FALSE))
  } else if (inherits(x, "cSEMResults")) {
    return(data.frame(time = NA_real_,
                      converged = x$Information$Weight_info$Convergence_status,
                      niter = x$Information$Weight_info$Number_iterations,
                      warnings = ifelse(x$Information$Information_resample$Number_of_admissibles <
                                          x$Information$Information_resample$Number_of_runs,
                                        paste0(x$Information$Information_resample$Number_of_runs -
                                                 x$Information$Information_resample$Number_of_admissibles,
                                               " out of ",
                                               x$Information$Information_resample$Number_of_runs,
                                               " bootstrap samples failed to converge."), ""),
                      stringsAsFactors = FALSE))
  } else {
    stop("Unsupported type of model.")
  }
}
#' @title Extracting results from the estimated models
#' @description
#' Gets parameters from estimated models.
#' @param x an object of one of the following classes: *lavaan*, *mplusObject*,
#' *miive* or *seminr_model*
#' @returns a list of two data frames:
#'
#' \strong{structural}
#' \describe{
#'   \item{DV}{dependent latent variable}
#'   \item{IV}{independent latent variable (predictor)}
#'   \item{est}{regression (slope) coefficient}
#'   \item{se}{standard error of `est`}
#'   \item{std.all}{standardized regression (slope) coefficient}
#' }
#' \strong{measurement}
#' \describe{
#'   \item{LV}{latent variable}
#'   \item{obsIndic}{observed indicator}
#'   \item{est}{factor loading}
#'   \item{se}{standard error of `est`}
#'   \item{lambda}{standardized factor loading}
#' }
#' @details
#' \itemize{
#'   \item{For *seminr_model* objects (PLS models) only standardized structural
#'         coefficients are reported in model results so in `structural` values
#'         in `est` always equal those in `std.all` for these models.}
#'   \item{For *miive* objects no loadings are reported for first observed
#'         indicators of each latent variable as these serve as so-called
#'         *scaling indicators* and play a special role in model parameters
#'         estimation.}
#' }
get_model_pars <- function(x) {
  # stop("Należałoby tu przenieść obliczanie BS parametrów wystandaryzowanych (obecnie robię to już przy analizie wyników), bo w przypadku Mplusa, z jego ograniczoną precyzją zapisu oszacowań, robienie tego post factum w oparciu o iloraz parametru standaryzowanego i niestandaryzowanego jest zawodne, gdy oszacowanie jest bardzo bliskie 0.")
  if (inherits(x, "try-error")) {
    structural <- data.frame(DV = vector(mode = "character", length = 0L),
                             IV = vector(mode = "character", length = 0L),
                             est =  vector(mode = "numeric", length = 0L),
                             se =  vector(mode = "numeric", length = 0L),
                             std.all =  vector(mode = "numeric", length = 0L),
                             std.all.se = vector(mode = "numeric", length = 0L))
    measurement <- data.frame(LV = vector(mode = "character", length = 0L),
                              obsIndic = vector(mode = "character", length = 0L),
                              est = vector(mode = "numeric", length = 0L),
                              se = vector(mode = "numeric", length = 0L),
                              lambda = vector(mode = "numeric", length = 0L),
                              lanbda.se = vector(mode = "numeric", length = 0L))
  } else if (inherits(x, "lavaan")) {
    x <- lavaan::summary(x, standardized = TRUE)$pe[, c("lhs", "op", "rhs",
                                                        "est", "se", "std.all")]
    structural <- x[grepl("^y[[:digit:]]*$", x$lhs) & x$op == "~",
                    colnames(x) != "op"]
    names(structural) <- c("DV", "IV", "est", "se", "std.all")
    structural$std.all.se <- structural$se * structural$std.all / structural$est
    structural$IV[grep(":", structural$IV)] <-
      paste0("xi", gsub("[^[:digit:]]", "", structural$IV[grep(":", structural$IV)]))
    measurement <- x[x$op == "=~" & !grepl("^xi", x$lhs), colnames(x) != "op"]
    names(measurement) <- c("LV", "obsIndic", "est", "se", "lambda")
    measurement$lambda.se <- measurement$se * measurement$lambda / measurement$est
  } else if (inherits(x, "miive")) {
    stdDevs <- data.frame(DV = sub("~~.*$", "", names(x$v$coefficients)),
                          IV = sub("^.*~~", "", names(x$v$coefficients)),
                          sd = suppressWarnings(sqrt(x$v$coefficients)),
                          stringsAsFactors = FALSE)
    stdDevs <- stdDevs[stdDevs$DV == stdDevs$IV, ]
    x <- data.frame(DV = sub("~.*$", "", names(stats::coef(x))),
                    IV = sub("^.*~", "", names(stats::coef(x))),
                    est = stats::coef(x),
                    se = sqrt(diag(x$coefCov)),
                    std.all = NA_real_,
                    stringsAsFactors = FALSE)
    structural <- x[grepl("^y[[:digit:]]*$", x$DV) & x$IV != "1", ]
    structural <- merge(structural, stdDevs[, c("DV", "sd")], by = "DV",
                        all.x = TRUE)
    structural <- merge(structural, stdDevs[, c("IV", "sd")], by = "IV",
                        all.x = TRUE)
    structural$std.all <- structural$est * structural$sd.y / structural$sd.x
    structural <- structural[, c("IV", "DV", "est", "se", "std.all")]
    structural$std.all.se <- structural$se * structural$std.all / structural$est

    measurement <- x[grepl("^[xy][[:digit:]]*_[[:digit:]]+$", x$DV) & x$IV != "1", ]
    measurement <- merge(measurement, stdDevs[, c("DV", "sd")], by = "DV",
                         all.x = TRUE)
    measurement <- merge(measurement, stdDevs[, c("IV", "sd")], by = "IV",
                         all.x = TRUE)
    measurement$std.all <- measurement$est * measurement$sd.y / measurement$sd.x
    measurement <- measurement[, c("IV", "DV", "est", "se", "std.all")]
    names(measurement) <- c("LV", "obsIndic", "est", "se", "lambda")
    measurement$lambda.se <- measurement$se * measurement$lambda / measurement$est
  } else if (inherits(x, "mplusObject")) {
    if (any(c("se", "posterior_sd") %in%
            names(x$results$parameters$unstandardized))) {
      if ("posterior_sd" %in% names(x$results$parameters$stdyx.standardized)) {
        names(x$results$parameters$stdyx.standardized) <-
          sub("^posterior_sd$", "se", names(x$results$parameters$stdyx.standardized))
        names(x$results$parameters$unstandardized) <-
          sub("^posterior_sd$", "se", names(x$results$parameters$unstandardized))
      }
      stdX <- x$results$parameters$stdyx.standardized[, c("paramHeader", "param", "est", "se")]
      names(stdX) <- c("paramHeader", "param", "std.all", "std.all.se")
      x <- merge(x$results$parameters$unstandardized,
                 stdX, by = c("paramHeader", "param"), all.x = TRUE)
      structural <- x[grep("\\.ON$", x$paramHeader),
                      c("paramHeader", "param", "est", "se", "std.all", "std.all.se")]
      structural$paramHeader <- tolower(sub("\\.ON", "", structural$paramHeader))
      structural$param <- tolower(structural$param)
      names(structural) <- c("DV", "IV", "est", "se", "std.all", "std.all.se")
      measurement <- x[grep("\\.BY$", x$paramHeader),
                       c("paramHeader", "param", "est", "se", "std.all", "std.all.se")]
      measurement$paramHeader <- tolower(sub("\\.BY", "", measurement$paramHeader))
      measurement$param <- tolower(measurement$param)
      names(measurement) <- c("LV", "obsIndic", "est", "se", "lambda", "lambda.se")
    } else {
      structural <- data.frame(DV = vector(mode = "character", length = 0L),
                               IV = vector(mode = "character", length = 0L),
                               est = vector(mode = "numeric", length = 0L),
                               se = vector(mode = "numeric", length = 0L),
                               std.all = vector(mode = "numeric", length = 0L),
                               std.all.se = vector(mode = "numeric", length = 0L))
      measurement <- data.frame(LV = vector(mode = "character", length = 0L),
                                obsIndic = vector(mode = "character", length = 0L),
                                est = vector(mode = "numeric", length = 0L),
                                se = vector(mode = "numeric", length = 0L),
                                lambda = vector(mode = "numeric", length = 0L),
                                lanbda.se = vector(mode = "numeric", length = 0L))
    }
  } else if (inherits(x, "boot_seminr_model")) {
    stopifnot("Models with cross-loadings are not supported for PLS models." =
                all(rowSums(x$outer_weights != 0) <= 1))
    structural <- data.frame(DV = rep(sub("^(y[[:digit:]]*) PLS Est\\.$", "\\1",
                                          grep("^y[[:digit:]]* PLS Est\\.$",
                                               colnames(x$paths_descriptives), value = TRUE)),
                                      each = nrow(x$paths_descriptives)),
                             IV = rep(rownames(x$paths_descriptives),
                                      sum(grepl("^y[[:digit:]]* PLS Est\\.$",
                                               colnames(x$paths_descriptives)))),
                             est = as.vector(x$paths_descriptives[, grep("^y[[:digit:]]* PLS Est\\.$",
                                                                         colnames(x$paths_descriptives))]),
                             se = as.vector(x$paths_descriptives[, grep("^y[[:digit:]]* Boot SD$",
                                                                          colnames(x$paths_descriptives))]))
    structural$std.all <- structural$est
    structural$std.all.se <- structural$se
    structural$IV <- sub("^x([[:digit:]]+)\\*x([[:digit:]]+)$",
                        "xi\\1\\2", structural$IV)
    measurement <- data.frame(LV = apply(x$outer_weights, 1,
                                         function(x, colnames) colnames[x != 0],
                                         colnames = colnames(x$outer_weights)),
                              obsIndic = rownames(x$outer_weights),
                              est = rowSums(x$outer_weights),
                              se = rowSums(x$weights_descriptives[, grep(" SD$",
                                                                         colnames(x$weights_descriptives))]),
                              lambda = rowSums(x$outer_loadings))
    measurement$lambda.se <- measurement$se * measurement$lambda / measurement$est
    measurement <- measurement[grep("^[xy][[:digit:]]*$", measurement$LV), ]
  } else if (inherits(x, "cSEMResults")) {
    x <- cSEM::infer(x)
    structural <- data.frame(DV = sub(" ~.*$", "", names(x$Path_estimates$mean)),
                             IV = sub("^.* ~ ", "", names(x$Path_estimates$mean)),
                             est = x$Path_estimates$mean - x$Path_estimates$bias,
                             se = x$Path_estimates$sd)
    structural$IV[grep("\\.", structural$IV)] <-
      paste0("xi", gsub("[^[:digit:]]", "", structural$IV[grep("\\.", structural$IV)]))
    structural$std.all <- structural$est
    structural$std.all.se <- structural$se
    measurement <- data.frame(LV = sub(" =~.*$", "", names(x$Loading_estimates$mean)),
                              obsIndic = sub("^.*=~ ", "", names(x$Loading_estimates$mean)),
                              est = x$Weight_estimates$mean - x$Weight_estimates$bias,
                              se = x$Weight_estimates$sd,
                              lambda = x$Loading_estimates$mean - x$Loading_estimates$bias)
    measurement$lambda.se <- measurement$se * measurement$lambda / measurement$est
  } else {
    stop("Unsupported type of model.")
  }
  return(list(structural = structural,
              measurement = measurement))
}
#' @title Auxiliary functions
#' @description
#' Extracts estimation time from mplus output.
#' @param x `mplusObject$results$output`
#' @returns a scalar (time in seconds)
get_mplus_time <- function(x) {
  x <- grep("Elapsed Time", x, value = TRUE)
  if (length(x) != 1L) return(NA_real_)
  x <- strsplit(sub("^.*Elapsed Time: *([[:digit:]:]+)$", "\\1", x), ":")[[1]]
  return(sum(as.numeric(x) * c(60^2, 60, 1)))
}
#' @title Auxiliary functions
#' @description
#' Extracts number of iterations from mplus output.
#' @param x `mplusObject$results$output`
#' @returns a scalar (time in seconds)
#' @details
#' Works only for models for which TECH8 reporting was turned on.
get_mplus_niter <- function(x) {
  x <- x[grep("^TECHNICAL 8 OUTPUT$", x):(grep("Beginning Time", x) - 1L)]
  x <- x[grep("^ +[[:digit:]]+ ", x)]
  return(as.integer(sub("^ +([[:digit:]]+).*$", "\\1", x[length(x)])))
}
