#' @title Preparing simulation conditions
#' @description
#' Prepares a data frame storing simulation conditions. Its arguments should be
#' data frames giving combinations of simulation parameters.
#' @param constant a one-row data frame containing values of parameters that
#' should be hold constant across all the simulation conditions
#' @param ... data frames with combinations of values of the other parameters
#' @details
#' Each parameter in the list below should be given as a column in exactly one
#' of the data frames provided as function arguments:
#'
#' \describe{
#'   \item{nObs}{ (an integer) - number of observations}
#'   \item{nExogLVs}{ (an integer) - number of exogenous latent variables)}
#'   \item{corExogLVs}{ (a number) - value of pairwise correlations between exogenous latent variables)}
#'   \item{nLIs}{ (an integer) - number of latent (first-order) interactions}
#'   \item{nIndic}{ (an integer) - number of observed indicators for each exogenous latent variable and for the dependent latent variable}
#'   \item{lambda}{ (a number) - value of standardized factor loadings (the same for all observed indicators of all latent variables)}
#'   \item{partCorMain}{ (a number) - value of partial correlation between (each) exogenous latent variable and the dependent latent variable, given all the other exogenous latent variables (and latent interactions)}
#'   \item{partCorInt}{ (a number) - value of partial correlation between (each) latent interaction variable and the dependent latent variable, given all the other latent interaction variables (and exogenous latent variables)}
#' }
#' @returns a data frame
#' @seealso [check_conditions], [check_conditions_names]
#' @examples
#' prepare_conditions(constant = data.frame(nExogLVs = 3),
#'                    data.frame(nObs = c(500, 1000)),
#'                    data.frame(corExogLVs = c(0.3, 0.5)),
#'                    data.frame(nLIs = c(1, 2, 3)),
#'                    data.frame(nIndic = c(3, 5)),
#'                    data.frame(lambda = c(0.4, 0.8)),
#'                    data.frame(partCorMain = c(0.1, 0.3, 0.5),
#'                               partCorInt = c(0.1, 0.3, 0.5)))
#' @export
prepare_conditions <- function(constant, ...) {
  stopifnot(is.data.frame(constant) | is.null(constant))
  if (!is.null(constant)) {
    stopifnot(nrow(constant) == 1L)
  }
  varying <- list(...)
  stopifnot("All arguments must be data frames" =
              all(sapply(varying, is.data.frame)))
  check_conditions_names(c(names(constant), unlist(lapply(varying, names))))

  conditions <- do.call(tidyr::expand_grid,
                        append(list(constant), varying))
  check_conditions(conditions)
  return(conditions)
}
