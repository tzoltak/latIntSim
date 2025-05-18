#' @title Preparing simulation conditions
#' @description
#' Returns names of parameters required to specify a simulation condition.
#' @returns a character vector
#' @seealso [check_conditions_names], [check_conditions]
#' @export
get_required_conditions_names <- function() {
  return(c("nObs", "nExogLVs", "corExogLVs", "nLIs", "nIndic", "lambda",
           "partCorMain", "partCorInt"))
}
#' @title Preparing simulation conditions
#' @description
#' Checks whether a given vector of names contains names of all the parameters
#' required to specify a simulation condition
#' @param parNames a character vector
#' @returns its argument (or throws an error)
#' @seealso [check_conditions], [prepare_conditions]
check_conditions_names <- function(parNames) {
  stopifnot(is.character(parNames),
            !anyNA(parNames))
  reqNames <- get_required_conditions_names()
  if (!all(reqNames %in% parNames)) {
    stop(paste0("Some required parameters has not been specified: '",
                paste(setdiff(reqNames, parNames), collapse = "', '"), "'."))
  }
  if (any(!(parNames %in% reqNames))) {
    message(paste0("Some additional parameters have been specified: '",
                   paste(setdiff(parNames, reqNames), collapse = "', '"), "'. ",
                   "These ones won't affect the simulation design."))
  }
  if (any(duplicated(parNames))) {
    stop(paste0("Each argument must be specified as a column in exactly one data frame (duplicated names: '",
                paste(parNames[duplicated(parNames)], collapse = "', '"), "')."))
  }
  return(parNames)
}
#' @title Preparing simulation conditions
#' @description
#' Checks whether columns in the data frame storing specifications of simulation
#' conditions have a correct set o values
#' @param conditions a data frame
#' @returns its argument (or throws an error)
#' @seealso [check_conditions_names], [prepare_conditions]
#' @export
check_conditions <- function(conditions) {
  stopifnot(is.data.frame(conditions))
  check_conditions_names(names(conditions))

  stopifnot(!anyNA(conditions[, get_required_conditions_names()]),
            is.numeric(conditions$nObs), all(conditions$nObs > 0),
            all(as.integer(conditions$nObs) == conditions$nObs),
            is.numeric(conditions$nExogLVs), all(conditions$nExogLVs > 0),
            all(conditions$nExogLVs < 10), # get_model_pars() won't work with more than 9
            all(as.integer(conditions$nExogLVs) == conditions$nExogLVs),
            is.numeric(conditions$nLIs),
            all(conditions$nLIs <= (conditions$nExogLVs * (conditions$nExogLVs - 1) / 2)),
            all(as.integer(conditions$nLIs) == conditions$nLIs),
            is.numeric(conditions$nIndic), all(conditions$nIndic > 1),
            all(as.integer(conditions$nIndic) == conditions$nIndic),
            is.numeric(conditions$corExogLVs),
            all(conditions$corExogLVs > -1), all(conditions$corExogLVs < 1),
            is.numeric(conditions$lambda),
            all(conditions$lambda > -1), all(conditions$lambda < 1),
            is.numeric(conditions$partCorMain),
            all(conditions$partCorMain > -1), all(conditions$partCorMain < 1),
            is.numeric(conditions$partCorInt),
            all(conditions$partCorInt > -1), all(conditions$partCorInt < 1))
  return(conditions)
}
