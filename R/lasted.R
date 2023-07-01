#' @title Auxiliary functions
#' @description
#' Computes time difference between moment of calling it and a given timestamp.
#' @param start a timestamp, typically result of calling [Sys.time]
#' @returns a string with information about the time difference
lasted <- function(start) {
  stopifnot(inherits(start, "POSIXt"), length(start) == 1L)
  lasted <- difftime(Sys.time(), start, units = "mins")
  lastedMins <- floor(lasted)
  lastedSecs <- paste0("00", round(60 * (lasted - lastedMins), 0))
  return(paste0("(lasted ", lastedMins, ":",
                 substr(lastedSecs, nchar(lastedSecs) - 1L, nchar(lastedSecs)),
                " mins.)"))
}
