#' @title Data generation
#' @description
#' Generates a matrix with the sampled values of latent variables
#' @param nObs number of observations to be generated
#' @param nLV number of exogenous latent variables
#' @param nLI number of latent interactions
#' @param r correlation between (each pair of) exogenous latent variables
#' @param partR partial correlation between (each) exogenous latent variable and
#' dependent latent variable
#' @param partRI partial correlation between (each) latent interaction variable
#' and dependent latent variable
#' @returns a numeric matrix with `nObs` rows and `nLV + nLI + 1` columns
#' (values of `r`, `partR`, `partRI`, as well as mapping of latent variables
#' on interaction variables, slope parameters and error term standard deviation
#' used to generate the dependent latent variable, and covariance matrix of
#' latent variables are returned as attributes)
#' @seealso [compute_error_std_dev], [compute_interaction_slope],
#' [compute_dep_std_var]
#' @examples
#' x <- generate_latent(nObs = 10^5, nLV = 3, nLI = 2,
#'                      r = 0.3, partR = 0.3, partRI = 0.2)
#' # a matrix to store partial correlations
#' res = diag(ncol(x))
#' colnames(res) = rownames(res) = colnames(x)
#' for (i in seq_len(ncol(res))) {
#'   for (j in seq_len(i - 1L)) {
#'     res[i, j] = res[j, i] = cor(resid(lm(x[, i] ~ 0 + x[, -c(i, j)])),
#'                                 resid(lm(x[, j] ~ 0 + x[, -c(i, j)])))
#'   }
#' }
#' # correlations
#' round(cor(x), 2)
#' # partial correlations
#' round(res, 3)
#' # comparison of theoretical and empirical covariances
#' round(attributes(x)$cov, 2)
#' round(cov(x), 2)
#' @export
generate_latent <- function(nObs, nLV, nLI, r, partR, partRI = partR) {
  stopifnot(is.numeric(nObs), length(nObs) == 1,
            as.integer(nObs) == nObs, nObs > 0,
            is.numeric(nLV), length(nLV) == 1,
            as.integer(nLV) == nLV, nLV > 0,
            is.numeric(nLI), length(nLI) == 1,
            as.integer(nLI) == nLI, nLI >= 0, nLI <= nLV*(nLV - 1L)/2L,
            is.numeric(r), length(r) == 1, !is.na(r),
            r > -1, r < 1,
            is.numeric(partR), length(partR) == 1, !is.na(partR),
            partR > -1, partR < 1,
            is.numeric(partRI), length(partRI) == 1, !is.na(partRI),
            partRI > -1, partRI < 1)

  vc <- matrix(r, nrow = nLV, ncol = nLV)
  diag(vc) <- 1
  rownames(vc) <- colnames(vc) <- paste0("x", 1L:nLV)
  interactions <- utils::combn(seq_len(nLV), 2)[, seq_len(nLI), drop = FALSE]
  mapping <- matrix(0L, nrow = ncol(interactions), ncol = nLV,
                    dimnames = list(paste0(rep("xi", nLI),
                                           apply(interactions, 2, paste,
                                                 collapse = "")),
                                    rownames(vc)))
  x <- cbind(mnormt::rmnorm(nObs, mean = rep(0, nLV), varcov = vc),
             matrix(rep(NA_real_, nObs*nLI), nrow = nObs, ncol = nLI,
                    dimnames = list(NULL, rownames(mapping))))
  for (i in seq_len(nLI)) {
    x[, nLV + i] <- x[, interactions[1L, i]]*x[, interactions[2L, i]]
    mapping[i, interactions[, i]] <- 1L
  }
  errorStdDev <- compute_error_std_dev(nLV, r, partR, 1)
  slopeI <- compute_interaction_slope(nLI, r, partRI, errorStdDev)
  slopes <- c(rep(1, nLV), rep(slopeI, nLI))
  names(slopes) <- colnames(x)
  yStdDev <- compute_dep_std_var(nLV, nLI, r, 1, slopeI, errorStdDev)
  expectedValues <- c(rep(0, nLV), rep(r, nLI))
  y <- (as.vector(x %*% slopes) + stats::rnorm(nObs, 0, errorStdDev) -
          sum(expectedValues * slopes)) / yStdDev

  cov <- matrix(0, nrow = nLV + nLI + 1, ncol = nLV + nLI + 1,
                dimnames = list(c(colnames(x), "y"), c(colnames(x), "y")))
  cov[seq_len(nLV), seq_len(nLV)] <- r
  cov[nLV + seq_len(nLI), nLV + seq_len(nLI)] <- r + r^2
  diag(cov) <- c(rep(1, nLV), rep(1 + r^2, nLI), 1)
  covY <- matrix(rep(slopes / yStdDev, length(slopes)),
                  nrow = length(slopes)) * cov[-nrow(cov), -ncol(cov)]
  cov[nrow(cov), -ncol(cov)] <- cov[-nrow(cov), ncol(cov)] <- colSums(covY)

  return(structure(cbind(x, y = y),
                   r = r,
                   rPart = partR,
                   rPartI = partRI,
                   mapping = mapping,
                   slopes = slopes / yStdDev,
                   errorStdDev = errorStdDev / yStdDev,
                   cov = cov))
}
#' @title Data generation
#' @description
#' Computes the (non-normalized, compare [compute_dep_std_var]) standard
#' deviation of the latent dependent variable's error term (used internally by
#' the [generate_latent])
#' @inheritParams generate_latent
#' @inheritParams compute_dep_std_var
#' @returns a scalar
compute_error_std_dev <- function(nLV, r, partR, slope = 1) {
  stopifnot(is.numeric(nLV), length(nLV) == 1,
            as.integer(nLV) == nLV, nLV > 0,
            is.numeric(r), length(r) == 1, !is.na(r),
            r > -1, r < 1,
            is.numeric(partR), length(partR) == 1, !is.na(partR),
            partR > -1, partR < 1,
            is.numeric(slope), length(slope) == 1, !is.na(slope))
  p11 <- 1 + (nLV - 1)*r
  s1 <- p11 / (1 + (nLV - 2)*r)
  return(slope * sqrt(p11*(nLV - (nLV - 1)*s1)*(1 - partR^2)) / partR)
}
#' @title Data generation
#' @description
#' Computes the value of the slope coefficient for the interaction latent
#' variables that should be used while generating the dependent latent variable
#' (used internally by the [generate_latent])
#' @inheritParams generate_latent
#' @inheritParams compute_dep_std_var
#' @returns a scalar
compute_interaction_slope <- function(nLI, r, partRI, errorStdDev) {
  stopifnot(is.numeric(nLI), length(nLI) == 1,
            as.integer(nLI) == nLI, nLI >= 0,
            is.numeric(r), length(r) == 1, !is.na(r),
            r > -1, r < 1,
            is.numeric(partRI), length(partRI) == 1, !is.na(partRI),
            partRI > -1, partRI < 1,
            is.numeric(errorStdDev), length(errorStdDev) == 1,
            !is.na(errorStdDev), errorStdDev > 0)
  if (nLI == 0) return(NA_real_)
  pi1 <- 1 + (nLI - 1)*r + nLI*r^2
  si <- pi1 / (1 + (nLI - 2)*r + (nLI - 1)*r^2)
  return(partRI * errorStdDev / sqrt((pi1*(nLI - (nLI - 1)*si) * (1 - partRI^2))))
}
#' @title Data generation
#' @description
#' Computes standard deviation of the latent dependent variable so it can be
#' standardized to have variance of 1 before returning it by the [generate_latent]
#' @inheritParams generate_latent
#' @param slope a scalar - value of the slope coefficient for the main effects
#' (typically set to 1)
#' @param slopeI a scalar - value of the slope coefficient for the interaction
#' effects, typically returned by the [compute_interaction_slope]
#' @param errorStdDev a scalar - standard deviation of the dependent latent
#' variable error term, typically returned by the [compute_error_std_dev]
#' @returns a scalar
compute_dep_std_var <- function(nLV, nLI, r, slope, slopeI,
                                errorStdDev) {
  stopifnot(is.numeric(nLV), length(nLV) == 1,
            as.integer(nLV) == nLV, nLV > 0,
            is.numeric(nLI), length(nLI) == 1,
            as.integer(nLI) == nLI, nLI >= 0,
            is.numeric(r), length(r) == 1, !is.na(r),
            r > -1, r < 1,
            is.numeric(slope), length(slope) == 1, !is.na(slope),
            is.numeric(slopeI), length(slopeI) == 1,
            is.numeric(errorStdDev), length(errorStdDev) == 1,
            !is.na(errorStdDev), errorStdDev > 0)
  v1 <- slope^2 * nLV * (1 + (nLV - 1)*r)
  if (nLI == 0) vi <- 0
  else vi <- slopeI^2 * nLI * (1 + (nLI - 1)*r + nLI*r^2)
  return(sqrt(v1 + vi + errorStdDev^2))
}
