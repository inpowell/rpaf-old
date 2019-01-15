#' Calculate hazards for dpaf study
#'
#' @param z model matrix
#' @param cf_l named list of vectors \eqn{\gamma^D, \gamma^M} of parameter
#'   estimates
#'
#' @return named list of hazards
#' @keywords internal
dpaf_lambda <- function(z, cf_l) {
  lapply(cf_l, function(cf) as.vector(exp(z %*% -gamma)))
}

#' Calculate survivals for mpaf study
#'
#' @param lambda_l named list of vectors of hazards corresponding to \code{ID}s
#'   and \code{PERIOD}s
#' @param ID vector of IDs
#' @param PERIOD vector of periods
#' @param dt delta-times -- length of each period
#'
#' @return named list of survivals
#' @keywords internal
dpaf_S <- function(lambda_l, ID, PERIOD, dt) {
  if (any(tapply(PERIOD, ID, is.unsorted)))
    stop("Periods must be in ascending order for each ID to calculate survival")

  lapply(lambda_l, function(lambda)
    exp(-stats::ave(lambda, ID, FUN = function(lda) cumsum(dt * lda))))
}

#' Calculate disease-free survival (**S**urvival **p**roduct)
#'
#' @param S_l names list of survivals \eqn{S^D, S^M}
#'
#' @return vector of survival products \eqn{S = S^D S^M}
#' @keywords internal
dpaf_Sp <- function(S_l) {
  do.call(`*`, S_l)
}

#' Calculate survival differences for dpaf study
#'
#' @param Sp vector of survival products
#' @param ID vector of corresponding IDs
#' @param PERIOD vector of corresponding periods
#'
#' @return vector of survival product differences over the course of each period
#' @keywords internal
dpaf_delta_Sp <- function(Sp, ID, PERIOD) {
  if (any(tapply(PERIOD, ID, is.unsorted)))
    stop("Periods must be in ascending order for each ID to calculate survival")
  # note that survival at time zero is 1
  stats::ave(Sp, ID, FUN = function(S_i.) -diff(c(1, S_i.)))
}

#' Calculate average/expected disease-free survival \eqn{I} for dpaf study
#'
#' @param lambda_l named list of hazards
#' @param dSp vector of survival differences -- if \eqn{I} over time
#'   \eqn{(0, a_j]} is desired, then \code{dS} should be \code{1-S}, but if
#'   mortality over a period \eqn{(a_{j-1}, a_j]} is desired, then \code{dSp}
#'   should be \eqn{\Delta{S}} as output by \code{\link{dpaf_delta_Sp}}.
#' @param ID vector of corresponding IDs
#' @param PERIOD vector of corresponding periods
#'
#' @section Warning: the names of the vector are exactly the levels of
#'   \code{PERIOD}, however this is a misnomer if calculating cumulative
#'   mortalities.
#'
#' @return named vector of mortalities \eqn{I} (see Warning on names)
#' @keywords internal
dpaf_I <- function(lambda_l, dSp, ID, PERIOD) {
  if (any(tapply(PERIOD, ID, is.unsorted)))
    stop("Periods must be in ascending order for each ID to calculate I")

  as.vector(tapply(
    lambda_l[["disease"]] / do.call(`+`, lambda_l) * dSp,
    PERIOD, mean
  ))
}
