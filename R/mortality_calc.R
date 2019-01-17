#' Calculate hazards for mpaf study
#'
#' @param z model matrix
#' @param cf named vector \eqn{\gamma} of parameter estimates
#'
#' @return named matrix of hazards
#' @keywords internal
mpaf_lambda <- function(z, cf) {
  as.vector(exp(z %*% -cf))
}

#' Calculate survivals for mpaf study
#'
#' @param lambda vector of hazards corresponding to \code{ID}s and
#'   \code{PERIOD}s
#' @param ID vector of IDs
#' @param PERIOD vector of periods
#' @param dt delta-times -- length of each period
#'
#' @return named matrix of survivals
#' @keywords internal
mpaf_S <- function(lambda, ID, PERIOD, dt) {
  if (any(tapply(PERIOD, ID, is.unsorted)))
    stop("Periods must be in ascending order for each ID to calculate survival")

  exp(-stats::ave(lambda, ID, FUN = function(lda) cumsum(dt * lda)))
}

#' Calculate survival differences for mpaf study
#'
#' Assumes arguments are sorted as in \code{\link{mpaf_S}}.
#'
#' @param S vector of survivals
#' @param ID vector of corresponding IDs
#' @param PERIOD vector of corresponding periods
#'
#' @return vector of survival differences over the course of each period
#' @keywords internal
mpaf_delta_S <- function(S, ID, PERIOD) {
  # note that survival at time zero is 1
  stats::ave(S, ID, FUN = function(S_i.) -diff(c(1, S_i.)))
}

#' Calculate average mortality $I$ for mpaf study
#'
#' @param dS vector of survival differences -- if cumulative mortality over time
#'   \eqn{(0, a_j]} is desired, then \code{dS} should be \code{1-S}, but if
#'   mortality over a period \eqn{(a_{j-1}, a_j]} is desired, then \code{dS}
#'   should be \eqn{\Delta{S}} as output by \code{\link{mpaf_delta_S}}.
#' @param ID vector of corresponding IDs
#' @param PERIOD vector of corresponding periods
#'
#' @section Warning: the names of the vector are exactly the levels of
#'   \code{PERIOD}, however this is a misnomer if calculating cumulative
#'   mortalities.
#'
#' @return named vector of mortalities $I$ (see Warning on names)
#' @keywords internal
mpaf_I <- function(dS, PERIOD) {
  as.vector(tapply(dS, PERIOD, mean))
}

#' Calculate gradient of hazard for mpaf study
#'
#' @param z model matrix
#' @param lambda vector of hazards
#'
#' @return matrix of gradients, where each row is a gradient for the
#'   corresponding hazard
#' @keywords internal
mpaf_grad_lambda <- function(z, lambda) {
  lambda * z # sic elementwise multiplication
}

#' Calculate gradient of survivals for mpaf study
#'
#' Assumes arguments are sorted as in \code{\link{mpaf_S}}.
#'
#' @param grad_lambda matrix of hazard gradients
#' @param S vector of survivals
#' @param ID corresponding vector of IDs
#' @param PERIOD corresponding vector of periods
#' @param dt delta-times -- length of each period
#'
#' @return matrix of gradients, where each row is a gradient for the
#'   corresponding survival
#' @keywords internal
mpaf_grad_S <- function(grad_lambda, S, ID, PERIOD, dt) {
  grad_lambda_sum <- apply(grad_lambda, 2, function(gl_..r)
    stats::ave(gl_..r, ID, FUN = function(gl_i.r) cumsum(gl_i.r * dt))
  )

  -S * grad_lambda_sum
}

#' Calculate differences of gradient of survivals for mpaf study
#'
#' Assumes arguments are sorted as in \code{\link{mpaf_S}}.
#'
#' @param grad_S matrix of survival gradients
#' @param ID corresponding vector of IDs
#' @param PERIOD corresponding vector of periods
#'
#' @return matrix of gradient differences, where each row is difference for the
#'   corresponding survival gradient
#' @keywords internal
mpaf_grad_delta_S <- function(grad_S, ID, PERIOD) {
  # note grad_survival at time zero is zero
  apply(grad_S, 2, function(grad_S_..r)
    stats::ave(grad_S_..r, ID,
               FUN = function(grad_S_i.r) -diff(c(0, grad_S_i.r)))
  )
}

#' Calculate average mortality gradient \eqn{\nabla{I}} for mpaf study
#'
#' @param grad_dS matrix of gradient survival differences -- if cumulative
#'   mortality over time $(0, a_j]$ is desired, then \code{dS} should be
#'   \code{-grad_S}, but if mortality over a period $(a_{j-1}, a_j]$ is desired,
#'   then \code{dS} should be the differences as output by
#'   \code{\link{mpaf_grad_delta_S}}.
#' @param ID vector of corresponding IDs
#' @param PERIOD vector of corresponding periods
#'
#' @section Warning: the names of the vector are exactly the levels of
#'   \code{PERIOD}, however this is a misnomer if calculating cumulative
#'   mortalities.
#'
#' @return named matrix of mortalities $I$ (see Warning on names)
#' @keywords internal
mpaf_grad_I <- function(grad_dS, ID, PERIOD) {
  apply(grad_dS, 2, function(grad_dS_..r) tapply(grad_dS_..r, PERIOD, mean))
}
