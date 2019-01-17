#' Calculate hazards for dpaf study
#'
#' @param z model matrix
#' @param cf_l named list of vectors \eqn{\gamma^D, \gamma^M} of parameter
#'   estimates
#'
#' @return named list of hazards
#' @keywords internal
dpaf_lambda <- function(z, cf_l) {
  lapply(cf_l, function(cf) as.vector(exp(z %*% -cf)))
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
#' Assumes arguments are sorted as in \code{\link{dpaf_S}}.
#'
#' @param Sp vector of survival products
#' @param ID vector of corresponding IDs
#' @param PERIOD vector of corresponding periods
#'
#' @return vector of survival product differences over the course of each period
#' @keywords internal
dpaf_delta_Sp <- function(Sp, ID, PERIOD) {
  # note that survival at time zero is 1
  stats::ave(Sp, ID, FUN = function(S_i.) -diff(c(1, S_i.)))
}

#' Calculate average/expected disease incidence \eqn{I} for dpaf study
#'
#' @param lambda_l named list of hazards
#' @param dSp vector of survival differences -- if expected disease incidence
#'   over a period \eqn{(a_{j-1}, a_j]} is desired, then \code{dSp} should be
#'   \eqn{\Delta{S}} as output by \code{\link{dpaf_delta_Sp}}. If expected
#'   disease incidence over a period \eqn{(0, a_j]} is desired, then calculate
#'   \eqn{I} as before and then take the sum.
#' @param ID vector of corresponding IDs
#' @param PERIOD vector of corresponding periods
#'
#' @section Warning: the names of the vector are exactly the levels of
#'   \code{PERIOD}, however this is a misnomer if calculating cumulative
#'   mortalities.
#'
#' @return named vector of mortalities \eqn{I} (see Warning on names)
#' @keywords internal
dpaf_I <- function(lambda_l, dSp, PERIOD) {
  as.vector(tapply(
    lambda_l[["disease"]] / do.call(`+`, lambda_l) * dSp,
    PERIOD, mean
  ))
}

#' Calculate hazard gradients for dpaf study
#'
#' @param z model matrix
#' @param cf_l named list of vectors \eqn{\gamma^D, \gamma^M} of parameter
#'   estimates
#'
#' @return named list of hazard gradients (in matrix form)
#' @keywords internal
dpaf_grad_lambda <- function(z, cf_l) {
  lapply(cf_l, `*`, z) # sic elementwise multiplication
}

#' Calculate gradient of survivals for mpaf study
#'
#' Assumes arguments are sorted as in \code{\link{dpaf_S}}.
#'
#' @param grad_lambda_l named list of matrices of hazard gradients
#' @param S_l named list of vector of survivals
#' @param ID corresponding vector of IDs
#' @param PERIOD corresponding vector of periods
#' @param dt delta-times -- length of each period
#'
#' @return names list of matrix of gradients, where each row in each list item
#'   is the gradient of the corresponding survival
#' @keywords internal
dpaf_grad_S <- function(grad_lambda_l, S_l, ID, PERIOD, dt) {
  grad_lambda_sum <- lapply(grad_lambda_l, function(grad_lambda)
    apply(grad_lambda, 2, function(gl_..r)
      stats::ave(gl_..r, ID, FUN = function(gl_i.r) cumsum(gl_i.r * dt))
  ))

  lapply(mapply(`*`, S_l, grad_lambda_sum, SIMPLIFY = FALSE), `-`)
}

#' Calculate gradient of survival products
#'
#' Assumes arguments are sorted as in \code{\link{dpaf_S}}.
#'
#' @param grad_S_l named list of survival gradient matrices
#' @param Sp vector of survival products
#' @param ID corresponding vector of IDs
#' @param PERIOD corresponding vector of PERIODs
#'
#' @return named list of survival product gradient matrices
#' @keywords internal
dpaf_grad_delta_Sp <- function(grad_S_l, Sp, ID, PERIOD) {
  # TODO check gradient vector calculus
  # sic elementwise multiplication
  grad_Sp_l <- lapply(grad_S_l, `*`, Sp)

  # note gradient survival at time zero is zero
  lapply(grad_Sp_l, function(grad_Sp)
    apply(grad_Sp, 2, function(grad_Sp_..r)
      stats::ave(grad_Sp_..r, ID,
                 FUN = function(grad_Sp_i.r) -diff(c(0, grad_Sp_i.r))))
  )
}

#' Calculate average/expected disease incidence \eqn{I} for dpaf study
#'
#' @param grad_lambda_l named list of hazard gradient matrices
#' @param grad_dSp named list of survival product difference matrices
#' @param lambda_l named list of hazards
#' @param dSp vector of survival differences
#' @param ID vector of corresponding IDs
#' @param PERIOD vector of corresponding periods
#'
#' @return named vector of mortalities \eqn{I} (see Warning on names)
#' @keywords internal
dpaf_grad_I <- function(grad_lambda_l, grad_dSp_l, lambda_l, dSp, PERIOD) {
  summands_d <- grad_lambda_l[["disease"]] *
    lambda_l[["mortality"]] / do.call(`+`, lambda_l)^2 * dSp +
    grad_dSp_l[["disease"]] * lambda_l[["disease"]] / do.call(`+`, lambda_l)

  summands_m <- -grad_lambda_l[["mortality"]] *
    lambda_l[["disease"]] / do.call(`+`, lambda_l)^2 * dSp +
    grad_dSp_l[["mortality"]] * lambda_l[["disease"]] / do.call(`+`, lambda_l)

  lapply(
    list("disease" = summands_d, "mortality" = summands_m),
    function(smnd) apply(smnd, 2, function(smnd_..r)
      as.vector(tapply(smnd_..r, PERIOD, mean)))
  )
}

#' Calculate gradients of \eqn{\log(1-PAF)} wrt disease and mortality
#' coefficients
#'
#' @param gi_mod named list of gradients of \eqn{I^*}
#' @param i_mod \eqn{I^*}
#' @param gi_raw named list of gradients of \eqn{I}
#' @param i_raw \eqn{I}
#'
#' @return a named list of gradient vectors of the iPAF
#' @keywords internal
dpaf_gipaf <- function(gi_mod, i_mod, gi_raw, i_raw) {
  mapply(
    `-`,
    lapply(gi_mod, `/`, i_mod),
    lapply(gi_raw, `/`, i_raw),
    SIMPLIFY = FALSE
  )
}

#' Calculate gradients of PAF wrt disease and mortality coefficients
#'
#' @param gi_mod named list of gradients of \eqn{I^*}
#' @param i_mod \eqn{I^*}
#' @param gi_raw named list of gradients of \eqn{I}
#' @param i_raw \eqn{I}
#'
#' @return a named list of gradient vectors of the PAF
#' @keywords internal
dpaf_gpaf <- function(gi_mod, i_mod, gi_raw, i_raw) {
  mapply(
    function(gi_mod_x, gi_raw_x)
      (gi_raw_x * i_mod - gi_mod_x * i_raw) / i_raw**2,
    gi_mod, gi_raw,
    SIMPLIFY = FALSE
  )
}
