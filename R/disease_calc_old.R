###### Low-level, intermediate calculation functions

### Key: (R symbol -- meaning -- paper algebraic describer)
# hz -- "hazard"     -- lambda
# sv -- "survival"   -- S
# i  -- "mortality"  -- I
# g* -- "gradient"   -- \frac{\partial (*)}{\partial \gamma}


#' Calculate hazards for dpaf study
#'
#' @param z model matrix
#' @param cf named matrix of coefficients, as in \code{\link{dpaf}}
#'
#' @return named matrix of hazards
#' @keywords internal
dpaf_hz <- function(z, cf) {
  exp(z %*% -cf)
}

#' Calculate survivals for dpaf study
#'
#' @param hz named matrix of hazards
#' @param ID vector of IDs
#' @param PERIOD vector of periods
#' @param dt delta-times -- length of each period
#'
#' @return named matrix of survivals
#' @keywords internal
dpaf_sv <- function(hz, ID, PERIOD, dt) {
  if (any(tapply(PERIOD, ID, is.unsorted)))
    stop("Periods must be in ascending order for each ID to calculate survival")

  apply(hz, 2, function(hz_col)
    exp(-stats::ave(hz_col, ID, FUN = function(l) cumsum(dt * l)))
  )
}

#' Calculate disease-free survival (**s**ur**v**ival **p**roducts)
#'
#' @param sv named matrix of survivals
#'
#' @return vector of disease-free survivals
#' @keywords internal
dpaf_svp <- function(sv) {
  apply(sv, 1, prod)
}

#' Calculate healthy mortality \eqn{I}
#'
#' @param hz **named** matrix of hazards (must have a column named "disease")
#' @param svp vector of disease-free survivals
#' @param ID vector of IDs
#' @param PERIOD factor vector of PERIODs
#'
#' @return (numeric, scalar) healthy mortality, perhaps to some constant factor
#' @keywords internal

dpaf_i <- function(hz, svp, ID, PERIOD) {
  if (any(tapply(PERIOD, ID, is.unsorted)))
    stop("Periods must be in ascending order for each ID to calculate I")

  svd <- stats::ave(svp, ID, FUN = function(s) -diff(c(1, s)))
  ## (survival at time 0 is 1)
  sum(hz[,"disease"] / rowSums(hz) * svd)
}

#' Calculate gradient of hazard in disease PAF studies
#'
#' @param z model matrix
#' @param hz named matrix of hazards
#'
#' @return a named list of matrices containing hazard gradients
#' @keywords internal
dpaf_ghz <- function(z, hz) {
  lapply(split(hz, col(hz, as.factor = TRUE)), `*`, z)
}

#' Calculate gradient of survival in disease PAF studies
#'
#' @param ghz a named list of hazard gradients
#' @inheritParams dpaf_sv
#'
#' @return a named list of matrices containing survival gradients
#' @keywords internal
dpaf_gsv <- function(ghz, sv, ID, PERIOD, dt) {
  if (any(tapply(PERIOD, ID, is.unsorted)))
    stop("Periods must be in ascending order for each ID")

  sv <- split(sv, col(sv, as.factor = TRUE))
  mapply(
    function(ghz_x, sv_x) -sv_x * apply(
      ghz_x, 2, function(col) stats::ave(
        col, ID, FUN = function(ghz_ij) cumsum(ghz_ij * dt)
      )
    ),
    ghz, sv, SIMPLIFY = FALSE
  )
}

#' Calculate gradient of healthy mortality \eqn{I}
#'
#' @param ghz a named list of hazard gradients
#' @param gsv a named list of survival gradients
#' @param sv a named matrix of survivals
#' @inheritParams dpaf_i
#'
#' @return a named list of vectors for \eqn{\grad{I}} for each of disease and
#'   mortality
#' @keywords internal
dpaf_gi <- function(ghz, gsv, hz, sv, ID, PERIOD) {
  if (any(tapply(PERIOD, ID, is.unsorted)))
    stop("Periods must be in ascending order for each ID")

  svd <- stats::ave(apply(sv, 1, prod), ID, FUN = function(s) -diff(c(1, s)))
  gsvd_d <- apply(
    apply(sv, 1, prod) * gsv$disease, 2,
    # sv[,"mortality"] * gsv$disease, 2,
    function(gsvp_col) stats::ave(gsvp_col, ID, FUN = function(s) -diff(c(0,s)))
  )
  gsvd_m <- apply(
    apply(sv, 1, prod) * gsv$mortality, 2,
    # sv[,"disease"] * gsv$mortality, 2,
    function(gsvp_col) stats::ave(gsvp_col, ID, FUN = function(s) -diff(c(0,s)))
  )

  smnd_d <- hz[,"mortality"] / rowSums(hz)**2 * svd * ghz$disease +
    hz[,"disease"] / rowSums(hz) * gsvd_d
  smnd_m <- -hz[,"disease"] / rowSums(hz)**2 * svd * ghz$mortality +
    hz[,"disease"] / rowSums(hz) * gsvd_m

  lapply(
    list(disease = smnd_d, mortality = smnd_m),
    function(smnd) apply(smnd, 2, sum)
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
