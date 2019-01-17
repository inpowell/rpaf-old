#' Estimate disease PAFs
#'
#' @param fit_d the \code{est_matrix} object fitted to disease incidence
#' @param fit_m the \code{est_matrix} object fitted to mortality
#' @param paf_data an object of class \code{paf_data}, used for the
#'   \code{\link{est_matrix}} call
#' @param newdata new prevalences, as an object of class \code{paf_data}
#' @param level width of confidence interval, default 0.95
#'
#' @return a list containing the following elements:
#'
#'   \item{paf, paf0}{the PAFs for disease over individual and cumulative
#'   periods, respectively, with their confidence intervals as specified by the
#'   \code{level} argument}
#'
#'   \item{se_ipaf, se_ipaf0}{the standard errors for \eqn{\log(1-PAF)}, as
#'   above}
#'
#'   \item{grad_paf, grad_paf0}{the gradients of the PAFs as above, for
#'   difference calculations}
#'
#' @export
dpaf_est_paf <- function(fit_d, fit_m, paf_data, newdata, level = 0.95) {
  if (!missing(newdata)) {
    # equivalence assertions
    stopifnot(inherits(newdata, "paf_data"))
    if (!identical(paf_data$data_call$ft_breaks, newdata$data_call$ft_breaks))
      stop("Original periods and new periods are incompatible")
    if (!identical(paf_data$data_call$variables, newdata$data_call$variables))
      stop("Original variables and new variables are different")
    if (!identical(paf_data$data_call$period_factor,
                   newdata$data_call$period_factor))
      stop(paste("Name of period factor columns are not equal.",
                 "Ensure gen_data is called with identical",
                 "period_factor arguments"))
    # setting ID and PERIOD
    ID <- newdata$ID
    PERIOD <- newdata$PERIOD
  } else {
    # setting ID and PERIOD
    ID <- paf_data$ID
    PERIOD <- paf_data$PERIOD
  }
  if (is_period_unsorted(ID, PERIOD))
    stop("Periods must be in ascending order by ID for PAF calculations")

  # ensure design frames are the same for each fit
  ## TEST this -- the correct answer may occur even if the following conditions
  ## don't hold (below is sufficient, but not necessary)
  stopifnot(isTRUE(all.equal(fit_d$terms, fit_m$terms)))
  stopifnot(identical(fit_d$design, fit_m$design))
  stopifnot(identical(fit_d$modified, fit_m$modified))
  stopifnot(identical(fit_d$modifications, fit_m$modifications))

  tm <- fit_d$terms
  cf_l <- list("disease" = fit_d$coefficients,
               "mortality" = fit_m$coefficients)
  vv_l <- list("disease" = fit_d$var,
               "mortality" = fit_m$var)
  dt <- diff(paf_data$breaks)

  # Point estimate calculations ---------------------------------------------

  if (missing(newdata)) {
    z <-      stats::model.matrix(tm, fit_m$design)
    z_star <- stats::model.matrix(tm, fit_m$modified)
  } else {
    newframes <- design_frames(newdata$data, tm, fit_m$modifications,
                               fit_m$survreg$xlevels)
    z <-      stats::model.matrix(tm, newframes$design)
    z_star <- stats::model.matrix(tm, newframes$modified)
  }

  lambda <-      dpaf_lambda(z,      cf_l)
  lambda_star <- dpaf_lambda(z_star, cf_l)

  S <-      dpaf_S(lambda,      ID, PERIOD, dt)
  S_star <- dpaf_S(lambda_star, ID, PERIOD, dt)

  Sp <-      dpaf_Sp(S)
  Sp_star <- dpaf_Sp(S_star)

  dSp <-      dpaf_delta_Sp(Sp,      ID, PERIOD)
  dSp_star <- dpaf_delta_Sp(Sp_star, ID, PERIOD)

  # I_(t, t+dt]^(*)
  I <-      dpaf_I(lambda,      dSp,      PERIOD)
  I_star <- dpaf_I(lambda_star, dSp_star, PERIOD)
  names(I) <- names(I_star) <- paste0(
    '(', utils::head(paf_data$breaks, -1), ',',
    utils::tail(paf_data$breaks, -1), ']'
  )

  # I_(0, t]^(*)
  I_0 <-      cumsum(I)
  I_0_star <- cumsum(I_star)
  names(I_0) <- names(I_0_star) <- paste0(
    "(", paf_data$breaks[1], ",", utils::tail(paf_data$breaks, -1), "]"
  )

  ipaf0 <- log(I_0_star) - log(I_0)
  ipaf <- log(I_star) - log(I)

  # Gradient calculations ---------------------------------------------------

  grad_lambda <-      dpaf_grad_lambda(z,      lambda)
  grad_lambda_star <- dpaf_grad_lambda(z_star, lambda_star)

  grad_S <-      dpaf_grad_S(grad_lambda,      S,      ID, PERIOD, dt)
  grad_S_star <- dpaf_grad_S(grad_lambda_star, S_star, ID, PERIOD, dt)


  grad_dS <-      dpaf_grad_delta_Sp(grad_S,      Sp,      ID, PERIOD)
  grad_dS_star <- dpaf_grad_delta_Sp(grad_S_star, Sp_star, ID, PERIOD)

  grad_I <-      dpaf_grad_I(grad_lambda,      grad_dS,      lambda,
                             dSp,      PERIOD)
  grad_I_star <- dpaf_grad_I(grad_lambda_star, grad_dS_star, lambda_star,
                             dSp_star, PERIOD)

  # need gradients to be row vectors before applying cumsum if there is only one
  # period
  if (length(levels(PERIOD)) == 1) {
    grad_I <- lapply(grad_I, rbind)
    grad_I_star <- lapply(grad_I_star, rbind)
  }

  grad_I_0 <-      lapply(grad_I,      function(gI) apply(gI, 2, cumsum))
  grad_I_0_star <- lapply(grad_I_star, function(gI) apply(gI, 2, cumsum))

  # also need these to be row vectors because magic
  if (length(levels(PERIOD)) == 1) {
    grad_I_0 <- lapply(grad_I, rbind)
    grad_I_0_star <- lapply(grad_I_star, rbind)
  }

  grad_ipaf <- dpaf_gipaf(grad_I_star, I_star, grad_I, I)
  grad_ipaf0 <- dpaf_gipaf(grad_I_0_star, I_0_star, grad_I_0, I_0)

  grad_paf <- dpaf_gpaf(grad_I_star, I_star, grad_I, I)
  grad_paf0 <- dpaf_gpaf(grad_I_0_star, I_0_star, grad_I_0, I_0)

  var_ipaf0 <- do.call(`+`, mapply(function(gipaf0, vv)
    gipaf0 %*% vv %*% t(gipaf0), grad_ipaf0, vv_l, SIMPLIFY = FALSE))
  var_ipaf <- do.call(`+`, mapply(function(gipaf, vv)
    gipaf %*% vv %*% t(gipaf), grad_ipaf, vv_l, SIMPLIFY = FALSE))

  # Standard error and confint calculations ---------------------------------

  a <- (1 - level)/2
  a <- c(a, 1 - a)

  se_ipaf0 <- sqrt(diag(var_ipaf0))
  ci0 <- ipaf0 + se_ipaf0 %o% stats::qnorm(1 - a)
  colnames(ci0) <- paste(format(100*a, trim = TRUE, scientific = FALSE,
                                digits = 3), "%")
  paf0 <- -expm1(cbind("PAF" = ipaf0, ci0))

  se_ipaf <- sqrt(diag(var_ipaf))
  ci <- ipaf + se_ipaf %o% stats::qnorm(1 - a)
  colnames(ci) <- paste(format(100*a, trim = TRUE, scientific = FALSE,
                               digits = 3), "%")
  paf <- -expm1(cbind("PAF" = ipaf, ci))

  list(
    paf0 = paf0,
    paf = paf,
    se_ipaf0 = se_ipaf0,
    se_ipaf = se_ipaf,
    grad_paf0 = grad_paf0,
    grad_paf = grad_paf
  )
}

#' Estimate significance of groupwise differences in disease PAFs
#'
#' The elements of the output from this function should be able to be passed to
#' \code{\link[stats]{printCoefmat}}.
#'
#' @param dpaf1,dpaf2 objects from \code{\link{dpaf_est_paf}} to be compared
#' @param vv_l named list of the covariance matrices of parameters for disease and
#'   mortality
#'
#' @return a matrix containing the PAF differences, their standard errors, Z
#'   values and p values.
#' @export
dpaf_est_diff <- function(dpaf1, dpaf2, vv_l) {
  paf_diff <- c(dpaf1$paf[,1] - dpaf2$paf[,1], dpaf1$paf0[,1] - dpaf2$paf0[,1])

  gpaf_diff <- mapply(
    rbind,
    mapply(`-`, dpaf1$grad_paf, dpaf2$grad_paf, SIMPLIFY = FALSE),
    mapply(`-`, dpaf1$grad_paf0, dpaf2$grad_paf0, SIMPLIFY = FALSE),
    SIMPLIFY = FALSE
  )

  se_paf_diff <- sqrt(
    do.call(`+`, mapply(function(gpd, vv) diag(gpd %*% vv %*% t(gpd)),
                        gpaf_diff, vv_l, SIMPLIFY = FALSE))
  )
  Z <- paf_diff / se_paf_diff

  cbind(
    "PAF Diff" = paf_diff,
    "SE(PAF Diff)" = se_paf_diff,
    "Z value" = Z,
    "Pr(>|Z|)" = 2 * stats::pnorm(-abs(Z))
  )
}
