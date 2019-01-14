#' @export
#' @rdname est_matrix
mpaf_est_matrix <- function(sr_formula, mpaf_data, modifications,
                            covar_model, level = 0.95, ...) {
  .Deprecated("est_matrix")
  est_matrix(sr_formula, mpaf_data, modifications, covar_model,
             level = 0.95, ...)
}

#' Estimate mortality PAFs
#'
#' @param mpaf_fit an object of class \code{paf_est_matrix}
#' @param newdata (not yet implemented) a new prevalence matrix, in some form
#' @param level width of confidence interval, default 0.95
#'
#' @return a list of class \code{mpaf} containing the following elements:
#'
#'   \item{paf, paf0}{the mortality PAFs for individual and cumulative periods,
#'   respectively, with their confidence intervals as specified by the
#'   \code{level} argument}
#'
#'   \item{se_ipaf, se_ipaf0}{the standard errors for \eqn{\log(1-PAF)}, as
#'   above}
#'
#'   \item{grad_paf, grad_paf0}{the gradients of the PAFs as above, for
#'   difference calculations}
#'
#' @export
mpaf_est_paf <- function(mpaf_fit, newdata, level = 0.95) {
  if (!missing(newdata)) {
    stopifnot(inherits(newdata, "mpaf_data"))
    if (!identical(mpaf_fit$data_call$ft_breaks, newdata$data_call$ft_breaks))
      stop("Original periods and new periods are incompatible")
    if (!identical(mpaf_fit$data_call$variables, newdata$data_call$variables))
      stop("Original variables and new variables are different")
    if (!identical(mpaf_fit$data_call$period_factor,
                   newdata$data_call$period_factor))
      stop(paste("Name of period factor columns are not equal.",
                 "Ensure mpaf_gen_data is called with identical",
                 "period_factor arguments"))
  }

  tm <- mpaf_fit$terms
  cf <- mpaf_fit$coefficients
  vv <- mpaf_fit$var
  dt <- diff(mpaf_fit$breaks)

  # Point estimate calculations ---------------------------------------------

  if (missing(newdata)) {
    z <-      stats::model.matrix(tm, mpaf_fit$design)
    z_star <- stats::model.matrix(tm, mpaf_fit$modified)
    ID <- mpaf_fit$ID
    PERIOD <- mpaf_fit$PERIOD
  } else {
    newframes <- design_frames(newdata$data, tm, mpaf_fit$modifications,
                               mpaf_fit$survreg$xlevels)
    z <-      stats::model.matrix(tm, newframes$design)
    z_star <- stats::model.matrix(tm, newframes$modified)
    ID <- newdata$ID
    PERIOD <- newdata$PERIOD
  }

  lambda <-      mpaf_lambda(z,      cf)
  lambda_star <- mpaf_lambda(z_star, cf)

  S <-      mpaf_S(lambda,      ID, PERIOD, dt)
  S_star <- mpaf_S(lambda_star, ID, PERIOD, dt)

  dS <-      mpaf_delta_S(S,      ID, PERIOD)
  dS_star <- mpaf_delta_S(S_star, ID, PERIOD)

  # I_(0, t]^(*)
  I_0 <-      mpaf_I(1-S,      PERIOD)
  I_0_star <- mpaf_I(1-S_star, PERIOD)
  names(I_0) <- names(I_0_star) <- paste0(
    "(", mpaf_fit$breaks[1], ",", utils::tail(mpaf_fit$breaks, -1), "]"
  )

  # I_(t, t+dt]^(*)
  I <-      mpaf_I(dS,      PERIOD)
  I_star <- mpaf_I(dS_star, PERIOD)

  ipaf0 <- log(I_0_star) - log(I_0)
  ipaf <- log(I_star) - log(I)

  # Gradient calculations ---------------------------------------------------

  grad_lambda <-      mpaf_grad_lambda(z,      lambda)
  grad_lambda_star <- mpaf_grad_lambda(z_star, lambda_star)

  grad_S <-      mpaf_grad_S(grad_lambda,      S,      ID, PERIOD, dt)
  grad_S_star <- mpaf_grad_S(grad_lambda_star, S_star, ID, PERIOD, dt)

  grad_dS <-      mpaf_grad_delta_S(grad_S,      ID, PERIOD)
  grad_dS_star <- mpaf_grad_delta_S(grad_S_star, ID, PERIOD)

  grad_I <-      mpaf_grad_I(grad_dS,      ID, PERIOD)
  grad_I_star <- mpaf_grad_I(grad_dS_star, ID, PERIOD)

  grad_I_0 <-      mpaf_grad_I(-grad_S,      ID, PERIOD)
  grad_I_0_star <- mpaf_grad_I(-grad_S_star, ID, PERIOD)
  dimnames(grad_I_0)[[1]] <- dimnames(grad_I_0_star)[[1]] <- paste0(
    "(", mpaf_fit$breaks[1], ",", utils::tail(mpaf_fit$breaks, -1), "]"
  )

  grad_ipaf0 <- grad_I_0_star / I_0_star - grad_I_0 / I_0
  grad_ipaf <- grad_I_star / I_star - grad_I / I

  var_ipaf0 <- grad_ipaf0 %*% vv %*% t(grad_ipaf0)
  var_ipaf <- grad_ipaf %*% vv %*% t(grad_ipaf)

  # standard error and confint calculations ---------------------------------

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
    grad_paf0 = (grad_I_0 * I_0_star - grad_I_0_star * I_0) / I_0 ^ 2,
    grad_paf = (grad_I * I_star - grad_I_star * I) / I ^ 2
  )
}

#' Estimate significance of groupwise differences in mortality PAFs
#'
#' The elements of the output from this function should be able to be passed to
#' \code{\link[stats]{printCoefmat}}.
#'
#' @param mpaf1,mpaf2 objects of class \code{mpaf} (from
#'   \code{\link{mpaf_est_paf}}) to be compared
#' @param vv the covariance matrix of parameters
#'
#' @return a matrix containing the PAF differences, their standard errors, Z
#'   values and p values.
#' @export
mpaf_est_diff <- function(mpaf1, mpaf2, vv) {
  dpaf <- c(mpaf1$paf[,1] - mpaf2$paf[,1], mpaf1$paf0[,1] - mpaf2$paf0[,1])

  gdpaf <- rbind(mpaf1$grad_paf - mpaf2$grad_paf,
                 mpaf1$grad_paf0 - mpaf2$grad_paf0)

  se_dpaf <- sqrt(diag(gdpaf %*% vv %*% t(gdpaf)))
  Z <- dpaf / se_dpaf

  cbind(
    "PAF Diff" = dpaf,
    "SE(PAF Diff)" = se_dpaf,
    "Z value" = Z,
    "Pr(>|Z|)" = 2 * stats::pnorm(-abs(Z))
  )
}
