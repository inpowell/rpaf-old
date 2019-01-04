#' Calculate model frames, parameter estimates and hazard ratios
#'
#' This function fits a survival regression and returns model frames (without
#' response columns) with and without modifications applied, parameter estimates
#' and their standard errors, and hazard ratios with their confidence intervals.
#' Note that hazard ratios reported as \code{NA} are reference levels.
#'
#' If \code{covar_model} is not specified, the "hazard ratios" reported include
#' all first order terms.
#'
#' @section Warnings:
#'
#'   A hazard ratio is not well-defined when the term of interest is an
#'   interaction. Use caution when interpreting such results.
#'
#'   When using interaction terms in \code{covar_model}, the wrong order of
#'   covariate names will result in an error. Only terms that appear in
#'   \code{attr(terms(sr_formula), "term.labels")} can be used.
#'
#' @param sr_formula a formula to be passed to \code{\link[survival]{survreg}}
#' @param mpaf_data a response object as output by \code{\link{mpaf_gen_data}}
#' @param modifications a list of modifications to apply for PAF calculation;
#'   see Details for more information
#' @param covar_model an optional character vector describing variables of
#'   interest for hazard ratio calculation
#' @param ... extra parameters to be passed to \code{survreg} and \code{confint}
#'
#' @return a list of type \code{mpaf_est_matrix} with the following items:
#'
#'   \item{est_matrix_call}{the call to this function}
#'
#'   \item{survreg}{the survival regression object}
#'
#'   \item{coefficients}{the parameter estimates as in \code{survreg}}
#'
#'   \item{var}{the covariance matrix as in \code{survreg}}
#'
#'   \item{HR}{the hazard ratios and their confidence intervals; see Warning}
#'
#'   \item{design}{the model frame of predictors}
#'
#'   \item{modified}{the model frame of predictors after modification}
#'
#' @export
mpaf_est_matrix <- function(sr_formula, mpaf_data, modifications,
                            covar_model, level = 0.95, ...) {
  matrix_data <- list(est_matrix_call = match.call(),
                      modifications = modifications)

  # fit the survival regression
  matrix_data$survreg <- survival::survreg(sr_formula, data = mpaf_data$data,
                                           dist = "exponential", y = FALSE, ...)

  matrix_data$coefficients <- matrix_data$survreg$coefficients
  matrix_data$var <- matrix_data$survreg$var

  # terms object for use later
  Terms <- matrix_data$survreg$terms
  if (!inherits(Terms, "terms"))
    # I don't know why this is needed -- it's from survival:::predict.survreg
    stop("invalid terms component of survreg in object")
  Terms <- stats::delete.response(Terms)
  matrix_data$terms <- Terms

  # handle missing covariate model
  if (missing(covar_model)) {
    covar_model <- attr(Terms, "term.labels")[attr(Terms, "order") == 1]
  }

  # select which hazard ratios to calculate and report. This is a mess.
  params <- unlist(lapply(covar_model, function(label) {
    inc <- which(attr(Terms, "factors")[,label] > 0)
    vars <- dimnames(attr(Terms, "factors"))[[1]][inc]
    vars <- lapply(vars, function(var)
      paste0(var, matrix_data$survreg$xlevels[[var]]))
    vars <- do.call(expand.grid, c(vars, KEEP.OUT.ATTRS = FALSE,
                                   stringsAsFactors = FALSE))
    apply(vars, 1, paste, collapse = ":")
  }))

  # show hazard ratios, with reference levels as NA
  matrix_data$HR <- cbind(
    "Hazard ratio" = exp(-matrix_data$coefficients)[params],
    exp(-confint(matrix_data$survreg, parm = params, level = level))[,2:1]
  )
  dimnames(matrix_data$HR)[[1]] <- params
  dimnames(matrix_data$HR)[[2]] <- dimnames(matrix_data$HR)[[2]][c(1,3,2)]

  matrix_data$design <- stats::model.frame(Terms, data = mpaf_data$data,
                                           na.action = na.pass,
                                           xlev = matrix_data$survreg$xlevels)
  mod_df <- apply_modifications(mpaf_data$data, modifications)
  matrix_data$modified <- stats::model.frame(Terms, data = mod_df,
                                             na.action = na.pass,
                                             xlev = matrix_data$survreg$xlevels)

  retlist <- c(mpaf_data, matrix_data)
  class(retlist) <- "mpaf_est_matrix"

  retlist
}

#' Estimate mortality PAFs
#'
#' @param mpaf_fit an object of class \code{mpaf_est_matrix}
#' @param newdata (not yet implemented) a new prevalence matrix, in some form
#'
#' @return
#' @export
mpaf_est_paf <- function(mpaf_fit, newdata) {
  if (!missing(newdata))
    # New data might need to be implemented in mpaf_est_matrix?
    stop("Use of new data is not yet implemented.")

  tm <- mpaf_fit$terms
  cf <- mpaf_fit$coefficients
  vv <- mpaf_fit$var
  ID <- mpaf_fit$ID
  PERIOD <- mpaf_fit$PERIOD
  dt <- diff(mpaf_fit$breaks)

  z <-      stats::model.matrix(tm, mpaf_fit$design)
  z_star <- stats::model.matrix(tm, mpaf_fit$modified)

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
    "(", mpaf_fit$breaks[1], ",", tail(mpaf_fit$breaks, -1), "]"
  )

  # I_(t, t+dt]^(*)
  I <-      mpaf_I(dS,      PERIOD)
  I_star <- mpaf_I(dS_star, PERIOD)

  ipaf0 <- log(I_0_star) - log(I_0)
  ipaf <- log(I_star) - log(I)

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
    "(", mpaf_fit$breaks[1], ",", tail(mpaf_fit$breaks, -1), "]"
  )

  grad_ipaf0 <- grad_I_0_star / I_0_star - grad_I_0 / I_0
  grad_ipaf <- grad_I_star / I_star - grad_I / I

  var_ipaf0 <- grad_ipaf0 %*% vv %*% t(grad_ipaf0)
  var_ipaf <- grad_ipaf %*% vv %*% t(grad_ipaf)

  list(
    lambda = lambda,
    S = S,
    I = I,
    grad_lambda = grad_lambda,
    grad_S = grad_S,
    grad_I = grad_I,
    ipaf0 = ipaf0,
    var_ipaf0 = var_ipaf0,
    ipaf = ipaf,
    var_ipaf = var_ipaf
  )
}
