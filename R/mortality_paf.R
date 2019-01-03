#' Calculate model frames, parameter estimates and hazard ratios (with
#' confidence intervals)
#'
#' This function fits a survival regression and returns model frames (without
#' response) with and without modifications applied.
#'
#' If \code{covar_model} is not specified, the "hazard ratios" reported include
#' all terms, including interactions.
#'
#' @section Warning: A hazard ratio is not well-defined when the term of
#'   interest is an interaction. Use caution when interpreting such results.
#'
#' @param sr_formula a formula to be passed to \code{\link[survival]{survreg}}
#' @param mpaf_data a response object as output by \code{\link{mpaf_gen_data}}
#' @param modifications a list of modifications to apply for PAF calculation;
#'   see Details for more information
#' @param covar_model an optional character vector describing variables
#'   of interest for hazard ratio calculation
#' @param ... extra parameters to be passed to \code{survreg} and \code{confint}
#'
#' @return a list of type \code{mpaf_est_matrix} with the following items:
#'
#' \item{est_matrix_call}{the call to this function}
#'
#' \item{survreg}{the survival regression object}
#'
#' \item{coefficients}{the parameter estimates as in \code{survreg}}
#'
#' \item{var}{the covariance matrix as in \code{survreg}}
#'
#' \item{HR}{the hazard ratios and their confidence intervals; see Warning}
#'
#' \item{design}{the model frame of predictors}
#'
#' \item{modified}{the model frame of predictors after modification}
#'
#' @export
mpaf_est_matrix <- function(sr_formula, mpaf_data, modifications,
                            covar_model, ...) {
  matrix_data <- list(est_matrix_call = match.call())

  matrix_data$survreg <- survival::survreg(sr_formula, data = mpaf_data$data,
                                           dist = "exponential", y = FALSE, ...)

  matrix_data$coefficients <- matrix_data$survreg$coefficients
  matrix_data$var <- matrix_data$survreg$var

  if (missing(covar_model)) {
    covar_model <- variable.names(matrix_data$survreg)
  }
  matrix_data$HR <- cbind(
    "Hazard ratio" = exp(-matrix_data$coefficients),
    exp(-confint(matrix_data$survreg, parm = covar_model, ...))[,2:1]
  )
  dimnames(matrix_data$HR)[[2]] <- dimnames(matrix_data$HR)[[2]][c(1,3,2)]

  Terms <- sr_formula$survreg$terms
  if (!inherits(Terms, "terms"))
    # I don't know why this is needed -- it's from survival:::predict.survreg
    stop("invalid terms component of survreg in object")
  Terms <- stats::delete.response(Terms)

  matrix_data$design <- stats::model.frame(Terms, data = mpaf_data$data,
                                           na.action = na.pass,
                                           xlev = matrix_data$survreg$xlevels)
  mod_df <- apply_modifications(mpaf_data$data, modifications)
  matrix_data$modified <- stats::model.frame(Terms, data = mod_df,
                                             na.action = na.pass,
                                             xlev = matrix_data$survreg$xlevels)

  matrix_data
}
