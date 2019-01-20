#' Calculate model frames, parameter estimates and hazard ratios
#'
#' This function fits a survival regression and returns model frames (without
#' response columns) with and without modifications applied, parameter estimates
#' and their standard errors, and hazard ratios with their confidence intervals.
#' Note that hazard ratios reported as \code{NA} are reference levels.
#'
#' \code{modifications} should be a named list of scalars or
#' vectors. Each name should correspond to a column in \code{df}. The
#' first element in each vector should be the desired output value. If
#' there are subsequent elements in the vector, these identify which
#' values to convert from.
#'
#' If \code{covar_model} is not specified,
#' the "hazard ratios" reported include all first order terms.
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
#' @param paf_response a response object as output by \code{\link{gen_data}}
#' @param modifications a list of modifications to apply for PAF calculation;
#'   see Details for more information
#' @param covar_model an optional character vector describing variables of
#'   interest for hazard ratio calculation
#' @param ... extra parameters to be passed to \code{survreg}
#' @param level width of the confidence intervals for hazard ratios, default
#'   0.95.
#'
#' @return a list of type \code{paf_est_matrix} with the following items:
#'
#'   \item{est_matrix_call}{the call to this function}
#'
#'   \item{modifications}{the list of modifications applied}
#'
#'   \item{survreg}{the survival regression object}
#'
#'   \item{coefficients}{the parameter estimates as in \code{survreg}}
#'
#'   \item{var}{the covariance matrix as in \code{survreg}}
#'
#'   \item{terms}{a \code{terms} object for the survival regression with
#'   responses removed}
#'
#'   \item{HR}{the hazard ratios and their confidence intervals; see Warning}
#'
#'   \item{design}{the model frame of predictors}
#'
#'   \item{modified}{the model frame of predictors after modification}
#'
#' @export
est_matrix <- function(sr_formula, paf_response, modifications,
                       covar_model, level = 0.95, ...) {
  stopifnot(inherits(paf_response, c("mpaf_response", "dpaf_response")))

  # store call and modifications for user reference later if desired
  matrix_data <- list(est_matrix_call = match.call(),
                      modifications = modifications)

  # fit the survival regression
  matrix_data$survreg <- survival::survreg(sr_formula, data = paf_response$data,
                                           dist = "exponential", y = FALSE, ...)

  # store the coefficients and variance in parent list, for convenience
  matrix_data$coefficients <- matrix_data$survreg$coefficients
  matrix_data$var <- matrix_data$survreg$var

  # store the terms object, removing response variables -- the result should
  # leave us with something that can create a predictor model.matrix or
  # model.frame
  Terms <- matrix_data$survreg$terms
  if (!inherits(Terms, "terms"))
    # I don't know why this is needed -- it's from survival:::predict.survreg
    stop("invalid terms component of survreg in object")
  Terms <- stats::delete.response(Terms)
  matrix_data$terms <- Terms

  # handle missing covariate model -- get all first-order effects
  if (missing(covar_model)) {
    covar_model <- attr(Terms, "term.labels")[attr(Terms, "order") == 1]
  }

  # select which hazard ratios to calculate and report. This is a mess.
  params <- unlist(lapply(covar_model, function(label) {
    # given a label (like "COVAR" or "factor1:factor2"), determine which factors
    # need to be included
    inc <- which(attr(Terms, "factors")[,label] > 0)

    # get the names of individual factors
    vars <- dimnames(attr(Terms, "factors"))[[1]][inc]

    # get all possible levels for each factor
    vars <- lapply(vars, function(var)
      paste0(var, matrix_data$survreg$xlevels[[var]]))

    # combine "FACTORlevel" names as appropriate, with colon separator
    vars <- do.call(expand.grid, c(vars, KEEP.OUT.ATTRS = FALSE,
                                   stringsAsFactors = FALSE))
    apply(vars, 1, paste, collapse = ":")
  }))

  # show hazard ratios, with reference levels as NA -- note we need to reorder
  # confidence intervals (and names) because we are taking a decreasing
  # transformation
  matrix_data$HR <- cbind(
    "Hazard ratio" = exp(-matrix_data$coefficients)[params],
    exp(-stats::confint(matrix_data$survreg,
                        parm = params,
                        level = level))[,2:1]
  )
  dimnames(matrix_data$HR)[[1]] <- params
  dimnames(matrix_data$HR)[[2]] <- dimnames(matrix_data$HR)[[2]][c(1,3,2)]

  retlist <- c(matrix_data,
               design_frames(paf_response$data, Terms, modifications,
                             matrix_data$survreg$xlevels))
  class(retlist) <- "paf_est_matrix"

  retlist
}
