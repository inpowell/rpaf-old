#' Fit all survival regression models for PAF data
#'
#' This function fits a survival regression object for every event specified in
#' a PAF data frame. It returns a list object with names determined by the
#' column names of the events. The primary event will always occur first in the
#' list.
#'
#' @param predictors a formula with empty LHS containing the predictors of interest
#' @param paf_data an extended data frame as from \code{\link{gen_data_fun}}
#' @param ... further objects to be passed to \code{\link[survival]{survreg}}
#'
#' @return a list of \code{survreg} objects from the \code{survival} package
#' @export
#'
#' @example exec/minifhs-prep.R
fit_models <- function(predictors, paf_data, ...) {
  primary <- attr(paf_data, "primary")
  secondaries <- attr(paf_data, "secondary")
  time <- attr(paf_data, "time")

  # response columns
  model_cols <- list(primary, secondaries)
  models <- lapply(model_cols, function(col) {
    lhs <- survival::Surv(paf_data[[time]], paf_data[[col]])
    fm <- stats::update(lhs ~ ., predictors)

    survival::survreg(fm, data = paf_data, dist = "exponential", ...)
  })

  structure(models, names = model_cols)
}
