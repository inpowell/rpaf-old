#' Calculate mortality PAF estimates, confidence intervals and groupwise
#' differences
#'
#'
#'
#' If \code{covar_model} is not specified,
#' the "hazard ratios" reported include all first order terms.
#'
#' \code{modifications} should be a named list of scalars or vectors. Each name
#' should correspond to a column in \code{df}. The first element in each vector
#' should be the desired output value. If there are subsequent elements in the
#' vector, these identify which values to convert from.
#' @param disease_resp formula with response for disease regression (note RHS
#'   should be \code{.})
#' @param death_resp formula with response for mortality regression (note RHS
#'   should be \code{.})
#' @param predictors formula with RHS containing predictors for both survival
#'   regressions (LHS should be empty)
#' @param dpaf_data an object of class \code{dpaf_response}, from
#'   \code{\link{gen_data}}
#' @inheritParams est_matrix
#' @param prevalence_data optional matrix of external prevalences as class
#'   \code{paf_data}
#' @param group_vars names of columns to group by in PAF difference estimations
#' @param level width of confidence intervals, default 0.95
#' @param ... extra parameters to be passed to \code{\link{est_matrix}}
#'
#' @return a list of potentially interesting hazard ratio and PAF estimates
#' @export
dpaf_summary <- function(disease_resp, death_resp, predictors, dpaf_data,
                         modifications, covar_model, prevalence_data,
                         group_vars, level = 0.95, ...) {
  # construct survival regression formulae
  disease_fm <- stats::update(disease_resp, predictors)
  mortality_fm <- stats::update(death_resp, predictors)

  fit_d <- est_matrix(disease_fm, dpaf_data, modifications, covar_model,
                      level = level, ...)
  fit_m <- est_matrix(mortality_fm, dpaf_data, modifications, covar_model,
                      level = level, ...)

  vv_l <- list("disease" = fit_d$var, "mortality" = fit_m$var)

  if (missing(prevalence_data))
    dpaf_all <- dpaf_est_paf(fit_d, fit_m, dpaf_data, level = level)
  else
    dpaf_all <- dpaf_est_paf(fit_d, fit_m, dpaf_data, prevalence_data,
                             level = level)

  paf_groups <- NULL
  paf_diffs <- NULL
  if (!missing(group_vars)) {
    pafdat <- if (missing(prevalence_data)) dpaf_data else prevalence_data
    pd_spl <- paf_data_split(pafdat, pafdat$data[, group_vars, drop = FALSE])
    paf_groups <- lapply(pd_spl, dpaf_est_paf,
                         fit_d = fit_d, fit_m = fit_m,
                         paf_data = dpaf_data,
                         level = level)
    names(paf_groups) <- names(pd_spl)

    paf_diffs <- utils::combn(paf_groups, 2, FUN = function(grp_lst)
      dpaf_est_diff(grp_lst[[1]], grp_lst[[2]], vv_l),
      simplify = FALSE
    )
    names(paf_diffs) <-
      utils::combn(names(pd_spl), 2, paste, collapse = " - ")
  }

  names(fit_d) <- paste0(names(fit_d), "_d")
  names(fit_m) <- paste0(names(fit_m), "_m")

  retlist <- c(list(call = match.call()),
               list("nobs" = dpaf_data$nobs, "na.action" = dpaf_data$na.action),
               fit_d, fit_m, dpaf_all,
               group_pafs = list(paf_groups),
               paf_diffs = list(paf_diffs))
  class(retlist) <- "dpaf_summary"
  retlist
}
