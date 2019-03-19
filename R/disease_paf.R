#' Calculate mortality PAF estimates, confidence intervals and groupwise
#' differences
#'
#'
#'
#' If \code{hr_out} is not specified, the "hazard ratios" reported include
#' all first order terms.
#'
#' \code{modifications} should be a named list of scalars or vectors. Each name
#' should correspond to a column in \code{df}. The first element in each vector
#' should be the desired output value. If there are subsequent elements in the
#' vector, these identify which values to convert from.
#' @param primary_resp formula with response for primary outcome regression (note RHS
#'   should be \code{.})
#' @param secondary_resp formula with response for secondary outcome regression
#'   (note RHS should be \code{.})
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
dpaf_summary <- function(primary_resp, secondary_resp, predictors, dpaf_data,
                         modifications, hr_out, prevalence_data,
                         group_vars, level = 0.95, ...) {
  # construct survival regression formulae
  primary_fm <- stats::update(primary_resp, predictors)
  secondary_fm <- stats::update(secondary_resp, predictors)

  # fit disease and mortality models using est_matrix -- see
  # est_matrix.R
  fit1 <- est_matrix(primary_fm, dpaf_data, modifications, hr_out,
                      level = level, ...)
  fit2 <- est_matrix(secondary_fm, dpaf_data, modifications, hr_out,
                      level = level, ...)

  vv_l <- list("primary" = fit1$var, "secondary" = fit2$var)

  # calculate PAFs, possibly using new prevalence data -- see
  # disease_est.R
  if (missing(prevalence_data))
    dpaf_all <- dpaf_est_paf(fit1, fit2, dpaf_data, level = level)
  else
    dpaf_all <- dpaf_est_paf(fit1, fit2, dpaf_data, prevalence_data,
                             level = level)
  # begin groupwise PAF estimates
  paf_groups <- NULL
  paf_diffs <- NULL
  if (!missing(group_vars)) {
    pafdat <- if (missing(prevalence_data)) dpaf_data else prevalence_data

    # split PAF data -- see paf_utils.R
    # add drop = FALSE to keep useful names
    pd_spl <- paf_data_split(pafdat, pafdat$data[, group_vars, drop = FALSE])

    # call dpaf_est_paf for each subgroup -- again see disease_est.R
    paf_groups <- lapply(pd_spl, dpaf_est_paf,
                         fit1 = fit1, fit2 = fit2,
                         paf_data = dpaf_data,
                         level = level)
    names(paf_groups) <- names(pd_spl)

    # call dpaf_est_diff for each pair of PAF estimates -- see disease_est.R
    paf_diffs <- utils::combn(paf_groups, 2, FUN = function(grp_lst)
      dpaf_est_diff(grp_lst[[1]], grp_lst[[2]], vv_l),
      simplify = FALSE
    )
    names(paf_diffs) <-
      utils::combn(names(pd_spl), 2, paste, collapse = " - ")
  }

  # prevent name collisions between fit1 and fit2
  names(fit1) <- paste0(names(fit1), "_1")
  names(fit2) <- paste0(names(fit2), "_2")

  retlist <- c(list(call = match.call()),
               list("nobs" = dpaf_data$nobs, "na.action" = dpaf_data$na.action),
               fit1, fit2, dpaf_all,
               group_pafs = list(paf_groups),
               paf_diffs = list(paf_diffs))
  class(retlist) <- "dpaf_summary"
  retlist
}
