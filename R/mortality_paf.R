#' Calculate mortality PAF estimates, confidence intervals and groupwise
#' differences
#'
#' @inheritParams mpaf_est_matrix
#' @param prevalence_data matrix of prevalences as class \code{mpaf_data}
#' @param group_vars names of columns to group by in PAF difference estimations
#' @param level width of confidence intervals, default 0.95
#' @param ... extra parameters to be passed to \code{\link{mpaf_est_matrix}}
#'
#' @return a list of potenitally interesting hazard ratio and PAF estimates
#' @export
mpaf_summary <- function(sr_formula, mpaf_data, modifications, covar_model,
                         prevalence_data, group_vars, level = 0.95, ...) {
  mpaf_fit <- mpaf_est_matrix(sr_formula, mpaf_data, modifications,
                              covar_model, level = level, ...)

  if (missing(prevalence_data))
    mpaf_all <- mpaf_est_paf(mpaf_fit, level = level)
  else
    mpaf_all <- mpaf_est_paf(mpaf_fit, prevalence_data, level = level)

  mpaf_groups <- NULL
  mpaf_diffs <- NULL
  if (!missing(group_vars)) {
    pafdat <- if (missing(prevalence_data)) mpaf_data else prevalence_data
    pd_spl <- paf_data_split(pafdat, pafdat$data[, group_vars, drop = FALSE])
    mpaf_groups <- lapply(pd_spl, mpaf_est_paf,
                          mpaf_fit = mpaf_fit, level = level)
    names(mpaf_groups) <- names(pd_spl)

    mpaf_diffs <- utils::combn(mpaf_groups, 2, FUN = function(grp_lst)
      mpaf_est_diff(grp_lst[[1]], grp_lst[[2]], mpaf_fit$var),
      simplify = FALSE
    )
    names(mpaf_diffs) <-
      utils::combn(names(pd_spl), 2, paste, collapse = " - ")
  }

  retlist <- c(list(call = match.call()), mpaf_fit, mpaf_all,
               group_pafs = list(mpaf_groups),
               paf_diffs = list(mpaf_diffs))
  class(retlist) <- "mpaf_summary"
  retlist
}
