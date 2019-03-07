#' Calculate mortality PAF estimates, confidence intervals and groupwise
#' differences
#'
#' @inheritParams est_matrix
#' @param paf_data an object containing the response as in \code{\link{gen_data}}
#' @param prevalence_data matrix of prevalences as class \code{paf_data}
#' @param group_vars names of columns to group by in PAF difference estimations
#' @param level width of confidence intervals, default 0.95
#' @param ... extra parameters to be passed to \code{\link{est_matrix}}
#'
#' @return a list of potenitally interesting hazard ratio and PAF estimates
#' @export
mpaf_summary <- function(sr_formula, paf_data, modifications, covar_model,
                         prevalence_data, group_vars, level = 0.95, ...) {
  # fit model using est_matrix -- see est_matrix.R
  mpaf_fit <- est_matrix(sr_formula, paf_data, modifications,
                         covar_model, level = level, ...)

  # calculate PAFs, possibly using new prevalence data -- see mortality_est.R
  if (missing(prevalence_data))
    mpaf_all <- mpaf_est_paf(mpaf_fit, paf_data, level = level)
  else
    mpaf_all <- mpaf_est_paf(mpaf_fit, paf_data, prevalence_data,
                             level = level)

  # begin groupwise PAF estimates
  mpaf_groups <- NULL
  mpaf_diffs <- NULL
  if (!missing(group_vars)) {
    pafdat <- if (missing(prevalence_data)) paf_data else prevalence_data

    # split PAF data -- see paf_utils.R
    # include drop = FALSE to keep useful names
    pd_spl <- paf_data_split(pafdat, pafdat$data[, group_vars, drop = FALSE])

    # call mpaf_est_paf for each subgroup -- see mortality_est.R again
    mpaf_groups <- lapply(pd_spl, mpaf_est_paf,
                          mpaf_fit = mpaf_fit, level = level)
    names(mpaf_groups) <- names(pd_spl)

    # call mpaf_est_diff for every pair of PAF estimates -- see mortality_est.R
    mpaf_diffs <- utils::combn(mpaf_groups, 2, FUN = function(grp_lst)
      mpaf_est_diff(grp_lst[[1]], grp_lst[[2]], mpaf_fit$var),
      simplify = FALSE
    )
    names(mpaf_diffs) <-
      utils::combn(names(pd_spl), 2, paste, collapse = " - ")
  }

  retlist <- c(list(call = match.call()), mpaf_fit, mpaf_all,
               list("nobs" = paf_data$nobs, "na.action" = paf_data$na.action),
               group_pafs = list(mpaf_groups),
               paf_diffs = list(mpaf_diffs))
  class(retlist) <- "mpaf_summary"
  retlist
}
