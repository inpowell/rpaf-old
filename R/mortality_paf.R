mpaf_summary <- function(sr_formula, mpaf_data, modifications, covar_model,
                         prevalence_data, group_var, level = 0.95, ...) {
  mpaf_fit <- mpaf_est_matrix(sr_formula, mpaf_data, modifications,
                              covar_model, level = level, ...)

  if (missing(prevalence_data))
    mpaf_all <- mpaf_est_paf(mpaf_fit, level = level)
  else
    mpaf_all <- mpaf_est_paf(mpaf_fit, prevalence_data, level = level)

  mpaf_groups <- list()
  if (!missing(group_var)) {

  }
}
