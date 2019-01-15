dpaf_est_paf <- function(fit_d, fit_m, paf_data, newdata, level = 0.95) {
  if (!missing(newdata)) {
    stopifnot(inherits(newdata, "paf_data"))
    if (!identical(mpaf_data$data_call$ft_breaks, newdata$data_call$ft_breaks))
      stop("Original periods and new periods are incompatible")
    if (!identical(mpaf_data$data_call$variables, newdata$data_call$variables))
      stop("Original variables and new variables are different")
    if (!identical(mpaf_data$data_call$period_factor,
                   newdata$data_call$period_factor))
      stop(paste("Name of period factor columns are not equal.",
                 "Ensure mpaf_gen_data is called with identical",
                 "period_factor arguments"))
  }

  # ensure design frames are the same for each fit
  ## TEST this -- the correct answer may occur even if the following conditions
  ## don't hold (below is sufficient, but not necessary)
  stopifnot(identical(fit_d$terms, fit_m$terms))
  stopifnot(identical(fit_d$design, fit_m$design))
  stopifnot(identical(fit_d$modified, fit_m$modified))
  stopifnot(identical(fit_d$modifications, fit_m$modifications))

  tm <- fit_d$terms
  cf <- list("disease" = fit_d$coefficients,
             "mortality" = fit_m$coefficients)
  vv <- list("disease" = fit_d$var,
             "mortality" = fit_m$var)
  dt <- diff(paf_data$breaks)

  if (missing(newdata)) {
    z <-      stats::model.matrix(tm, fit_m$design)
    z_star <- stats::model.matrix(tm, fit_m$modified)
    ID <- paf_data$ID
    PERIOD <- paf_data$PERIOD
  } else {
    newframes <- design_frames(newdata$data, tm, fit_m$modifications,
                               fit_m$survreg$xlevels)
    z <-      stats::model.matrix(tm, newframes$design)
    z_star <- stats::model.matrix(tm, newframes$modified)
    ID <- newdata$ID
    PERIOD <- newdata$PERIOD
  }

}
