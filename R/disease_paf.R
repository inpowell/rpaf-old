#' Prepare data for disease PAF calculations
#'
#' @param rawdata data frame from cohort study
#' @param variables character vector of columns of interest
#' @param id_var name of column containing unique individual IDs
#' @param mortality_time name of column for time of death or censoring
#' @param mortality_ind name of logical column for death (assumed \code{TRUE} if
#'   death occurred during follow-up)
#' @param disease_time name of column for time of disease or censoring
#' @param disease_ind name of logical column for disease (assumed \code{TRUE} if
#'   disease occurred during follow-up)
#' @param ft_length numeric length of follow-up time
#' @param ft_delta numeric length of constant hazard period
#' @param period_start_time
#' @param period_factor
#' @param period_censor_time name for new column describing time
#' @param na.action
#'
#' @return
#' @export
#'
#' @examples
dpaf_data <- function(rawdata, variables, id_var,
                      mortality_time, mortality_ind,
                      disease_time, disease_ind,
                      ft_length, ft_delta,
                      period_start_time = "p_start",
                      period_factor = "period",
                      period_censor_time = "f_end",
                      na.action = getOption("na.action")) {
  # deal with missing values
  data <- rawdata[,c(id_var, mortality_time, mortality_ind,
                     disease_time, disease_ind, variables)]

  data <- match.fun(na.action)(data)

  # disease must always occur before death
  data[data[,disease_ind],mortality_ind] <- FALSE

  data[,period_censor_time] <-
    with(data, pmin(get(mortality_time), get(disease_time)))

  # handle time periods
  n_fperiod <- ft_length %/% ft_delta
  ## (note we're considering right-closed intervals
  if (ft_length %% ft_delta == 0)
    n_fperiod <- n_fperiod - 1

  period_df <- data.frame((0:n_fperiod) * ft_delta)
  names(period_df) <- period_start_time

  long_data <- merge(period_df, data, by = NULL)

  # determine if period is before, during, or after period where censoring takes
  # place
  post_censor <- with(
    long_data,
    get(period_censor_time) < get(period_start_time)
  )
  pre_censor <- with(
    long_data,
    get(period_censor_time) > ft_delta + get(period_start_time)
  )
  during_censor <- !post_censor & !pre_censor

  # modify death and disease indicators wrt relationship with period
  long_data[post_censor, c(period_censor_time, mortality_ind, disease_ind)] <-
    NA
  long_data[pre_censor, period_censor_time] <- ft_delta
  long_data[pre_censor, c(mortality_ind, disease_ind)] <- FALSE
  long_data[c(during_censor), period_censor_time] <-
    with(long_data[during_censor,],
         get(period_censor_time) - get(period_start_time))

  # label periods appropriately
  long_data[,period_factor] <- with(
    long_data,
    factor(get(period_start_time), levels = period_df[, period_start_time],
           ordered = TRUE)
  )
  if (length(levels(long_data[,period_factor])) > 1)
    contrasts(long_data[,period_factor]) <- "contr.treatment"

  levels(long_data[,period_factor]) <-
    paste0("(", levels(long_data[,period_factor]), ",",
           as.numeric(levels(long_data[,period_factor])) + ft_delta, "]")

  # return the data
  list(data =
         long_data[,c(
           id_var,
           period_censor_time, # time of death, disease, or censoring
           disease_ind, # did death occur first?
           mortality_ind, # did disease occur first?
           period_factor, period_start_time, # factor, numeric representation of period
           variables # everything else
         )],
       nobs = nrow(data),
       na.action = na.action(data),
       post.censor = which(post_censor)
  )
}

#' Internal fitting function to calculate PAFs for disease censored by death
#'
#' Note that the survivals output have been adjusted for for a smaller final
#' period \eqn{((J-1) \mathtt{ft_delta}, \mathtt{ft_length}]}.
#'
#' @param x model matrix
#' @param cf_d \code{survreg} coefficients for disease
#' @param cf_m \code{survreg} coefficients for mortality
#' @param id numeric vector of individual ID's, corresponding to rows of
#'   \code{x}
#' @param period factor describing period of interest
#' @param breaks vector of times of breaks -- should have length one greater
#'   than the number of periods
#'
#' @return
#' @export
#' @keywords internal
#'
#' @section Warning: this function assumes that \code{period} is internally
#'   sorted by \code{id} (for \code{cumsum} to work as expected).
#'
dpaf_fit <- function(x, cf_d, cf_m,
                     id, period, breaks) {

  unsorted <- any(unlist(
    tapply(period, id, is.unsorted)
  ))

  if (unsorted) {
    warning(paste("Model matrix is not correctly sorted.",
                  "`dpaf_fit` will sort, and results will",
                  "not match id and period arguments given."))
    ord <- order(period)
    x <- x[ord,]
    id <- id[ord]
    period <- period[ord]
  }

  # calculate hazards
  hz_d <- as.vector(exp(x %*% -cf_d))
  hz_m <- as.vector(exp(x %*% -cf_m))

  # calculate time differences between periods
  deltas <- diff(breaks)

  # calculate survivals
  sv_d <- ave(hz_d, id, FUN = function(hz) exp(-cumsum(deltas * hz)))
  sv_m <- ave(hz_m, id, FUN = function(hz) exp(-cumsum(deltas * hz)))

  # cumulative disease-free survival
  sv <- sv_d * sv_m

  # differences in each survival type
  dsv <- ave(sv, id, FUN = function(svl) -diff(c(1, svl)))
  ##    (note survival at time zero is 1)

  # corrected survivals diffs for disease morbidity
  cdsv <- dsv * hz_d / (hz_d + hz_m)

  # morbidity -- mean of disease survivals for each time period
  morb <- mean(tapply(cdsv, id, sum))

  list(
    "hazard.disease" = hz_d,
    "hazard.mortality" = hz_m,
    "survival.disease" = sv_d,
    "survival.mortality" = sv_m,
    "morbidity" = morb
  )
}

#' Fit PAFs for disease data
#'
#' @param disease_resp formula with response for disease regression (note RHS
#'   should be \code{.})
#' @param death_resp formula with response for mortality regression (note RHS
#'   should be \code{.})
#' @param predictors formula with RHS containing predictors for both survival
#'   regressions (LHS should be empty)
#' @param modlist list of modifications for PAF calculation (see
#'   \code{\link{apply_modifications}} and \code{\link{mpafreg}} for more
#'   information)
#' @param data_list a list containing data frames and other information, from
#'   \code{\link{dpaf_data}}
#' @param ft_length length of follow-up time
#' @param ft_delta length of periods of constant hazard
#' @param id_var name of column containing individual ID, as for
#'   \code{\link{dpaf_data}}
#' @param period_var name of column containing period factor, as for
#'   \code{\link{dpaf_data}}
#' @param ... other arguments to be passed to \code{\link[survival]{survreg}}.
#'
#' @return
#' @export
#'
#' @examples
dpaf <- function(disease_resp, death_resp, predictors,
                 modlist, data_list, ft_length, ft_delta,
                 id_var, period_var, ...) {
  formula_d <- update(predictors, disease_resp)
  formula_m <- update(predictors, death_resp)

  if (length(data_list$post.censor) != 0)
    mf <- data_list$data[-data_list$post.censor,]
  else
    mf <- data_list$data

  survreg_d <- survival::survreg(formula_d, mf, dist = "exponential", ...)
  survreg_m <- survival::survreg(formula_m, mf, dist = "exponential", ...)

  x_pred <- model.matrix(predictors, data_list$data, na.action = na.fail,
                         xlev = survreg_d$xlevels)
  x_mdfd <- model.matrix(predictors,
                         apply_modifications(data_list$data, modlist),
                         na.action = na.fail, xlev = survreg_d$xlevels)

  breaks <- seq(from = 0, to = ft_length, by = ft_delta)
  if (tail(breaks, 1) != ft_length)
    breaks <- c(breaks, ft_length)

  cf_d <- survreg_d$coefficients
  cf_m <- survreg_m$coefficients
  ID <- data_list$data[[id_var]]
  PERIOD <- data_list$data[[period_var]]

  fit_pred <- dpaf_fit(x_pred, cf_d, cf_m, ID, PERIOD, breaks)
  fit_mdfd <- dpaf_fit(x_mdfd, cf_d, cf_m, ID, PERIOD, breaks)

  retlist <- list(
    call = match.call(),
    nobs = data_list$nobs,
    na.action = na.action(data_list$data),
    model = model.frame(predictors, data_list$data),
    disease.survreg = survreg_d,
    mortality.survreg = survreg_m,
    modlist = modlist,
    breaks = breaks,
    ID = ID,
    PERIOD = PERIOD,
    fitted.raw = fit_pred,
    fitted.modified = fit_mdfd,
    paf = 1 - fit_mdfd$morbidity / fit_pred$morbidity
  )

  class(retlist) <- "dpaf"

  retlist
}

#' Calculate gradients of morbidity for PAF calculations
#'
#' @param x
#' @param hz_d
#' @param hz_m
#' @param sv_d
#' @param sv_m
#' @param breaks
#' @param ID
#' @param PERIOD
#'
#' @return
#' @export
#' @keywords internal
dpaf_grad <- function(x, hz_d, hz_m, sv_d, sv_m, breaks, ID, PERIOD) {
  ghz_d <- hz_d * x # yes, row-wise multiplication
  ghz_m <- hz_m * x

  gsv_d <- sv_d * apply(ghz_d, 2, function(ghz_col)
    ave(ghz_col, ID, FUN = function(ghz_i) - cumsum(diff(breaks) * ghz_i))
  )
  gsv_m <- sv_m * apply(ghz_m, 2, function(ghz_col)
    ave(ghz_col, ID, FUN = function(ghz_i) - cumsum(diff(breaks) * ghz_i))
  )

  # calculate differences for ease of later use
  ## (note survival is 1 at t=0 and gradient is zero at t=0)
  dsv <- ave(sv_d * sv_m, ID, FUN = function(sv) -diff(c(1, sv)))
  gdsv_d <- apply(sv_d * sv_m * gsv_d, 2, function(gsv_col)
    ave(gsv_col, ID, FUN = function(gsv) -diff(c(0, gsv)))
  )
  gdsv_m <- apply(sv_d * sv_m * gsv_m, 2, function(gsv_col)
    ave(gsv_col, ID, FUN = function(gsv) -diff(c(0, gsv)))
  )

  gmorb_d <- apply(
    ghz_d * hz_m / (hz_d + hz_m)**2 * dsv + hz_d / (hz_d + hz_m) * gdsv_d,
    2, function(col) mean(tapply(col, ID, sum))
  )
  gmorb_m <- apply(
    -ghz_m * hz_d / (hz_d + hz_m)**2 * dsv + hz_d / (hz_d + hz_m) * gdsv_m,
    2, function(col) mean(tapply(col, ID, sum))
  )

  list(
    hazard.disease = ghz_d,
    hazard.mortality = ghz_m,
    survival.disease = gsv_d,
    survival.mortality = gsv_m,
    morbidity.disease = gmorb_d,
    morbidity.mortality = gmorb_m
  )
}

#' Title
#'
#' @param object
#' @param parm
#' @param level
#' @param ...
#'
#' @return
#' @export
#' @method confint dpaf
#' @keywords internal
#'
#' @examples
confint.dpaf <- function(object, parm, level = 0.95, ...) {
  Terms = object$disease.survreg$terms
  Terms <- stats::delete.response(Terms)

  x <- stats::model.matrix(Terms, object$model,
                           xlev = object$disease.survreg$xlevels)
  x_mod <- stats::model.matrix(Terms, apply_modifications(object$model, object$modlist),
                               xlev = object$disease.survreg$xlevels)

  grad <- with(object$fitted.raw,
               dpaf_grad(x, hazard.disease, hazard.mortality,
                         survival.disease, survival.mortality,
                         object$breaks, object$ID, object$PERIOD))
  grad_mod <- with(object$fitted.modified,
                   dpaf_grad(x_mod, hazard.disease, hazard.mortality,
                             survival.disease, survival.mortality,
                             object$breaks, object$ID, object$PERIOD))

  glmorb_d <- grad$morbidity.disease / object$fitted.raw$morbidity
  glmorb_m <- grad$morbidity.mortality / object$fitted.raw$morbidity
  glmorb_mod_d <- grad_mod$morbidity.disease / object$fitted.modified$morbidity
  glmorb_mod_m <- grad_mod$morbidity.mortality / object$fitted.modified$morbidity

  vv_d <- stats::vcov(object$disease.survreg)
  vv_m <- stats::vcov(object$mortality.survreg)

  ipaf <- log1p(-object$paf)

  gipaf_d <- glmorb_mod_d - glmorb_d
  gipaf_m <- glmorb_mod_m - glmorb_m

  var_ipaf <- diag(t(gipaf_d) %*% vv_d %*% gipaf_d +
                     t(gipaf_m) %*% vv_m %*% gipaf_m)

  alpha <- (1 - level) / 2
  ci <- ipaf + sqrt(var_ipaf) * qnorm(c("lwr" = 1 - alpha, "upr" = alpha),
                                      lower.tail = TRUE)
  ci <- -expm1(ci)
  c(ci, "se" = sqrt(var_ipaf))
}
