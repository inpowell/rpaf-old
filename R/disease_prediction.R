#### Mid-level functions to calculate quantities of interest
# - `predict.gpaf` (for PAFs and HRs)
# - `dpaf_groups` (for group differences)


#' Predicting PAFs in disease-death cohort studies
#'
#' This function currently only allows for the calculation of PAFs for disease,
#' with their standard errors and confidence intervals, as well as hazard ratios.
#'
#' \code{modlist} should be a named list of scalars or vectors. Each name should
#' correspond to a column in the object's data frame. The first element in each
#' vector should be the desired output value. If there are subsequent elements
#' in the vector, these identify which values to convert from.
#'
#' The \code{type} argument can be:
#'
#' \itemize{
#'
#' \item{\code{'paf'}, where this function will calculate the population
#' attributable fraction, and other quantities requested}
#'
#' \item{\code{'hr'}, where this function will calculate the hazard ratios. In
#' this case, arguments \code{modlist}, \code{newdata}, and \code{gradient} are
#' ignored.}
#'
#' }
#'
#' @param object an object of class \code{\link{dpaf}}
#' @param modlist an named list of replacement vectors (see Details)
#' @param newdata an object of class \code{\link{dpaf_data}}
#' @param type the quantity to be estimated
#' @param confint if \code{TRUE}, calculate confidence intervals
#' @param level confidence level for interval calculations, default 0.95
#' @param se.trans if \code{TRUE}, calculate standard errors in the space
#'   assumed to be normal
#' @param gradient if \code{TRUE}, calculate the gradient of the PAF and return
#'   as a list
#' @param ... other arguments
#'
#' @return A named vector or matrix containing the prediction, its standard
#'   error in normal space, and its confidence interval, when the latter two are
#'   requested.
#'
#'   If \code{gradient = TRUE}, then the return value will be a named list
#'   containing the above vector (in \code{estimate}) and gradients with respect
#'   to disease and mortality coefficients (in \code{disease} and
#'   \code{mortality}, respectively)
#'
#' @export
#' @method predict dpaf
predict.dpaf <- function(object, modlist, newdata,
                         type = c('paf', 'hr'),
                         confint = TRUE, level = 0.95,
                         se.trans = FALSE, gradient = FALSE,
                         ...) {
  type <- match.arg(type)

  cf <- object$coefficients
  vv_l <- object$var

  if (type == 'paf') {
    Terms <- object$survreg_d$terms
    if (!inherits(Terms, "terms"))
      # I don't know why this is needed -- it's from survival:::predict.survreg
      stop("invalid terms component of survreg in object")
    Terms <- stats::delete.response(Terms)

    if (missing(newdata))
      dpdta <- object
    else
      dpdta <- newdata

    mod_df <- apply_modifications(dpdta$data, modlist)
    raw_frame <- stats::model.frame(Terms, data = dpdta$data,
                                    na.action = stats::na.pass,
                                    xlev = object$survreg_d$xlevels)
    mod_frame <- stats::model.frame(Terms, data = mod_df,
                                    na.action = stats::na.pass,
                                    xlev = object$survreg_d$xlevels)

    z_raw <- stats::model.matrix(Terms, raw_frame)
    z_mod <- stats::model.matrix(Terms, mod_frame)

    hz_raw <- dpaf_hz(z_raw, cf)
    sv_raw <- dpaf_sv(hz_raw, dpdta$ID, dpdta$PERIOD, diff(dpdta$breaks))
    svp_raw <- dpaf_svp(sv_raw)
    i_raw <- dpaf_i(hz_raw, svp_raw, dpdta$ID, dpdta$PERIOD)

    hz_mod <- dpaf_hz(z_mod, cf)
    sv_mod <- dpaf_sv(hz_mod, dpdta$ID, dpdta$PERIOD, diff(dpdta$breaks))
    svp_mod <- dpaf_svp(sv_mod)
    i_mod <- dpaf_i(hz_mod, svp_mod, dpdta$ID, dpdta$PERIOD)

    if (confint || se.trans || gradient) {
      ghz_raw <- dpaf_ghz(z_raw, hz_raw)
      gsv_raw <- dpaf_gsv(ghz_raw, sv_raw, dpdta$ID, dpdta$PERIOD, diff(dpdta$breaks))
      gi_raw <- dpaf_gi(ghz_raw, gsv_raw, hz_raw, sv_raw, dpdta$ID, dpdta$PERIOD)

      ghz_mod <- dpaf_ghz(z_mod, hz_mod)
      gsv_mod <- dpaf_gsv(ghz_mod, sv_mod, dpdta$ID, dpdta$PERIOD, diff(dpdta$breaks))
      gi_mod <- dpaf_gi(ghz_mod, gsv_mod, hz_mod, sv_mod, dpdta$ID, dpdta$PERIOD)

      gipaf <- dpaf_gipaf(gi_mod, i_mod, gi_raw, i_raw)

      v_ipaf <- sum(mapply(function(grad, vv) t(grad) %*% vv %*% grad,
                           gipaf, vv_l))
    }

    pred <- c("PAF" = 1 - i_mod / i_raw)

    if (se.trans)
      pred <- c(pred, "SE iPAF"=sqrt(v_ipaf))
    if (confint) {
      halpha <- (1 - level) / 2
      ipaf <- log(i_mod) - log(i_raw)
      ci <- ipaf + sqrt(v_ipaf) *
        stats::qnorm(c("lwr" = 1 - halpha, "upr" = halpha))
      ci <- -expm1(ci) # convert to PAF space
      pred <- c(pred, ci)
    }
    if (gradient)
      pred <- c("estimate" = list(pred), dpaf_gpaf(gi_mod, i_mod, gi_raw, i_raw))
  }
  if (type == 'hr') {
    pred <- split(-cf, col(cf, as.factor = TRUE)) # log hazard ratios
    names(pred$disease) <- row.names(cf)
    names(pred$mortality) <- row.names(cf)

    # calculations in normal space
    if (confint || se.trans)
      se <- lapply(vv_l, function(vv) sqrt(diag(vv)))
    if (confint) {
      halpha <- (1 - level) / 2
      ci <- mapply(
        function(lhr, se) lhr + se %o%
          stats::qnorm(c("lwr" = halpha, "upr" = 1 - halpha)),
        pred, se, SIMPLIFY = FALSE
      )
    }

    # transformations to HR space
    pred <- lapply(pred, exp)
    if (confint)
      ci <- lapply(ci, exp)

    # return data
    pred <- lapply(pred, function(hr) cbind("Hazard Ratio" = hr))
    if (se.trans)
      pred <- mapply(function(pred, se) cbind(pred, "SE(logHR)" = se),
                     pred, se, SIMPLIFY = FALSE)
    if (confint)
      pred <- mapply(cbind, pred, ci, SIMPLIFY = FALSE)
  }

  return(pred)
}

#' Analyse groupwise PAF estimates
#'
#' @param object an object of class \code{dpaf}
#' @param modlist a list of modifications to be passed to
#'   \code{\link{predict.dpaf}}
#' @param split_data a list of objects of class \code{dpaf_data}, from
#'   \code{\link{dpaf_split}}
#' @param gdiffs if \code{TRUE}, analyse the differences in PAFs between groups
#' @param ... other arguments to be passed to \code{predict.dpaf}
#'
#' @return a list of up to two elements
#'
#'   \item{group_pafs}{a matrix of calculated disease PAFs, as well as other
#'   information output by \code{predict.dpaf}}
#'
#'   \item{paf_diffs}{(optional) a matrix of differences between PAFs, with
#'   their standard errors, test statistics, and p-values}
#'
#' @export
dpaf_groups <- function(object, modlist, split_data, gdiffs = TRUE, ...) {
  stopifnot(inherits(object, "dpaf"))
  stopifnot(all(sapply(split_data, inherits, "dpaf_split")))

  gwise_preds <- lapply(split_data, stats::predict, object = object,
                        modlist = modlist, type = 'paf', gradient = gdiffs, ...)

  if (gdiffs) {
    pafs <- do.call(rbind, lapply(gwise_preds, `[[`, "estimate"))
    diffs <- utils::combn(gwise_preds, 2, FUN = function(preds) {
      pred1 <- preds[[1]]
      pred2 <- preds[[2]]

      dpaf <- pred1$estimate["PAF"] - pred2$estimate["PAF"]
      names(dpaf) <- NULL

      dgrad <- list(
        disease = pred1$disease - pred2$disease,
        mortality = pred1$mortality - pred2$mortality
      )

      se_dpaf <- sqrt(sum(mapply(function(g, vv) t(g) %*% vv %*% g,
                                 dgrad, object$var)))

      z <- dpaf / se_dpaf
      p <- 2 * stats::pnorm(-abs(z))
      c("PAF Diff" = dpaf, "SE(PAF Diff)" = se_dpaf, "Z value" = z, "Pr(>|Z|)" = p)
    }, simplify = FALSE)

    names(diffs) <- utils::combn(names(gwise_preds), 2, paste, collapse = " - ")
  } else {
    pafs <- do.call(rbind, gwise_preds)
  }
  if (gdiffs)
    list(
      group_pafs = pafs,
      paf_diffs = do.call(rbind, diffs)
    )
  else
    list(group_pafs = pafs)
}
