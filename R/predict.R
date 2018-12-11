#' Predict mortality PAFs and Hazard ratios
#'
#' This function allows for users to estimate various survival statistics (e.g.
#' hazard ratios) from an object of class \code{mpafreg}, as output by
#' \code{\link{mpafreg}}. Current options available are PAFs and
#' hazard ratios, and their transformations in normal space. Options exist to
#' also calculate asymptotic confidence intervals for each output.
#'
#' This function gives various outputs depending on the \code{type} input. If
#' \code{type} is \code{'paf'}, then the estimation of interest is the PAF. If
#' \code{type} is \code{'l1mpaf'}, then the variable of interest is \eqn{\log(1
#' - PAF)}, which is assumed to be asymptotically normal distributed. Standard
#' errors and variances requested for PAFs are delivered in this space.
#' Similarly, when \code{type} is \code{'hr'} or \code{'lhr'}, the outputs are
#' respectively the hazard ratios (sometimes akin to the relative risk) and
#' log-hazard ratios, the latter of which are assumed to be normal.
#'
#' \code{modlist} should be a named list of scalars or vectors. Each name should
#' correspond to a column in the object's data frame. The first element in each
#' vector should be the desired output value. If there are subsequent elements
#' in the vector, these identify which values to convert from.
#'
#' @param object An object of class \code{mpafreg}
#' @param modlist A named list of modifications used in PAF calculation. See
#'   Details for more information. Ignored if calculating hazard ratios.
#' @param newdata An optional data frame of predictors to estimate PAFs from.
#'   Ignored if calculating hazard ratios. If needed and not supplied, the
#'   \code{model} data frame from the \code{object} will be used.
#' @param type Type of prediction. See Details for options.
#' @param confint Include a confidence interval?
#' @param level Confidence level for interval calculations
#' @param se.trans Return standard error in normal (transformed) space?
#' @param vcov.trans Return variance-covariance matrix in normal (transformed)
#'   space?
#' @param na.action function determining what to do with missing values in
#'   \code{newdata}. Default is to predict NA.
#' @param ... other arguments (ignored)
#'
#' @return \code{predict.survreg} produces either a vector of predictions, or a
#'   matrix of predictions and confidence limits with names \code{fit},
#'   \code{lwr}, and \code{upr} (as for \code{\link[stats]{predict.lm}}). If
#'   either \code{se.trans} or \code{vcov.trans} are \code{TRUE}, then return a
#'   named list, with names "prediction" (for above), "se.trans" (for the
#'   standard error), and "vcov.trans" (for the variance matrix).
#'
#'   If \code{type = 'gpaf'}, then a further list element named
#'   \code{"gradient"} is included, containing an array of gradients for each
#'   time point.
#'
#' @method predict mpafreg
#' @export
#'
predict.mpafreg <- function(object, modlist, newdata,
                            type = c('paf', 'l1mpaf', 'hr', 'lhr', 'gpaf'),
                            confint = FALSE, level = 0.95,
                            se.trans = FALSE, vcov.trans = FALSE,
                            na.action = stats::na.pass, ...) {
  type <- match.arg(type)

  cf <- object$survreg$coefficients
  vv <- object$survreg$var

  # hazard ratios -----------------------------------------------------------
  if (type %in% c('lhr', 'hr')) {
    # do we need the standard error?
    if (confint || se.trans)
      se <- sqrt(diag(vv))

    pred <- -cf

    if (confint)
      pred <- cbind(
        "fit" = pred,
        "lwr" = pred + se * stats::qnorm((1 - level) / 2, lower.tail = TRUE),
        "upr" = pred + se * stats::qnorm((1 - level) / 2, lower.tail = FALSE)
      )

    # transform, if required
    if (type == 'hr')
      pred <- exp(pred)
  }


  # PAFs --------------------------------------------------------------------
  if (type %in% c('paf', 'l1mpaf', 'gpaf')) {
    # figure out the terms of interest
    Terms <- object$survreg$terms
    if (!inherits(Terms, "terms"))
      # I don't know why this is needed -- it's from survival:::predict.survreg
      stop("invalid terms component of survreg in object")
    Terms <- stats::delete.response(Terms)

    if (missing(newdata))
      df <- object$model
    else
      df <- newdata

    mod_df <- apply_modifications(df, modlist)
    raw_frame <- stats::model.frame(Terms, data = df, na.action = na.action,
                                    xlev = object$survreg$xlevels)
    mod_frame <- stats::model.frame(Terms, data = mod_df, na.action = na.action,
                                    xlev = object$survreg$xlevels)

    na.action.raw <- attr(raw_frame, "na.action")
    na.action.mod <- attr(mod_frame, "na.action")
    # error checking on missing value treatment
    stopifnot(identical(na.action.mod, na.action.raw))

    if (length(na.action.raw)) {
      ID <- df[-na.action.raw, object$id_var]
      PERIOD <- df[-na.action.raw, object$period_var]
    } else {
      ID <- df[,object$id_var]
      PERIOD <- df[,object$period_var]
    }
    x_raw <- stats::model.matrix(Terms, raw_frame)
    x_mod <- stats::model.matrix(Terms, mod_frame)

    # calculate survival statistics
    ss_raw <- survival_statistics(x_raw, cf, ID, PERIOD,
                                  object$ft_length, object$ft_delta,
                                  gradient = confint || se.trans ||
                                    vcov.trans || type == 'gpaf')
    ss_mod <- survival_statistics(x_mod, cf, ID, PERIOD,
                                  object$ft_length, object$ft_delta,
                                  gradient = confint || se.trans ||
                                    vcov.trans || type == 'gpaf')

    # calculate standard errors, if required
    if (confint || se.trans || vcov.trans) {
      Dlm <- ss_mod$grad_lmortality - ss_raw$grad_lmortality
      vv <- Dlm %*% vv %*% t(Dlm)
    }
    if (confint || se.trans)
      se <- sqrt(diag(vv))

    # calculate transformed PAFs
    pred <- log(ss_mod$mortality / ss_raw$mortality)
    if (confint) {
      pred <- cbind(
        pred,
        pred + se * stats::qnorm((1 - level) / 2, lower.tail = TRUE),
        pred + se * stats::qnorm((1 - level) / 2, lower.tail = FALSE)
      )
      # need to reverse upper and lower CLs if transforming to PAF space
      if (type %in% c('paf', 'gpaf)'))
        pred[,c(2,3)] <- pred[,c(3,2)]
      colnames(pred) <- c("fit", "lwr", "upr")
    }

    # transform PAFs, if required
    if (type %in% c('paf', 'gpaf'))
      pred <- -expm1(pred)
  }


  # PAFs with gradient ------------------------------------------------------

  if (type == 'gpaf') {
    grad <- ss_raw$grad_mortality * ss_mod$mortality -
      ss_mod$grad_mortality * ss_raw$mortality
    grad <- grad / (ss_raw$mortality ** 2)
  }

  # return values -----------------------------------------------------------
  if (confint)
    attr(pred, "level") <- level

  if (se.trans || vcov.trans || type == 'gpaf') {
    retlist <- list("prediction" = pred)
    if (se.trans)
      retlist <- append(retlist, list("se.trans" = se))
    if (vcov.trans)
      retlist <- append(retlist, list("vcov.trans" = vv))
    if (type == 'gpaf')
      retlist <- append(retlist, list("gradient" = grad))
    retlist
  } else {
    pred
  }
}



