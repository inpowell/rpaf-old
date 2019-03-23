#' Calculate hazard ratios from a survival regression object
#'
#' This function reports hazard ratios from a survival analysis where (possibly
#' piecewise) constant hazards are assumed. Where interaction terms are present,
#' the function automatically determines which groups to calculate hazard ratios
#' for. Note that approximate confidence intervals based off the asymptotic
#' normality of log-hazard ratios are reported.
#'
#' @param survreg a \code{\link[survival]{survreg}} object of the model to base
#'   hazard ratios upon
#' @param hr_out a list of variables to report hazard ratios for
#' @param level the desired widtch of the confidence interval, default 0.95
#'
#' @return A named nested list of hazard ratio tables. The top level list
#'   separates different risk factors requested in \code{hr_out}. The nested
#'   lists distinguish groups across which hazard ratios may differ.
#' @export
hazard_ratios <- function(survreg, hr_out, level = 0.95) {
  # coefficient vector and covariance matrix
  cf <- stats::coef(survreg)
  vv <- stats::vcov(survreg)

  # metastructural variables for survival regression
  tt <- stats::terms(survreg)
  Terms <- stats::delete.response(tt)
  ttf <- attr(Terms, "factors")

  # named list of variables with all possible desired values
  var.list <- lapply(rownames(ttf), function(var) {
    if (var %in% names(survreg$xlevels))
      survreg$xlevels[[var]] # factor, with levels given in xlevels
    else
      0:1 # numeric (hopefully), with just a unit change wanted
  })
  names(var.list) <- rownames(ttf)

  # all data combinations, including irrelevant ones
  comb.dat <- expand.grid(var.list)
  MM <- stats::model.matrix(survreg, data = comb.dat)
  AA <- attr(MM, "assign")

  # calculating hazard ratios for each given risk factor
  hr_list <- lapply(hr_out, function(label) {
    # variables which are relevant to the requested factor (catches interaction
    # terms), as rows of the terms factor table, and as a character vector
    fr <- apply(ttf[, ttf[label, ] > 0, drop = FALSE], 1, any)
    vars <- rownames(ttf)[fr]

    # corresponding columns in the model matrix
    mmc <- ttf[label, AA, drop = FALSE] > 0

    # account for possible intercept term
    if (AA[1] == 0)
      mmc <- c(FALSE, mmc)

    # remove duplicate rows
    dup <- duplicated(comb.dat[,vars])
    mm <- sweep(MM[!dup,], MARGIN = 2, mmc, `*`)

    # want to determine sensible groups to calculate hazard ratios
    gf <- comb.dat[!dup, setdiff(vars, label), drop = FALSE]

    # No groups present --> make a dummy "ALL" group for data-structural
    # consistency
    if (ncol(gf) == 0)
      gf <- data.frame(ALL = gl(1, nrow(gf), labels = ""))

    # row names of the hazard ratio table
    rn <- if (label %in% names(survreg$xlevels))
      survreg$xlevels[[label]]
    else
      c("Base", "Unit increase")

    # calculate hazard ratios for each group
    by(as.data.frame(mm), INDICES = gf, FUN = hr_group,
       cf = cf, vv = vv, rn = rn, level = level, simplify = FALSE)
  })
  names(hr_list) <- hr_out

  hr_list
}

#' Calculate hazard ratios within a particular group
#'
#' @param mm model matrix fragment of interest
#' @param cf coefficient vector of survival regression object
#' @param vv variance matrix from survival regression object
#' @param rn rownames of output matrix
#' @param level confidence level desired
#'
#' @return matrix of hazard ratios and confidence intervals
hr_group <- function(mm, cf, vv, rn, level) {
  # `by` coerces model matrix to a data frame
  mm <- as.matrix(mm)

  # log hazard ratios
  hr <- as.vector(mm %*% cf)

  # standard errors
  se <- sqrt(diag(mm %*% vv %*% t(mm)))

  # confints
  a <- (1 - level) / 2
  a <- c(1 - a, a)
  hr <- cbind(hr, hr + se %o% stats::qnorm(a))
  colnames(hr) <- c("Point", paste(
    format(100*rev(a), trim = TRUE, scientific = FALSE, digits = 3), "%"
  ))
  rownames(hr) <- rn

  # true hazard ratios
  exp(-hr)
}
