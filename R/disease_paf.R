
# Data preparation utilities ----------------------------------------------

#' Prepare data for PAF disease studies
#'
#' @param rawdata a data frame containing
#' @param variables character vector of variables of interest
#' @param time_breaks numeric vector of breaks between periods of constant
#'   hazard
#' @param id_var column of data frame identifying individuals (will be generated
#'   sequentially if missing or \code{NULL})
#' @param death_time column name, describing time of death or censoring
#' @param death_ind column name, which is \code{TRUE} if individual died during
#'   follow up
#' @param disease_time column name, describing time of disease or censoring
#' @param disease_ind column name, which is \code{TRUE} if individual contracted
#'   disease during follow up
#' @param period_factor new column name, describing the time period for baseline
#'   hazard calculation
#' @param ft_end new column name, describing the time of death, disease or
#'   censoring from the start of each period
#' @param na.action a function (found using \code{match.fun}) used to treat
#'   missing data
#' @param include_period logical, which is \code{TRUE} if the output data should
#'   contain a factor column for period
#'
#' @return a list of class \code{dpaf_data} containing at least the following
#'   elements:
#'
#'   \item{data_call}{The function call to generate data in this form}
#'
#'   \item{nobs}{The number of observations used from the initial data frame}
#'
#'   \item{na.action}{The IDs removed from the data due to missingness}
#'
#'   \item{breaks}{A vector of times where baseline hazard changes}
#'
#'   \item{ID}{A vector of IDs corresponding to rows of the data frame}
#'
#'   \item{PERIOD}{A factor vector of periods corresponding to rows of the data
#'   frame}
#'
#'   \item{data}{A data frame containing predictor columns as given in the
#'   \code{variables} argument, as well as a column for the period (if
#'   requested)}
#'
#'   If all of \code{death_time}, \code{death_ind}, \code{disease_time} and
#'   \code{disease_ind} are specified, then the data frame has response columns
#'   describing censoring status and follow-up time end, which can be passed to
#'   survival regression functions. In this case, the returned value has
#'   additional class \code{dpaf_response}.
#'
#' @export
#'
dpaf_data <- function(rawdata, variables, time_breaks, id_var,
                      death_time, death_ind,
                      disease_time, disease_ind,
                      period_factor = "period",
                      ft_end = "f_end",
                      na.action = getOption("na.action"),
                      include_period = TRUE) {
  retlist <- list(data_call = match.call())

  # ensure breaks are ordered
  if (is.unsorted(time_breaks))
    stop("Time breaks are not in order.")

  if (missing(id_var) || is.null(id_var)) {
    id_var <- "ID"
    rawdata[[id_var]] <- 1:nrow(rawdata)
  }

  # select columns
  if (missing(death_time) || missing(death_ind) ||
      missing(disease_time) || missing(disease_ind))
    data <- rawdata[,c(id_var, variables), drop = FALSE]
  else
    data <- rawdata[,c(id_var, variables, death_time, death_ind,
                       disease_time, disease_ind)]

  # handle missing values
  data <- match.fun(na.action)(data)
  retlist <- c(retlist, list(nobs = nrow(data), na.action = na.action(data)))


  # add in time breaks
  retlist <- c(retlist, list("breaks" = time_breaks))

  # lengthen data -----------------------------------------------------------
  # integer codes as in findInterval -- note not a factor yet
  intervals <- data.frame(1:(length(time_breaks) - 1))
  names(intervals) <- period_factor

  data <- merge(intervals, data, by = NULL)

  # correct indicators and times
  if (!missing(death_time) && !missing(death_ind) &&
      !missing(disease_time) && !missing(disease_ind)) {
    data[[ft_end]] <- pmin(data[[death_time]], data[[disease_time]])
    # note that breaks[data[[period_factor]]] gives the time at start of period
    # and breaks[-1][data[[period_factor]]] gives time at end of period
    ft_end_int <- findInterval(data[[ft_end]], time_breaks, left.open = TRUE)
    post_end <- ft_end_int < data[[period_factor]]
    pre_end <- ft_end_int > data[[period_factor]]
    is_end <- ft_end_int == data[[period_factor]]

    data[[death_ind]][post_end] <- data[[disease_ind]][post_end] <- NA
    data[[death_ind]][pre_end] <- data[[disease_ind]][pre_end] <- FALSE
    data[[death_ind]][is_end] <- data[[death_ind]][is_end] &
      data[[death_time]][is_end] <= data[[disease_time]][is_end]
    data[[disease_ind]][is_end] <- data[[disease_ind]][is_end] &
      data[[disease_time]][is_end] <= data[[death_time]][is_end]

    # adjust times
    data[[ft_end]][post_end] <- NA
    data[[ft_end]][pre_end] <- diff(time_breaks)[data[pre_end, period_factor]]
    data[[ft_end]][is_end] <- data[is_end, ft_end] -
      time_breaks[ft_end_int[is_end]]
  }


  # set period to factor
  data[[period_factor]] <- factor(data[[period_factor]],
                                  labels = levels(cut(numeric(0), time_breaks)))


  # add data to list
  retlist <- c(retlist, list(ID = data[[id_var]],
                             PERIOD = data[[period_factor]]))

  # include only relevant data now
  if (missing(death_time) || missing(death_ind) ||
      missing(disease_time) || missing(disease_ind))
    rel_data <- data[,c(variables), drop = FALSE]
  else
    rel_data <- data[,c(death_ind, disease_ind, ft_end, variables)]

  if (include_period)
    rel_data <- cbind(rel_data, data[,period_factor, drop = FALSE])

  retlist <- c(retlist, list(data = rel_data))


  if (missing(death_time) || missing(death_ind) ||
      missing(disease_time) || missing(disease_ind))
    class(retlist) <- "dpaf_data"
  else
    class(retlist) <- c("dpaf_data", "dpaf_response")

  retlist
}

#' Fit survival regressions for disease PAF calculations
#'
#' @param disease_resp formula with response for disease regression (note RHS
#'   should be \code{.})
#' @param death_resp formula with response for mortality regression (note RHS
#'   should be \code{.})
#' @param predictors formula with RHS containing predictors for both survival
#'   regressions (LHS should be empty)
#' @param dpaf_data an object of class \code{dpaf_response}, from
#'   \code{\link{dpaf_data}}
#' @param ... other arguments to be passed to \code{\link[survival]{survreg}}.
#'
#' @return an object of class \code{dpaf}, containing all the elements of the
#'   \code{\link{dpaf_data}}, plus the following elements:
#'
#'   \item{call}{the function call to fit the models}
#'
#'   \item{survreg_d, survreg_m}{the survival regression objects}
#'
#'   \item{coefficients}{a named matrix of coefficients from the survival
#'   regressions}
#'
#'   \item{var}{a named list of the variance-covariance matrices from the
#'   survival regressions}
#'
#'   Note that the order of names should always be \code{"disease", "mortality"}
#'   for downpipe functions to work.
#' @export
#'
dpaf <- function(disease_resp, death_resp, predictors, dpaf_data, ...) {
  model <- dpaf_data
  stopifnot(inherits(model, "dpaf_response"))

  formula_d <- update(predictors, disease_resp)
  formula_m <- update(predictors, death_resp)

  model$call <- match.call()

  model$survreg_d <- survival::survreg(formula_d, model$data,
                                       dist = "exponential", ...)
  model$survreg_m <- survival::survreg(formula_m, model$data,
                                       dist = "exponential", ...)

  model$coefficients <- do.call(cbind, list(
    disease = model$survreg_d$coefficients,
    mortality = model$survreg_m$coefficients
  ))

  model$var <- list(
    disease = model$survreg_d$var,
    mortality = model$survreg_m$var
  )

  class(model) <- "dpaf"
  model
}

#' Split disease PAF data into groups by a factor
#'
#' This function is a wrapper around R-base's \code{split} function designed for
#' \code{dpaf_data} structures.
#'
#' @param object An object of class \code{\link{dpaf_data}}.
#' @param INDICES A factor or list of factors that can be passed to
#'   \code{\link{split}}. Can optionally be named if a list.
#' @inheritParams name_by
#'
#' @return A named list of \code{\link{dpaf_data}} objects, with names
#'   corresponding to the names of the list output by \code{split}. If a named
#'   list was given in \code{INDICES}, then these names will be included.
#' @export
dpaf_split <- function(object, INDICES, lsep = ": ", fsep = ", ") {
  split_call <- match.call()

  if (!is.null(names(INDICES))) {
    INDICES <- mapply(function(fctr, nm) {
      levels(fctr) <- paste0(nm, lsep, levels(fctr))
      fctr
    }, INDICES, names(INDICES), SIMPLIFY = FALSE)
  }
  IDs <- split(object$ID, INDICES, sep = fsep)
  PERIODs <- split(object$PERIOD, INDICES)
  data <- split(object$data, INDICES)

  mapply(function(ids, periods, dfs) {
    ret <- list(data_call = object$data_call,
                split_call = split_call,
                na.action = object$na.action,
                breaks = object$breaks,
                ID = ids,
                PERIOD = periods,
                data = dfs)
    class(ret) <- c(class(object), "dpaf_split")
    ret
  }, IDs, PERIODs, data, SIMPLIFY = FALSE)
}

# End-user calculations ---------------------------------------------------

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


# Print and summary functions ---------------------------------------------

#' @export
print.dpaf <- function(x, ..., digits = max(getOption("digits") - 3, 3),
                       vsep = strrep('-', getOption("width") * 0.75)) {
  if (!is.null(cl <- x$call)) {
    cat("\nCall:\n")
    dput(cl)
  }
  if (!is.null(dcl <- x$data_call)) {
    cat("\nData preparation call:\n")
    dput(dcl)
  }
  cat(vsep)
  if (!is.null(srd <- x$survreg_d) && !is.null(srm <- x$survreg_m)) {
    cat("\nFormulas:")
    cat("\n  Disease:\n"); dput(stats::formula(srd))

    cat("\n  Mortality:\n"); dput(stats::formula(srm))
    cat("\n")
    cat(vsep)
  }

  if (length(x$na.action)) {
    cat("\nCohort size: ", x$nobs, " (", stats::naprint(x$na.action), ")\n",
        sep = "")
    cat(vsep, "\n")
  } else {
    cat("\nCohort size: ", x$nobs, "\n", sep = "")
    cat(vsep, "\n")
  }

  if (length(cf <- x$coefficients)) {
    cat("Coefficients:\n")
    print(cf, digits = digits, ...)
    cat(vsep, '\n')
  }
  invisible(x)
}

#' Calculate and summarise disease PAFs
#'
#' @param object an object of class \code{dpaf} to summarise
#' @param modlist an optional list of modifications for PAF calculations
#' @param group an optional character vector of column names to group by
#' @param newdata an optional \code{\link{dpaf_data}} object for external
#'   prevalences
#' @param comprehensive shorthand for \code{confint}, \code{group_diffs} and
#'   setting \code{survreg_summ} to \code{'both'} (if \code{TRUE}) or
#'   \code{'none'} (if \code{FALSE})
#' @param confint if \code{TRUE}, calculate confidence intervals for statistics
#' @param group_diffs if \code{TRUE} (and \code{group} is supplied), analyse the
#'   differences between groupwise PAFs
#' @param survreg_summ which survival regression summaries, if any, to include
#' @param ... optional arguments passed to \code{\link{predict.dpaf}},
#'   \code{\link[survival]{summary.survreg}}, \code{\link{dpaf_split}} and other
#'   functions
#'
#' @return a list of potentially interesting quantities of class
#'   \code{'summary.dpaf'}:
#'
#'   \item{call, data_call, nobs, na.action, breaks}{as for \code{dpaf_data} and
#'   \code{dpaf}}
#'
#'   \item{summary_d, summary_m}{summaries of the survival regressions for
#'   disease and mortality respectively}
#'
#'   \item{hazard_ratios}{a list of hazard ratios, possibly with confidence
#'   intervals and standard errors in log-space}
#'
#'   \item{level}{(optional) the level of the confidence interval, if used in
#'   calculations}
#'
#'   \item{modlist}{(optional) the list of modifications used for PAF
#'   calculation}
#'
#'   \item{paf}{(optional) a pointwise disease PAF estimate, possibly with
#'   confidence interval and standard error for the iPAF}
#'
#'   \item{group_pafs}{(optional) pointwise estimates for PAFs in each group, as
#'   well as other information requested}
#'
#'   \item{paf_diffs}{(optional) differences between PAF estimates, with their
#'   standard errors, Z-test statistics, and p-values}
#'
#' @export
summary.dpaf <- function(object, modlist, group, newdata,
                         comprehensive = FALSE,
                         confint = comprehensive, level = 0.95,
                         group_diffs = comprehensive,
                         survreg_summ = if (comprehensive)
                           c('both', 'disease', 'mortality', 'none')
                         else c('none', 'disease', 'mortality', 'both'),
                         ...) {
  x <- object[match(c("call", "data_call", "nobs", "na.action", "breaks"),
                    names(object), nomatch = 0)]

  survreg_summ <- match.arg(survreg_summ)
  if (survreg_summ %in% c('disease', 'both'))
    x$summary_d <- summary(object$survreg_d, ...)
  if (survreg_summ %in% c('mortality', 'both'))
    x$summary_m <- summary(object$survreg_m, ...)

  x$hazard_ratios <- predict.dpaf(object, type = 'hr', confint = confint,
                                  level = level, ...)
  if (confint) x$level <- level

  if (!missing(modlist)) {
    x$modlist <- modlist
    # calculate PAFs
    x$paf <- predict.dpaf(object, modlist, newdata, confint = confint,
                          level = level, gradient = FALSE, ...)

    if (!missing(group)) {
      dpdta <- if (missing(newdata)) object else newdata
      dp_spl <- dpaf_split(dpdta, dpdta$data[, group, drop = FALSE], ...)
      x <- c(x, dpaf_groups(object, modlist, dp_spl, gdiffs = group_diffs, ...))
    }
  }

  class(x) <- 'summary.dpaf'
  return(x)
}

#' @export
#' @method print summary.dpaf
print.summary.dpaf <- function(x, ..., coef_regex = "",
                               digits = max(3L, getOption("digits") - 3L),
                               vsep = strrep("-", getOption("width") * 0.75)) {
  if (!is.null(cl <- x$call)) {
    cat("\nCall:\n")
    dput(cl)
  }
  if (!is.null(dcl <- x$data_call)) {
    cat("\nData preparation call:\n")
    dput(dcl)
  }
  cat(vsep)

  if (!is.null(x$summary_d)) {
    cat("\nDisease regression summary:\n")
    print(x$summary_d, digits = digits, ...)
    cat(vsep)
  }
  if (!is.null(x$summary_m)) {
    cat("\nMortality regression summary:\n")
    print(x$summary_m, digits = digits, ...)
    cat(vsep)
  }


  if (length(x$na.action)) {
    cat("\nCohort size: ", x$nobs, " (", stats::naprint(x$na.action), ")\n",
        sep = "")
    cat(vsep)
  } else {
    cat("\nCohort size: ", x$nobs, "\n", sep = "")
    cat(vsep)
  }

  if (!is.null(x$hazard_ratios)) {
    cat("\nHazard ratios:\n")
    coef_rows <- grep(coef_regex, rownames(x$hazard_ratios$disease))

    cat("\n    Disease:\n")
    hr_d <- x$hazard_ratios$disease[coef_rows,]
    if (!is.null(x$level)) {
      colnames(hr_d)[which(colnames(hr_d) == "lwr")] <-
        paste("Lower", paste0(format(100 * x$level, digits = digits, ...), "%"),
              "CL")
      colnames(hr_d)[which(colnames(hr_d) == "upr")] <-
        paste("Upper", paste0(format(100 * x$level, digits = digits, ...), "%"),
              "CL")
    }
    print(hr_d, digits = digits, ...)

    cat("\n    Mortality:\n")
    hr_d <- x$hazard_ratios$mortality[coef_rows,]
    if (!is.null(x$level)) {
      colnames(hr_d)[which(colnames(hr_d) == "lwr")] <-
        paste("Lower", paste0(format(100 * x$level, digits = digits, ...), "%"),
              "CL")
      colnames(hr_d)[which(colnames(hr_d) == "upr")] <-
        paste("Upper", paste0(format(100 * x$level, digits = digits, ...), "%"),
              "CL")
    }
    print(hr_d, digits = digits, ...)

    cat(vsep)
  }

  if (!is.null(x$modlist)) {
    cat("\nModifications:\n")
    lapply(X = seq_along(x$modlist), FUN = function(i, lst) {
      if (length(lst[[i]]) > 1)
        cat(names(lst)[i], ": ",  paste(lst[[i]][-1], collapse = ", "),
            " -> ", lst[[i]][1], "\n", sep = "")
      else
        cat(names(lst)[i], ": ", ".",
            " -> ", lst[[i]][1], "\n", sep = "")
    }, lst = x$modlist)
    cat(vsep)
  }

  if (!is.null(x$paf)) {
    cat("\nPAF for disease incidence, censored by death:\n")
    if (!is.null(x$level)) {
      names(x$paf)[which(names(x$paf) == "lwr")] <-
        paste("Lower", paste0(format(100 * x$level, digits = digits, ...), "%"),
              "CL")
      names(x$paf)[which(names(x$paf) == "upr")] <-
        paste("Upper", paste0(format(100 * x$level, digits = digits, ...), "%"),
              "CL")
    }
    print(x$paf, digits = digits, ...)
    cat(vsep)
  }

  if (!is.null(pafs <- x$group_pafs)) {
    if (!is.null(x$level)) {
      colnames(pafs)[which(colnames(pafs) == "lwr")] <-
        paste("Lower", paste0(format(100 * x$level, digits = digits, ...), "%"),
              "CL")
      colnames(pafs)[which(colnames(pafs) == "upr")] <-
        paste("Upper", paste0(format(100 * x$level, digits = digits, ...), "%"),
              "CL")
    }

    cat("\nGroupwise PAF estimates:\n")
    print(pafs, digits = digits, ...)
    cat(vsep)
  }

  if (!is.null(x$paf_diffs)) {
    cat("\nAnalysis of differences between groupwise PAFs:\n")
    printCoefmat(x$paf_diffs, digits = digits, k = 2)
    cat(vsep)
  }
}

# Intermediate calculation functions --------------------------------------

#' Calculate hazards for dpaf study
#'
#' @param z model matrix
#' @param cf named matrix of coefficients, as in \code{\link{dpaf}}
#'
#' @return named matrix of hazards
#' @keywords internal
dpaf_hz <- function(z, cf) {
  exp(z %*% -cf)
}

#' Calculate survivals for dpaf study
#'
#' @param hz named matrix of hazards
#' @param ID vector of IDs
#' @param PERIOD vector of periods
#' @param dt delta-times -- length of each period
#'
#' @return named matrix of survivals
#' @keywords internal
dpaf_sv <- function(hz, ID, PERIOD, dt) {
  if (any(tapply(PERIOD, ID, is.unsorted)))
    stop("Periods must be in ascending order for each ID to calculate survival")

  apply(hz, 2, function(hz_col)
    exp(-stats::ave(hz_col, ID, FUN = function(l) cumsum(dt * l)))
  )
}

#' Calculate disease-free survival (**s**ur**v**ival **p**roducts)
#'
#' @param sv named matrix of survivals
#'
#' @return vector of disease-free survivals
#' @keywords internal
dpaf_svp <- function(sv) {
  apply(sv, 1, prod)
}

#' Calculate healthy mortality \eqn{I}
#'
#' @param hz **named** matrix of hazards (must have a column named "disease")
#' @param svp vector of disease-free survivals
#' @param ID vector of IDs
#' @param PERIOD factor vector of PERIODs
#'
#' @return (numeric, scalar) healthy mortality, perhaps to some constant factor
#' @keywords internal

dpaf_i <- function(hz, svp, ID, PERIOD) {
  if (any(tapply(PERIOD, ID, is.unsorted)))
    stop("Periods must be in ascending order for each ID to calculate I")

  svd <- stats::ave(svp, ID, FUN = function(s) -diff(c(1, s)))
  ## (survival at time 0 is 1)
  sum(hz[,"disease"] / rowSums(hz) * svd)
}

#' Calculate gradient of hazard in disease PAF studies
#'
#' @param z model matrix
#' @param hz named matrix of hazards
#'
#' @return a named list of matrices containing hazard gradients
#' @keywords internal
dpaf_ghz <- function(z, hz) {
  lapply(split(hz, col(hz, as.factor = TRUE)), `*`, z)
}

#' Calculate gradient of survival in disease PAF studies
#'
#' @param ghz a named list of hazard gradients
#' @inheritParams dpaf_sv
#'
#' @return a named list of matrices containing survival gradients
#' @keywords internal
dpaf_gsv <- function(ghz, sv, ID, PERIOD, dt) {
  if (any(tapply(PERIOD, ID, is.unsorted)))
    stop("Periods must be in ascending order for each ID")

  sv <- split(sv, col(sv, as.factor = TRUE))
  mapply(
    function(ghz_x, sv_x) -sv_x * apply(
      ghz_x, 2, function(col) stats::ave(
        col, ID, FUN = function(ghz_ij) cumsum(ghz_ij * dt)
      )
    ),
    ghz, sv, SIMPLIFY = FALSE
  )
}

#' Calculate gradient of healthy mortality \eqn{I}
#'
#' @param ghz a named list of hazard gradients
#' @param gsv a named list of survival gradients
#' @param sv a named matrix of survivals
#' @inheritParams dpaf_i
#'
#' @return a named list of vectors for \eqn{\grad{I}} for each of disease and
#'   mortality
#' @keywords internal
dpaf_gi <- function(ghz, gsv, hz, sv, ID, PERIOD) {
  if (any(tapply(PERIOD, ID, is.unsorted)))
    stop("Periods must be in ascending order for each ID")

  svd <- stats::ave(apply(sv, 1, prod), ID, FUN = function(s) -diff(c(1, s)))
  gsvd_d <- apply(
    apply(sv, 1, prod) * gsv$disease, 2,
    # sv[,"mortality"] * gsv$disease, 2,
    function(gsvp_col) stats::ave(gsvp_col, ID, FUN = function(s) -diff(c(0,s)))
  )
  gsvd_m <- apply(
    apply(sv, 1, prod) * gsv$mortality, 2,
    # sv[,"disease"] * gsv$mortality, 2,
    function(gsvp_col) stats::ave(gsvp_col, ID, FUN = function(s) -diff(c(0,s)))
  )

  smnd_d <- hz[,"mortality"] / rowSums(hz)**2 * svd * ghz$disease +
    hz[,"disease"] / rowSums(hz) * gsvd_d
  smnd_m <- -hz[,"disease"] / rowSums(hz)**2 * svd * ghz$mortality +
    hz[,"disease"] / rowSums(hz) * gsvd_m

  lapply(
    list(disease = smnd_d, mortality = smnd_m),
    function(smnd) apply(smnd, 2, sum)
  )
}

#' Calculate gradients of \eqn{\log(1-PAF)} wrt disease and mortality
#' coefficients
#'
#' @param gi_mod named list of gradients of \eqn{I^*}
#' @param i_mod \eqn{I^*}
#' @param gi_raw named list of gradients of \eqn{I}
#' @param i_raw \eqn{I}
#'
#' @return a named list of gradient vectors of the iPAF
#' @keywords internal
dpaf_gipaf <- function(gi_mod, i_mod, gi_raw, i_raw) {
  mapply(
    `-`,
    lapply(gi_mod, `/`, i_mod),
    lapply(gi_raw, `/`, i_raw),
    SIMPLIFY = FALSE
  )
}

#' Calculate gradients of PAF wrt disease and mortality coefficients
#'
#' @param gi_mod named list of gradients of \eqn{I^*}
#' @param i_mod \eqn{I^*}
#' @param gi_raw named list of gradients of \eqn{I}
#' @param i_raw \eqn{I}
#'
#' @return a named list of gradient vectors of the PAF
#' @keywords internal
dpaf_gpaf <- function(gi_mod, i_mod, gi_raw, i_raw) {
  mapply(
    function(gi_mod_x, gi_raw_x)
      (gi_raw_x * i_mod - gi_mod_x * i_raw) / i_raw**2,
    gi_mod, gi_raw,
    SIMPLIFY = FALSE
  )
}
