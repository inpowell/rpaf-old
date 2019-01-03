
#' Prepare data for mortality PAF calculation
#'
#' Fit a parametric survival model in a suitable manner for subsequent
#' calculations for mortality PAFs.
#'
#' Each row of \code{rawdata} should contain information about an individual
#' (with unique ID given in \code{id_var}, birth cohort in \code{cohort_var}),
#' including if they died before the end of follow-up (in \code{censor_var}) and
#' when they were followed up until (\code{censor_time}). Further variables of
#' interest should be included in \code{covariate_model} to be included in the
#' survival regression, with interactions given by column names separated by a
#' colon.
#'
#' @param rawdata a data frame containing raw data from cohort study
#' @param id_var name of column in \code{rawdata} for individual ID (integer,
#'   unique)
#' @param censor_var name of column in \code{rawdata} which is \code{TRUE} if
#'   censoring occurs before follow-up (logical)
#' @param censor_time name of column in \code{rawdata} for when censoring occurs
#'   (numeric)
#' @param cohort_var name of column in \code{rawdata} which gives birth cohort
#'   (factor)
#' @param covariate_model character string of covariate model of interest (to be
#'   passed to \code{reformulate})
#' @param ft_length numeric length of follow-up time
#' @param ft_delta numeric length of constant hazard intervals
#' @param period_var name of new column in output identifying the current
#'   follow-up period
#' @param na.action a missing-data filter function applied to the model frame
#' @param survreg.control control parameters passed to the survival regression
#'   function (see \code{\link[survival]{survreg.control}})
#'
#' @param period_time_var name of new column representing censor time in each
#'   follow-up period (note this is \code{NA} in any period after censoring
#'   occurs)
#' @param period_censor_var name of new column which indicates if censoring has
#'   occurred
#'
#' @importFrom stats na.fail reformulate terms update
#'
#' @return an object of class \code{mpafreg}
#' @export
#'
#' @family mortality cohort functions
mpafreg <- function(
  rawdata,
  id_var,
  censor_var,
  censor_time,
  cohort_var,
  covariate_model,
  ft_length,
  ft_delta,
  period_var,
  period_time_var = "f_time",
  period_censor_var = "censor",
  na.action = getOption("na.action"),
  survreg.control = survival::survreg.control()
) {
  # # TODO
  # warning("You ran `prepare_mortality`, which hasn't been tested yet")
  covars <- as.character(
    attr(terms(reformulate(covariate_model)), "variables"))[-1]

  # handle missing data off the bat
  data <- match.fun(na.action)(
    subset(rawdata, select = c(id_var, censor_var, censor_time,
                               cohort_var, covars))
  )

  n_fperiod <- ft_length %/% ft_delta
  period_df <- data.frame((0:n_fperiod) * ft_delta)
  names(period_df) <- period_var

  # block by time periods
  surv_data <- merge(period_df, data, by = NULL)

  # calculate periodic follow-up time
  f_time <- surv_data[,censor_time] - surv_data[,period_var]

  # remove non-positive follow-up times
  f_time[f_time <= 0] <- NA

  # limit follow-up times to be the follow-up time interval
  f_time[f_time >= ft_delta] <- ft_delta

  # build the censor variable
  censor <- surv_data[,censor_var] & f_time < ft_delta

  # attach columns to data frame
  surv_data[,period_time_var] <- f_time
  surv_data[,period_censor_var] <- censor

  # make the time period a factor
  surv_data[,period_var] <- as.factor(surv_data[,period_var])

  # predictors
  pred_formula <- reformulate(c(paste(cohort_var, period_var, sep = "*"),
                                covariate_model), intercept = FALSE)

  full_formula <- update(
    pred_formula,
    survival::Surv(get(period_time_var), get(period_censor_var)) ~ .
  )


  # perform the regression
  censor_reg <- survival::survreg(
    full_formula,
    stats::na.exclude(surv_data), # only removes observations post-censor
    dist = "exponential",
    na.action = na.fail, # there should be no NA's
    score = TRUE, y = FALSE,
    control = survreg.control
  )

  retlist <- list(
    call = match.call(),
    nobs = nrow(data),
    survreg = censor_reg,
    model = surv_data,
    id_var = id_var,
    period_var = period_var,
    period_time_var = period_time_var,
    period_censor_var = period_censor_var,
    censor_var = censor_var,
    censor_time = censor_time,
    cohort_var = cohort_var,
    covariate_model = covariate_model,
    ft_length = ft_length,
    ft_delta = ft_delta,
    na.action = na.action(data)
  )
  class(retlist) <- "mpafreg"

  retlist
}

#' Test for a significant difference between two PAFs
#'
#' @param pafs1,pafs2 numerical vectors of PAF estimates
#' @param gpafs1,gpafs2 numerical matrices of PAF gradients (row corresponds to
#'   time period, column to coefficient)
#' @param vcov variance-covariance matrix of \code{survreg} object
#'
#' @return a vector of four elements: difference, std err of difference, Z-test
#'   value, and P-value
#' @export
#'
compare_mpafs <- function(pafs1, gpafs1, pafs2, gpafs2, vcov) {
  dpafs <- pafs1 - pafs2
  g <- gpafs1 - gpafs2
  vv <- g %*% vcov %*% t(g)
  se <- sqrt(diag(vv))
  z <- dpafs / se
  p <- 2 * stats::pnorm(-abs(z))

  cbind("PAF Diff" = dpafs, "SE(PAF Diff)" = se, "Z value" = z, "Pr(>|Z|)" = p)
}

#' Convert a data frame to a type compatible with an \code{mpafreg} object.
#'
#' @param object An object of class \code{mpafreg}
#' @param rawdata A data frame to be expanded and coerced, similar to the
#'   \code{rawdata} argument in \code{\link{mpafreg}}.
#'
#' @return a data frame ready to be used with \code{\link{predict.mpafreg}}.
#' @export
#'
make_data.mpafreg <- function(object, rawdata) {

  covars <- as.character(
    attr(terms(reformulate(object$covariate_model)), "variables"))[-1]

  # no missing data treatment
  data <- subset(rawdata, select =
                   with(object, c(id_var, censor_var, censor_time,
                                     cohort_var, covars)))


  n_fperiod <- with(object, ft_length %/% ft_delta)
  period_df <- data.frame((0:n_fperiod) * object$ft_delta)
  names(period_df) <- object$period_var

  # block by time periods
  newdata <- merge(period_df, data, by = NULL)

  # calculate periodic follow-up time
  f_time <- with(object, newdata[,censor_time] - newdata[,period_var])

  # remove non-positive follow-up times
  f_time[f_time <= 0] <- NA

  # limit follow-up times to be the follow-up time interval
  f_time[f_time >= object$ft_delta] <- object$ft_delta

  # build the censor variable
  censor <- newdata[,object$censor_var] & f_time < object$ft_delta

  # attach columns to data frame
  newdata[,object$period_time_var] <- f_time
  newdata[,object$period_censor_var] <- censor

  # make the time period a factor
  newdata[,object$period_var] <- as.factor(newdata[,object$period_var])

  newdata
}

#' @export
print.mpafreg <- function(x, ...) {
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }

  sr <- x$survreg
  sr$call <- NULL

  form <- stats::formula(sr)
  attributes(form) <- NULL

  cat("\nSurvreg formula:\n")
  dput(form)
  cat("\n")

  cat(strrep("-", 0.75 * getOption("width")), "\n")
  print(sr, ...)
  cat(strrep("-", 0.75 * getOption("width")), "\n")

  if (length(x$na.action)) {
    cat("\nCohort size: ", x$nobs, " (", stats::naprint(x$na.action), ")\n",
        sep = "")
  } else {
    cat("\nCohort size: ", x$nobs, "\n", sep = "")
  }

  invisible(x)
}

#' Calculate and summarise mortality PAFs
#'
#' @param object an \code{mpafreg} object to summarise
#' @param modlist an optional list of modifications to make (see
#'   \code{\link{apply_modifications}} for structure)
#' @param group an optional character vector of columns to group PAF
#'   calculations by
#' @param newdata an optional data frame of the form of \code{object$model} to
#'   calculate PAFs using
#' @param comprehensive include detailed information at the cost of computing
#'   time? Sets the default value for \code{confint}, \code{vcov.trans},
#'   \code{group_diffs} and \code{survreg.summ}
#' @param confint calculate confidence intervals?
#' @param vcov.trans calculate variance-covariance matrices in normal space?
#' @param group_diffs test significance of differences between groups?
#' @param survreg.summ include summary of the underlying \code{survreg} object?
#' @param ... extra parameters to be passed to \code{predict.mpafreg}
#'
#' @return a list of summary statistics which might even be interesting
#' @export
#'
summary.mpafreg <- function(object, modlist, group, newdata,
                            comprehensive = FALSE,
                            confint = comprehensive,
                            vcov.trans = comprehensive,
                            group_diffs = comprehensive,
                            survreg.summ = comprehensive,
                            ...) {
  x <- object[match(c("call", "nobs", "na.action"),
                    names(object), nomatch = 0)]
  x <- c(x, list(confint = confint, covars = as.character(attr(terms(
    reformulate(object$covariate_model)), "variables"))[-1]))

  if (survreg.summ)
    x <- c(x, list(survreg.summ = summary(object$survreg)))

  if (missing(newdata))
    df <- object$model
  else
    df <- newdata

  # hazard ratios -----------------------------------------------------------

  hazard_ratios <- stats::predict(object, type = 'hr', confint = confint,
                                  vcov.trans = vcov.trans, ...)


  if (vcov.trans) {
    x <- c(x, list(hazard_ratios = hazard_ratios[["prediction"]],
                   var_hr = hazard_ratios[["vcov.trans"]]))
  } else if (is.list(hazard_ratios)) {
    x <- c(x, list(hazard_ratios = hazard_ratios[["prediction"]]))
  } else {
    x <- c(x, list(hazard_ratios = hazard_ratios))
  }

  # no modifications => no PAFs
  if (missing(modlist)) {
    class(x) <- "summary.mpafreg"
    return(x)
  }

  # PAFs --------------------------------------------------------------------

  x <- c(x, list(modlist = modlist))
  total_pafs <- stats::predict(object, modlist, df, type = 'paf',
                               confint = confint,
                               vcov.trans = vcov.trans, ...)

  if (vcov.trans) {
    x <- c(x, list(pafs = total_pafs[["prediction"]],
                   var_paf = total_pafs[["vcov.trans"]]))
  } else if (is.list(total_pafs)) {
    x <- c(x, list(pafs = total_pafs[["prediction"]]))
  } else {
    x <- c(x, list(pafs = total_pafs))
  }

  # no group variable => no groups
  if (missing(group)) {
    class(x) <- "summary.mpafreg"
    return(x)
  }

  # Groups ------------------------------------------------------------------

  if (!group_diffs) {
    group_paf_lists <- name_by(by(
      df, df[,group, drop = FALSE],
      function(group_df)
        stats::predict(object, modlist, group_df, type = 'paf',
                       confint = confint, vcov.trans = vcov.trans,
                       ...)
    ))

    if (vcov.trans) {
      group_pafs <- lapply(group_paf_lists, `[[`, "prediction")
      var_group_pafs <- lapply(group_paf_lists, `[[`, "vcov.trans")

      x <- c(x, list(group_pafs = group_pafs,
                     var_group_pafs = var_group_pafs))
    } else if (is.list(group_paf_lists[[1]])) {
      group_pafs <- lapply(group_paf_lists, `[[`, "prediction")
      x <- c(x, list(group_pafs = group_pafs))
    } else {
      x <- c(x, list(group_pafs = group_paf_lists))
    }
  } else {
    # Differences between groups ----------------------------------------------
    group_paf_lists <- name_by(by(
      df, df[,group, drop = FALSE],
      function(group_df)
        stats::predict(object, modlist, group_df, type = 'gpaf',
                       confint = confint, vcov.trans = vcov.trans,
                       ...)
    ))

    if (vcov.trans) {
      group_pafs <- lapply(group_paf_lists, `[[`, "prediction")
      var_group_pafs <- lapply(group_paf_lists, `[[`, "vcov.trans")

      x <- c(x, list(group_pafs = group_pafs,
                     var_group_pafs = var_group_pafs))
    } else { # will always have gradient element
      group_pafs <- lapply(group_paf_lists, `[[`, "prediction")
      x <- c(x, list(group_pafs = group_pafs))
    }

    group_pairs <- utils::combn(seq_along(group_paf_lists), 2)
    group_pairs <- split(group_pairs, col(group_pairs))

    group_comps <- lapply(X = group_pairs, FUN = function(ij, x, vcov) {
      x1 <- x[[ij[1]]]
      x2 <- x[[ij[2]]]

      pafs1 <- x1[["prediction"]][,"fit"]
      gpafs1 <- x1[["gradient"]]
      pafs2 <- x2[["prediction"]][,"fit"]
      gpafs2 <- x2[["gradient"]]

      compare_mpafs(pafs1, gpafs1, pafs2, gpafs2, vcov)
    }, x = group_paf_lists, vcov = object$survreg$var)

    names(group_comps) <- lapply(X = group_pairs, FUN = function(ij, nam) {
      paste(nam[ij], collapse = " - ")
    }, nam = names(group_paf_lists))

    x <- c(x, list(group_comparison = group_comps))
  }

  class(x) <- "summary.mpafreg"
  return(x)
}

#' @export
#' @method print summary.mpafreg
print.summary.mpafreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  vsep = strrep("-", getOption("width") * 0.75),
                                  ...) {
  cat("\nCall:\n")
  dput(x$call)

  if (!is.null(x$survreg.summ)) {
    cat(vsep, "\nSurvival regression:\n")
    print(x$survreg.summ)
  }

  if (!is.null(x$na.action)) {
    cat(vsep)
    cat("\nCohort size: ", x$nobs, " (", stats::naprint(x$na.action), ")\n",
        sep = "")
  } else {
    cat("\nCohort size: ", x$nobs, "\n", sep = "")
  }

  if (!is.null(x$hazard_ratios)) {
    cat(vsep)
    cat("\nHazard ratios:\n")

    hr_regex <- paste0("^", x$covars, collapse = "|")
    if (x$confint) {
      hr_table <- x$hazard_ratios[grep(hr_regex, rownames(x$hazard_ratios)),]
      clevel <- attr(x$hazard_ratios, "level")
      attr(x$hazard_ratios, "level") <- NULL
      colnames(hr_table) <- c(
        "HR",
        paste("Lower", paste0(format(100 * clevel, digits = digits, ...), "%"),
              "CL"),
        paste("Upper", paste0(format(100 * clevel, digits = digits, ...), "%"),
              "CL")
      )
      print(hr_table, digits = digits, ...)
    } else {
      print(x$hazard_ratios[grep(hr_regex, names(x$hazard_ratios))],
            digits = digits, ...)
    }
  }

  if (!is.null(x$modlist)) {
    cat(vsep)
    cat("\nModifications:\n")
    lapply(X = seq_along(x$modlist), FUN = function(i, lst) {
      if (length(lst[[i]]) > 1)
        cat(names(lst)[i], ": ",  paste(lst[[i]][-1], collapse = ", "),
            " -> ", lst[[i]][1], "\n", sep = "")
      else
        cat(names(lst)[i], ": ", ".",
            " -> ", lst[[i]][1], "\n", sep = "")
    }, lst = x$modlist)
  }

  if (!is.null(x$pafs)) {
    cat(vsep)
    cat("\nPAFs:\n")

    cat("\nWhole cohort\n------------\n")

    if (x$confint) {
      paf_table <- x$pafs
      clevel <- attr(x$pafs, "level")
      attr(paf_table, "level") <- NULL
      colnames(paf_table) <- c(
        "PAF",
        paste("Lower", paste0(format(100 * clevel, digits = digits, ...), "%"),
              "CL"),
        paste("Upper", paste0(format(100 * clevel, digits = digits, ...), "%"),
              "CL")
      )
      print(paf_table, digits = digits, ...)
    } else {
      print(x$pafs, digits = digits, ...)
    }
  }

  if (!is.null(x$group_pafs)) {
    cat("\n")
    cat("\n    By group:\n")

    lapply(X = seq_along(x$group_pafs), FUN = function(i, pafs_list) {
      name <- names(pafs_list)[i]
      cat("\n", name, "\n", rep("-", nchar(name)), "\n", sep = "")
      pafs <- pafs_list[[i]]
      if (x$confint) {
        paf_table <- pafs
        clevel <- attr(paf_table, "level")
        attr(paf_table, "level") <- NULL
        colnames(paf_table) <- c(
          "PAF",
          paste("Lower", paste0(format(100 * clevel, digits = digits, ...), "%"),
                "CL"),
          paste("Upper", paste0(format(100 * clevel, digits = digits, ...), "%"),
                "CL")
        )
        print(paf_table, digits = digits, ...)
      } else {
        print(pafs, digits = digits, ...)
      }

    }, pafs_list = x$group_pafs)
  }

  if (!is.null(x$group_comparison)) {
    cat("\n")
    cat(vsep)
    cat("\nPairwise comparisons between groups:\n")

    coefmat <- t(sapply(x$group_comparison, function(table) table[1,]))
    times <- sapply(x$group_comparison, function(table) row.names(table)[1])

    if (length(unique(times)) == 1)
      cat("\nTime interval:", times[1], "\n")
    else
      row.names(coefmat) <- paste(row.names(coefmat), "; t =", times)

    stats::printCoefmat(coefmat, digits = digits, cs.ind = 1:2, tst.ind = 3,
                        P.values = TRUE, has.Pvalue = TRUE)
  }

  invisible(x)
}
