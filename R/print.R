######### Functions to print PAF objects ##########

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
    stats::printCoefmat(x$paf_diffs, digits = digits, k = 2)
    cat(vsep)
  }
}

#' @export
print.mpaf_summary <- function(x, ..., coef_regex = "",
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

  if (!is.null(x$survreg)) {
    cat("\nSurvival regression summary:\n")
    print(summary(x$survreg), digits = digits, ...)
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

  if (!is.null(x$HR)) {
    cat("\nHazard ratios:\n")
    print(x$HR, digits = digits, ...)

    cat(vsep)
  }

  if (!is.null(x$modifications)) {
    cat("\nModifications:\n")
    lapply(X = seq_along(x$modifications), FUN = function(i, lst) {
      if (length(lst[[i]]) > 1)
        cat(names(lst)[i], ": ",  paste(lst[[i]][-1], collapse = ", "),
            " -> ", lst[[i]][1], "\n", sep = "")
      else
        cat(names(lst)[i], ": ", ".",
            " -> ", lst[[i]][1], "\n", sep = "")
    }, lst = x$modifications)
    cat(vsep)
  }

  if (!is.null(x$paf)) {
    cat("\nPAFs for death:\n")
    print(rbind(x$paf, x$paf0), digits = digits, ...)
    cat(vsep)
  }

  if (!is.null(x$group_pafs)) {
    cat("\nGroupwise PAF estimates:\n")
    lapply(seq_along(x$group_pafs), function(i, gpaf, nm) {
      cat(paste0("\n    ", nm[i], "\n    ", strrep('-', nchar(nm[i])), "\n"))
      print(rbind(gpaf[[i]]$paf, gpaf[[i]]$paf0), digits = digits, ...)
    }, gpaf = x$group_pafs, nm = names(x$group_pafs))

    cat(vsep)
  }

  if (!is.null(x$paf_diffs)) {
    cat("\nAnalysis of differences between groupwise PAFs:\n")
    lapply(seq_along(x$paf_diffs), function(i, dpaf, nm) {
      cat(paste0("\n    ", nm[i], "\n    ", strrep('-', nchar(nm[i])), "\n"))
      stats::printCoefmat(dpaf[[i]], digits = digits, k = 2)
    }, dpaf = x$paf_diffs, nm = names(x$paf_diffs))
    cat(vsep)
  }
}
