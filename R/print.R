######### Functions to print PAF objects ##########

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
    print(x$HR, digits = digits, na.print = "--", ...)

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
  cat('\n')
  invisible(x)
}

#' @export
print.dpaf_summary <- function(x, ..., coef_regex = "",
                               digits = max(3L, getOption("digits") - 3L),
                               vsep = strrep("-", getOption("width") * 0.75)) {
  if (!is.null(cl <- x$call)) {
    cat("\nCall:\n")
    dput(cl)
  }
  cat(vsep)

  if (!is.null(x$survreg_1) || !is.null(x$survreg_2)) {
    cat("\nSurvival regression summary:\n")

    if (!is.null(x$survreg_1)) {
      cat("\n    Disease:\n")
      print(summary(x$survreg_1), digits = digits, ...)
    }

    if (!is.null(x$survreg_2)) {
      cat("\n    Mortality:\n")
      print(summary(x$survreg_2), digits = digits, ...)
    }
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

  if (!is.null(x$HR_1) || !is.null(x$HR_2)) {
    cat("\nHazard ratios:\n")

    if (!is.null(x$HR_1)) {
      cat("\n    Disease:\n")
      print(x$HR_1, digits = digits, na.print = '--', ...)
    }

    if (!is.null(x$HR_2)) {
      cat("\n    Mortality:\n")
      print(x$HR_2, digits = digits, na.print = '--', ...)
    }

    cat(vsep)
  }

  if (!is.null(x$modifications_2)) {
    cat("\nModifications:\n")
    lapply(X = seq_along(x$modifications_2), FUN = function(i, lst) {
      if (length(lst[[i]]) > 1)
        cat(names(lst)[i], ": ",  paste(lst[[i]][-1], collapse = ", "),
            " -> ", lst[[i]][1], "\n", sep = "")
      else
        cat(names(lst)[i], ": ", ".",
            " -> ", lst[[i]][1], "\n", sep = "")
    }, lst = x$modifications_2)
    cat(vsep)
  }

  if (!is.null(x$paf)) {
    cat("\nPAFs for disease:\n")
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
  cat('\n')
  invisible(x)
}
