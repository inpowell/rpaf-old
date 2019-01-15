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

#' @export
print.dpaf_summary <- function(x, ..., coef_regex = "",
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

  if (!is.null(x$survreg_d) || !is.null(x$survreg_m)) {
    cat("\nSurvival regression summary:\n")

    if (!is.null(x$survreg_d)) {
      cat("\n    Disease:\n")
      print(summary(x$survreg_d), digits = digits, ...)
    }

    if (!is.null(x$survreg_m)) {
      cat("\n    Mortality:\n")
      print(summary(x$survreg_m), digits = digits, ...)
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

  if (!is.null(x$HR_d) || !is.null(x$HR_m)) {
    cat("\nHazard ratios:\n")

    if (!is.null(x$HR_d)) {
      cat("\n    Disease:\n")
      print(x$HR_d, digits = digits, ...)
    }

    if (!is.null(x$HR_m)) {
      cat("\n    Mortality:\n")
      print(x$HR_m, digits = digits, ...)
    }

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
}
