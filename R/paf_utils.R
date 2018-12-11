#' Hazard, survival and mortality statistics
#'
#' \code{survival_statistics} returns a list of statistics important in PAF
#' calculation, namely the hazards, survivals, and mortalities for each period
#' of interest, and optionally their gradients.
#'
#' This function assumes that \code{model == model[order(ID, PERIOD),]} and the
#' levels of \code{PERIOD} are in ascending order. If \code{tapply} does not
#' preserve order of arguments (for \code{cumsum}), then this function is
#' broken.
#'
#' Note that \code{as.numeric(levels(PERIOD)[PERIOD])} should return a numeric
#' vector of times at the start of each period.
#'
#' @section TODO: convert use of \code{cumsum} to \code{ave} for greater
#'   robustness to disorder
#'
#' @param model complete model matrix
#' @param survreg_coef the coefficients from a \code{survreg} object
#' @param ID a numerical vector of individual ID's corresponding to rows of
#'   \code{model}
#' @param PERIOD a factor of time period beginnings corresponding to rows of
#'   \code{model}
#' @inheritParams mpafreg
#' @param gradient if \code{TRUE}, also calculate gradients for the hazards,
#'   survivals and mortalities.
#'
#' @return a list with named elements \code{hazard}, \code{survival},
#'   \code{survival_ft}, and \code{mortality} (where \code{mortality} is named
#'   based on time period), and optionally their gradients which are matrices
#'   with names as previous, but prepended by \code{grad_}. Note that
#'   \code{survival_ft} corresponds to the survival at the end of follow-up as
#'   given by \code{ft_length}, which will be the first entry in calculated
#'   mortalities.
#' @export
#'
survival_statistics <- function(
  model,
  survreg_coef,
  ID,
  PERIOD,
  ft_length,
  ft_delta,
  gradient = TRUE
) {
  # # TODO
  # warning("You ran `survival_statistics`, which is currently unavailable")

  # calculate fitted hazards from model matrix and coefficients
  hazard <- as.vector(exp(model %*% -survreg_coef))
  # calculate survivals at the end of each period
  survival <- unlist(
    tapply(hazard, ID, function(id_hazard) exp(-ft_delta * cumsum(id_hazard)))
  )

  ## calculate differences in survival for each period
  survival_diffs <- unlist(
    tapply(survival, ID, function(id_surv) -diff(c(1, id_surv)))
  ) ## note that survival at time zero is 1

  mortality <- as.vector(
    tapply(survival_diffs, PERIOD, mean)
  )
  names(mortality) <- paste0(
    "(", levels(PERIOD), ",", as.numeric(levels(PERIOD)) + ft_delta, "]"
  )

  total_mortality <- as.vector(
    tapply(1 - survival, PERIOD, mean)
  )
  names(total_mortality) <- paste0(
    "(0,", as.numeric(levels(PERIOD)) + ft_delta, "]"
  )

  hazard_ft_contrib <- pmin(ft_length - as.numeric(levels(PERIOD)), ft_delta)
  survival_ft <- unlist(
    tapply(hazard, ID, function(hazard) exp(-sum(hazard_ft_contrib * hazard)))
  )

  mortality_ft <- mean(1 - survival_ft)
  names(mortality_ft) <- paste0("(0,", ft_length, "]")

  mortality_all <- c(mortality_ft, mortality, total_mortality)

  retlist <- list(
    ### scream test removal
    # hazard = hazard,
    # survival = survival,
    # survival_ft = survival_ft,
    mortality = mortality_all
  )

  # calculate gradients, if requested
  if (gradient) {
    grad_hazard <- hazard * model # sic elementwise multiplication
    grad_sumhazard <- apply(grad_hazard, 2, function(col) unlist(
      tapply(col, ID, cumsum)
    ))
    grad_survival <- apply(grad_sumhazard, 2, function(col)
      -col * ft_delta * survival
    )
    # diffs only used internally -- not returned in list
    grad_survival_diffs <- apply(grad_survival, 2, function(col) unlist(
      tapply(col, ID, function(surv) -diff(c(0, surv)))
    )) ## (note grad(survival) at time zero is 0 as survival is constant here)

    # calculate mortalities for each time period, and cumulatively
    grad_mortality <- apply(grad_survival_diffs, 2, function(col) unlist(
      tapply(col, PERIOD, mean, na.rm = TRUE)
    ))
    grad_total_mortality <- apply(-grad_survival, 2, function(col) unlist(
      tapply(col, PERIOD, mean, na.rm = TRUE)
    ))

    # calculate logarithmic gradients of mortality
    grad_lmortality <- grad_mortality / mortality
    grad_total_lmortality <- grad_total_mortality / total_mortality

    # calculate statistics over time period (0, ftlength]
    grad_survival_ft <- apply(grad_hazard, 2, function(col) unlist(
      tapply(col, ID, function(ghaz) -sum(ghaz * hazard_ft_contrib))
    ) * survival_ft)

    grad_mortality_ft <- apply(-grad_survival_ft, 2, mean, na.rm = TRUE)
    grad_lmortality_ft <- grad_mortality_ft / mortality_ft

    grad_mortality_all <- rbind(grad_mortality_ft, grad_mortality,
                                grad_total_mortality)
    grad_lmortality_all <- rbind(grad_lmortality_ft, grad_lmortality,
                                 grad_total_lmortality)
    row.names(grad_mortality_all) <- names(mortality_all)
    row.names(grad_lmortality_all) <- names(mortality_all)

    retlist <- append(retlist, list(
      ### scream test removal
      # grad_hazard = grad_hazard,
      # grad_survival = grad_survival,
      # grad_survival_ft = grad_survival_ft,
      grad_mortality = grad_mortality_all,
      grad_lmortality = grad_lmortality_all
    ))
  }
  retlist
}


#' Modify data frame for PAF calculations
#'
#' This function returns a data frame which has had the specified modifications
#' applied.
#'
#' \code{modlist} should be a named list of scalars or vectors. Each name should
#' correspond to a column in \code{df}. The first element in each vector should
#' be the desired output value. If there are subsequent elements in the vector,
#' these identify which values to convert from.
#'
#' @param df a data frame to modify
#' @param modlist named list of vectors or scalars to convert to
#'
#' @return a modified data frame
#' @keywords internal
#' @export
#'
#' @section Warning: Column types (particularly of factors) are prone to change,
#'   and might drop levels. Hence using the \code{xlev} argument of
#'   \code{model.frame} is wise.
#'
#'
#' @examples
#' # changes smoking status to 'Never' for smokers who haven't already quit.
#' apply_modifications(minifhs, list(SMOKE = c("Never", "<30/day", ">=30/day")))
#' ## (note no change to rows with value "Former")
#'
#' # changes BMI status from non-overweight, regardless of previous value
#' ( beware <- apply_modifications(minifhs, list(BMI_2  = "<25.0")) )
#'
#' # beware!
#' str(minifhs$BMI_2)
#' str(beware$BMI_2)
#'
apply_modifications <- function(df, modlist) {
  for (col in names(modlist)) {
    repl <- modlist[[col]]
    if (length(repl) > 1)
      df[df[,col] %in% repl[-1], col] <- repl[1]
    else
      df[,col] <- repl
  }
  df
}

#' Convert a \code{by} object to a named list
#'
#' @param x a \code{by} object
#' @param lsep character separator between a factor name and the level
#' @param fsep character separator between two factors
#'
#' @return a named, flattened list
#' @keywords internal
#'
name_by <- function(x, lsep = ": ", fsep = ", ") {
  ## produced with the help of the `print.by` base function
  d <- dim(x)
  dn <- dimnames(x)
  dnn <- names(dn)

  newnames <- sapply(X = seq_along(x), FUN = function(i, lsep, fsep) {
    ii <- i - 1L

    paste(
      sapply(X = seq_along(dn), FUN = function(j, lsep) {
        paste0(dnn[j], lsep, dn[[j]][ii %% d[j] + 1L])
      }, lsep),
      collapse = fsep
    )
  }, lsep, fsep)

  attributes(x) <- NULL
  names(x) <- newnames
  x
}
