
# Data preparation utilities ----------------------------------------------

#' Prepare data for PAF disease studies
#'
#' @param rawdata a data frame containing
#' @param variables character vector of variables of interest
#' @param time_breaks numeric vector of breaks between periods of constant
#'   hazard
#' @param id_var column of data frame identifying individuals (will be generated
#'   sequentially if missing)
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

  if (missing(id_var)) {
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
  retlist <- c(retlist, list("na.action" = na.action(data)))


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
#'   \code{dpaf_data}, plus the following elements:
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
dpaf_split <- function(object, INDICES, lsep = ": ", vsep = ", ") {
  split_call <- match.call()

  if (!is.null(names(INDICES))) {
    INDICES <- mapply(function(fctr, nm) {
      levels(fctr) <- paste0(nm, lsep, levels(fctr))
      fctr
    }, INDICES, names(INDICES), SIMPLIFY = FALSE)
  }
  IDs <- split(object$ID, INDICES, sep = vsep)
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
#'   \code{\link{dpaf_spl}}
#' @param ... other arguments to be passed to \code{predict.dpaf}
#'
#' @return a list of two elements
#'
#'   \item{group_pafs}{calculated disease PAFs, as well as other information
#'   output by \code{predict.dpaf}}
#'
#'   \item{paf_diffs}{differences between PAFs, with their standard errors, test
#'   statistics, and p-values}
#'
#' @export
dpaf_groups <- function(object, modlist, split_data, ...) {
  stopifnot(inherits(object, "dpaf"))
  stopifnot(any(sapply(inherits, split_data, "dpaf_split")))

  gwise_preds <- lapply(split_data, stats::predict, object = object,
                        modlist = modlist, type = 'paf', gradient = TRUE, ...)

  pafs <- do.call(rbind, lapply(gwise_preds, `[[`, "estimate"))
  # row.names(pafs) <- names(gwise_preds)

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

  list(
    group_pafs = pafs,
    paf_diffs = do.call(rbind, diffs)
  )
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
