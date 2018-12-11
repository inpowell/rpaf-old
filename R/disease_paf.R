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
#' @examples
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

  # select columns ----------------------------------------------------------
  if (missing(death_time) || missing(death_ind) ||
      missing(disease_time) || missing(disease_ind))
    data <- rawdata[,c(id_var, variables), drop = FALSE]
  else
    data <- rawdata[,c(id_var, variables, death_time, death_ind,
                       disease_time, disease_ind)]

  # handle missing values ---------------------------------------------------
  data <- match.fun(na.action)(data)
  retlist <- c(retlist, list("na.action" = na.action(data)))


  # add in time breaks ------------------------------------------------------

  retlist <- c(retlist, list("breaks" = time_breaks))

  # lengthen data -----------------------------------------------------------
  # integer codes as in findInterval -- note not a factor yet
  intervals <- data.frame(1:(length(time_breaks) - 1))
  names(intervals) <- period_factor

  data <- merge(intervals, data, by = NULL)

  # correct indicators and times ---------------------------------------------

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


  # set period to factor ----------------------------------------------------

  data[[period_factor]] <- factor(data[[period_factor]],
                                  labels = levels(cut(numeric(0), time_breaks)))


  # add data to list --------------------------------------------------------

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

#' Fit PAFs for disease data
#'
#' @param disease_resp formula with response for disease regression (note RHS
#'   should be \code{.})
#' @param death_resp formula with response for mortality regression (note RHS
#'   should be \code{.})
#' @param predictors formula with RHS containing predictors for both survival
#'   regressions (LHS should be empty)
#' @param modlist list of modifications for PAF calculation (see
#'   \code{\link{apply_modifications}} and \code{\link{mpafreg}} for more
#'   information)
#' @param data_list a list containing data frames and other information, from
#'   \code{\link{dpaf_data}}
#' @param ft_length length of follow-up time
#' @param ft_delta length of periods of constant hazard
#' @param id_var name of column containing individual ID, as for
#'   \code{\link{dpaf_data}}
#' @param period_var name of column containing period factor, as for
#'   \code{\link{dpaf_data}}
#' @param ... other arguments to be passed to \code{\link[survival]{survreg}}.
#'
#' @return
#' @export
#'
#' @examples
dpaf <- function(disease_resp, death_resp, predictors,
                 modlist, data_list, ft_length, ft_delta,
                 id_var, period_var, ...) {
}

#' Calculate gradients of morbidity for PAF calculations
#'
#' @param x
#' @param hz_d
#' @param hz_m
#' @param sv_d
#' @param sv_m
#' @param breaks
#' @param ID
#' @param PERIOD
#'
#' @return
#' @export
#' @keywords internal
dpaf_grad <- function(x, hz_d, hz_m, sv_d, sv_m, breaks, ID, PERIOD) {
}

#' Title
#'
#' @param object
#' @param parm
#' @param level
#' @param ...
#'
#' @return
#' @export
#' @method confint dpaf
#' @keywords internal
#'
#' @examples
confint.dpaf <- function(object, parm, level = 0.95, ...) {
}
