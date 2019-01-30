#' Generate mortality PAF data from cohort data.
#'
#' @param indata Original data frame containing data from cohort study
#' @param id_var Name of column containing individual IDs
#' @param ft_breaks Numeric vector of breaks between time periods
#' @param death_ind (Optional) Name of logical column indicating if death
#'   occurred before censoring (i.e. \code{TRUE} if an individual died before
#'   the end of follow-up)
#' @param death_time (Optional) Name of numerical column with individual time
#'   until death or censoring
#' @param disease_ind (Optional) Name of logical column indicating if disease
#'   occurred before censoring
#' @param disease_time (Optional) Name of numerical column with individual time
#'   until disease or censoring
#' @param na.action Function (found via call to \code{match.fun}) which
#'   indicates how this function should treat missing data
#' @param variables Character vector of predictor columns to keep from
#'   \code{indata}
#' @param period_factor,time_var Names of new columns in data frame representing
#'   individual follow-up period and follow-up time respectively
#'
#' @return A list of class \code{"paf_data"} containing these elements:
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
#'   \item{data}{A data frame containing predictor (and possibly response)
#'   columns as given in the \code{variables} argument, as well as a column for
#'   the period (if requested)}
#'
#'   If a response is included, the returned list also inherits
#'   \code{"dpaf_response"} or \code{"mpaf_response"} for disease or mortality
#'   PAF calculations respectively.
#'
#' @export
gen_data <- function(
  indata,
  id_var,
  ft_breaks,
  death_ind,
  death_time,
  disease_ind,
  disease_time,
  variables = character(0),
  na.action = getOption("na.action"),
  period_factor = "f_period",
  time_var = "f_end"
) {
  if (is.unsorted(ft_breaks)) {
    warning("Time breaks are not in order; reordering them.")
    ft_breaks <- sort(ft_breaks)
  }

  # first include the function call for user reference later
  paf_data <- list(data_call = match.call())

  # check to see if we have enough information for a disease PAF (1st case), a
  # mortality PAF (2nd case) or just predictors (3rd case)
  if (!missing(death_time) && !missing(death_ind) &&
      !missing(disease_time) && !missing(disease_ind))
    type <- c("paf_data", "dpaf_response")
  else if (!missing(death_time) && !missing(death_ind))
    type <- c("paf_data", "mpaf_response")
  else
    type <- c("paf_data")

  # select columns we want to keep from the data frame, which depends on the
  # information we have
  if ("dpaf_response" %in% type)
    df <- indata[,c(id_var, variables, death_time, death_ind,
                    disease_time, disease_ind)]
  else if ("mpaf_response" %in% type)
    df <- indata[,c(id_var, variables, death_time, death_ind)]
  else
    df <- indata[,c(id_var, variables), drop = FALSE]

  # handle missing data and include actions
  df <- match.fun(na.action)(df)
  paf_data$nobs <- nrow(df)
  paf_data$na.action <- stats::na.action(df)

  # include time break information
  paf_data$breaks <- ft_breaks

  # use integer codes compatible with findInterval to initially represent
  # period, and name the small, one-column data.frame appropriately
  intervals <- data.frame(1:(length(ft_breaks) - 1))
  names(intervals) <- period_factor

  # create a long data frame by combinatorial merge. This replicates each row
  # for each period, with a period indicator in the column named by
  # period_factor
  ldf <- merge(intervals, df, by = NULL)

  if ("mpaf_response" %in% type) {
    # move death_time column to the new time_var column
    ldf[[time_var]] <- ldf[[death_time]]

    # figure out which interval death (or censoring!) lies in
    ft_intervals <- findInterval(ldf[[time_var]], ft_breaks)

    # "partial time" to event is undefined after event has happened -- set to
    # missing in this case
    ldf[ldf[[period_factor]] > ft_intervals, c(time_var, death_ind)] <- NA

    # before event, death  hasn't happened
    ldf[ldf[[period_factor]] < ft_intervals, death_ind] <- FALSE

    # "partial times" if no event has happened is defined as the length of the
    # period
    ldf[ldf[[period_factor]] < ft_intervals, time_var] <-
      subset(diff(ft_breaks)[ldf[[period_factor]]],
             ft_intervals > ldf[[period_factor]])

    # partial times when event has happened needs to be measured from the start
    # of the period
    ldf[ft_intervals == ldf[[period_factor]], time_var] <-
      subset(ldf[[time_var]] - ft_breaks[ft_intervals],
             ft_intervals == ldf[[period_factor]])

    # remove original death_time column from the data frame
    ldf <- subset(ldf, select = setdiff(names(ldf), death_time))
  }
  if ("dpaf_response" %in% type) {
    # the time until censoring is the minimum time until event
    ldf[[time_var]] <- pmin(ldf[[death_time]], ldf[[disease_time]])

    # note that breaks[data[[period_factor]]] gives the time at start of period
    # and breaks[-1][data[[period_factor]]] gives time at end of period
    ft_intervals <- findInterval(ldf[[time_var]], ft_breaks, left.open = TRUE)
    post_end <- ft_intervals < ldf[[period_factor]]
    pre_end <- ft_intervals > ldf[[period_factor]]
    is_end <- ft_intervals == ldf[[period_factor]]

    ldf[[death_ind]][post_end] <- ldf[[disease_ind]][post_end] <- NA
    ldf[[death_ind]][pre_end] <- ldf[[disease_ind]][pre_end] <- FALSE
    ldf[[death_ind]][is_end] <- ldf[[death_ind]][is_end] &
      ldf[[death_time]][is_end] <= ldf[[disease_time]][is_end]
    ldf[[disease_ind]][is_end] <- ldf[[disease_ind]][is_end] &
      ldf[[disease_time]][is_end] <= ldf[[death_time]][is_end]

    # adjust times:
    #  - after event, partial time is undefined
    #  - before event, partial time is the length of the period of interest
    #  - when event occurs, need to measure time from the start of the period
    ldf[[time_var]][post_end] <- NA
    ldf[[time_var]][pre_end] <- diff(ft_breaks)[ldf[pre_end, period_factor]]
    ldf[[time_var]][is_end] <- ldf[is_end, time_var] -
      ft_breaks[ft_intervals[is_end]]
  }

  # convert period column to factor, with useful names
  ldf[[period_factor]] <- factor(
    ldf[[period_factor]],
    labels = paste0("(", utils::head(ft_breaks, -1), ",",
                    utils::tail(ft_breaks, -1), "]")
  )

  # store ID and PERIOD columns separately, because they may not be used in the
  # regression model
  paf_data$ID <- ldf[[id_var]]
  paf_data$PERIOD <- ldf[[period_factor]]

  # remove ID column from data frame
  paf_data$data <- subset(ldf, select = setdiff(names(ldf), c(id_var)))

  class(paf_data) <- type
  paf_data
}

#' @export
#' @rdname gen_data
mpaf_gen_data <- function(indata, id_var, ft_breaks, death_ind, death_time,
              variables = character(0), na.action = getOption("na.action"),
              period_factor = "f_period", time_var = "f_end") {
  .Deprecated("gen_data")
  gen_data(
    indata = indata,
    id_var = id_var,
    ft_breaks = ft_breaks,
    death_ind = death_ind,
    death_time = death_time,
    variables = variables,
    na.action = na.action,
    period_factor = period_factor,
    time_var = time_var
  )
}
