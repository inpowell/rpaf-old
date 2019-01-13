#' Generate mortality PAF data from cohort data.
#'
#' @param indata Original data frame containing data from cohort study
#' @param id_var Name of column containing individual IDs
#' @param ft_breaks Numeric vector of breaks between time periods
#' @param death_ind Name of logical column indicating if death occurred before
#'   censoring (i.e. \code{TRUE} if an individual died before the end of
#'   follow-up)
#' @param death_time Name of numerical column with individual time until death
#'   or censoring
#' @param na.action Function (found via call to \code{match.fun}) which
#'   indicates how this function should treat missing data
#' @param variables Character vector of predictor columns to keep from
#'   \code{indata}
#' @param period_factor,time_var Names of new columns in data frame representing
#'   individual follow-up period and follow-up time respectively
#'
#' @return A list of class \code{"mpaf_data"} containing these elements:
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
#' @export
mpaf_gen_data <- function(
  indata,
  id_var,
  ft_breaks,
  death_ind,
  death_time,
  variables = character(0),
  na.action = getOption("na.action"),
  period_factor = "f_period",
  time_var = "f_end"
) {
  if (is.unsorted(ft_breaks))
    stop("Time breaks are not in order.")

  mpaf <- list(data_call = match.call())

  # drop columns
  # select columns
  if (missing(death_time) || missing(death_ind))
    df <- indata[,c(id_var, variables), drop = FALSE]
  else
    df <- indata[,c(id_var, variables, death_time, death_ind)]

  # handle missing data
  df <- match.fun(na.action)(df)
  mpaf$nobs <- nrow(df)
  mpaf$na.action <- stats::na.action(df)

  # time break information
  mpaf$breaks <- ft_breaks

  # period information in data frame
  # integer codes as in findInterval -- note not a factor yet
  intervals <- data.frame(1:(length(ft_breaks) - 1))
  names(intervals) <- period_factor

  ldf <- merge(intervals, df, by = NULL)

  if (!missing(death_ind) && !missing(death_time)) {
    ldf[[time_var]] <- ldf[[death_time]]
    ft_intervals <- findInterval(ldf[[time_var]], ft_breaks)
    ldf[ldf[[period_factor]] > ft_intervals, c(time_var, death_ind)] <- NA
    ldf[ldf[[period_factor]] < ft_intervals, death_ind] <- FALSE
    ldf[ldf[[period_factor]] < ft_intervals, time_var] <-
      subset(diff(ft_breaks)[ldf[[period_factor]]],
             ft_intervals > ldf[[period_factor]])
    ldf[ft_intervals == ldf[[period_factor]], time_var] <-
      subset(ldf[[time_var]] - ft_breaks[ft_intervals],
             ft_intervals == ldf[[period_factor]])

    ldf <- subset(ldf, select = setdiff(names(ldf), death_time))
  }

  # convert period column to factor
  ldf[[period_factor]] <- factor(
    ldf[[period_factor]],
    labels = paste0("(", head(ft_breaks, -1), ",", tail(ft_breaks, -1), "]")
  )

  mpaf$ID <- ldf[[id_var]]
  mpaf$PERIOD <- ldf[[period_factor]]

  mpaf$data <- subset(ldf, select = setdiff(names(ldf), c(id_var)))

  class(mpaf) <- "mpaf_data"

  if (!missing(death_ind) && !missing(death_time))
    class(mpaf) <- c(class(mpaf), "mpaf_response")

  mpaf
}
