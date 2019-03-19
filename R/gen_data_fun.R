#' Generate time-varying cohort data
#'
#' This function takes time-to-event data and converts it to a format that
#' varies with time. The original data consist of a data frame with covariate
#' columns and response columns. The response columns should include one
#' numerical column for time, plus one logical column for each event indicating
#' if this event was the first to occur.
#'
#' @param indata data frame to be extended in time
#' @param ft_breaks numeric vector of breaks between time periods
#' @param time name of response column for time
#' @param primary name of response column for event of interest
#' @param secondary names of response columns for competing risk events
#' @param id_col optional name of column containing unique identifiers for each
#'   row
#' @param period_factor name of new column containing time period
#'
#' @return a data frame with the at least the \code{ft_breaks}, \code{time},
#'   \code{primary}, \code{secondary} and \code{period_factor} arguments as
#'   attributes. Further, the data frame has an \code{id} attribute that gives
#'   unique values for each individual, which may have been constructed if
#'   \code{id_col} was not given.
#'
#'   The class of the data frame has been prepended with \code{paf_data} if
#'   there are no competing events, or \code{cr_paf_data} if potential competing
#'   risks hav ebeen included.
#'
#' @export
#'
#' @example exec/minifhs-prep.R
gen_data_fun <- function(indata, ft_breaks, time,
                         primary, secondary = NULL,
                         id_col, period_factor = "period") {
  # check if id_col is given
  if (missing(id_col)) {
    id_col <- make.unique(c(names(indata), "ID"))[ncol(indata) + 1]
    indata[[id_col]] <- seq_len(nrow(indata))
  } else if (!id_col %in% names(id_col)) {
    indata[[id_col]] <- seq_len(nrow(indata))
  } else if (any(duplicated(indata[[id_col]]))) {
    stop("Duplicated ID values in indata")
  }

  # lengthen data frame by time
  ldf <- merge(
    structure(data.frame(1:(length(ft_breaks) - 1)), names = period_factor),
    indata
  )

  time_int <- findInterval(ldf[[time]], ft_breaks, rightmost.closed = TRUE)

  # set indicators and times to missing if event occurs afterwards
  ldf[ldf[[period_factor]] > time_int, c(primary, secondary, time)] <- NA
  # ... to FALSE/dt if event occurs before
  event_before <- ldf[[period_factor]] < time_int
  ldf[event_before, c(primary, secondary)] <- FALSE
  ldf[event_before, time] <- diff(ft_breaks)[ldf[event_before, period_factor]]

  # set time to partial time under which event occurs
  event_occurs <- ldf[[period_factor]] == time_int
  ldf[event_occurs, time] <- ldf[event_occurs, time] - ft_breaks[time_int[event_occurs]]

  # convert period factor to actual factor
  ldf[[period_factor]] <- factor(
    ldf[[period_factor]], levels = 1:(length(ft_breaks) - 1),
    labels = paste0("(", ft_breaks[-length(ft_breaks)], ",", ft_breaks[-1], "]")
  )

  structure(
    ldf,
    class = c(if (is.null(secondary)) NULL else "cr_paf_data", "paf_data", class(ldf)),
    ft_breaks = ft_breaks,
    primary = primary,
    secondary = secondary,
    time = time,
    period = period_factor,
    id = id_col
  )
}

#' Collapse time-to-event data from multiple events
#'
#' This function converts multiple time-to-event columns into one time column,
#' with indicator columns that are \code{TRUE} only for the first event (or
#' first events, if these occur simultaneously).
#'
#' @param indata data frame to be modified
#' @param ind1,time1,... columns representing indicator and time data
#'   (respectively) for each event, repeated for each event
#' @param time.out name of the output time column
#'
#' @return a data frame with a single column of event times, with modified
#'   indicators such that only the first events are recorded
#' @export
#'
#' @example exec/minifhs-prep.R
collapse_times <- function(indata, ind1, time1, ..., time.out = "time") {
  xargs <- c(ind1, time1, ...) # ellipsis arguments as character vectors

  if (length(xargs) %% 2 != 0) # should have an even number of arguments
    stop("Wrong number of secondary column names")

  inds <- xargs[seq(1, length(xargs), 2)]
  times <- xargs[seq(2, length(xargs), 2)]
  time_data <- indata[, times, drop = FALSE]

  # identify time indices for earliest event in each row
  first_times <- apply(time_data, 1, min)

  # set all other indicator columns to FALSE
  indata[, inds][time_data > first_times] <- FALSE

  cbind(
    structure(data.frame(first_times), names = time.out),
    indata[, !names(indata) %in% times, drop = FALSE]
  )
}

#' Expand time-to-event data from multiple events
#'
#' This function converts an event factor from time-to-event data into multiple
#' indicator columns. The indicator columns are \code{TRUE} if the corresponding
#' event was recorded.
#'
#' @param indata data frame to be modified
#' @param event column of factor describing the event that occurred
#' @param prefix prefix to append to new column names
#'
#' @return a data frame with the original event column removed, and replaced
#'   with an indicator column for each level of the event factor
#' @export
#'
#' @examples
#' # Our original data frame with liver transplant data has one
#' # of four events occurring, stored in the "event" column.
#' head(survival::transplant)
#'
#' # We want to convert the event column into multiple indicator
#' # columns
#' tsplt_ind <- expand_events(survival::transplant, "event")
#' head(tsplt_ind)
expand_events <- function(indata, event, prefix = paste0(event, ".")) {
  lv <- levels(indata[[event]])
  new_cols <- paste0(prefix, lv)

  tf_event <- structure(lapply(lv, `==`, indata[[event]]),
                        names = new_cols)
  tf_cols <- do.call(cbind, tf_event)

  cbind(
    tf_cols,
    indata[, !names(indata) == event, drop = FALSE]
  )
}
