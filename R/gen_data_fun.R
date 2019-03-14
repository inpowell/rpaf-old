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
#' @examples
#' # The original data frame we will work with has two time
#' # columns -- one for diabetes incidence (DIAB_FT), and one
#' # for death (DEATH_FT).
#' head(minifhs)
#'
#' # we are interested in the first of diabetes incidence and death
#' ct_minifhs <- collapse_times(minifhs, "DIAB", "DIAB_FT", "DEATH", "DEATH_FT")
#' head(col_minifhs)
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
