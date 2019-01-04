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

#' Create design frames from data
#'
#' @param df initial data frame
#' @param terms object of terms with response dropped
#' @param modifications list of modifications to make (see
#'   \code{apply_modifications} for format)
#' @param xlev xlevels object from the survival regression
#'
#' @return a list of two model frames, \code{design} and \code{modified},
#'   carrying predictor values
#' @keywords internal
design_frames <- function(df, terms, modifications, xlev) {
  design <- stats::model.frame(terms, df, na.action = na.pass, xlev = xlev)
  mdf <- apply_modifications(df, modifications)
  modified <- stats::model.frame(terms, mdf, na.action = na.pass, xlev = xlev)

  list(design = design, modified = modified)
}
