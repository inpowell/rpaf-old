% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paf_utils.R
\name{is_period_unsorted}
\alias{is_period_unsorted}
\title{Test if period factor is internally sorted by ID}
\usage{
is_period_unsorted(ID, PERIOD)
}
\arguments{
\item{ID}{vector of individuals' ID numbers}

\item{PERIOD}{factor of corresponding time periods}
}
\value{
\code{TRUE} if period unsorted (should raise error), \code{FALSE}
  otherwise
}
\description{
Test if period factor is internally sorted by ID
}
\keyword{internal}
