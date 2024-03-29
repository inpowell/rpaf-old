% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mortality_calc.R
\name{mpaf_S}
\alias{mpaf_S}
\title{Calculate survivals for mpaf study}
\usage{
mpaf_S(lambda, ID, PERIOD, dt)
}
\arguments{
\item{lambda}{vector of hazards corresponding to \code{ID}s and
\code{PERIOD}s}

\item{ID}{vector of IDs}

\item{PERIOD}{vector of periods}

\item{dt}{delta-times -- length of each period}
}
\value{
named matrix of survivals
}
\description{
Assumes arguments are sorted -- check with \code{\link{is_period_unsorted}}.
}
\keyword{internal}
