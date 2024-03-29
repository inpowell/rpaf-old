% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mortality_est.R
\name{mpaf_est_matrix}
\alias{mpaf_est_matrix}
\title{Calculations for mortality PAFs (deprecated)}
\usage{
mpaf_est_matrix(sr_formula, paf_data, modifications, hr_out,
  level = 0.95, ...)
}
\arguments{
\item{sr_formula}{a formula to be passed to \code{\link[survival]{survreg}}}

\item{paf_data}{an object of class \code{mpaf_response} from
\code{\link{gen_data}}}

\item{modifications}{a list of modifications to apply for PAF calculation;
see Details for more information}

\item{hr_out}{an optional character vector describing variables of
interest for hazard ratio calculation}

\item{level}{width of the confidence intervals for hazard ratios, default
0.95.}

\item{...}{extra parameters to be passed to \code{survreg}}
}
\description{
This function is deprecated, see \code{\link{est_matrix}} for current usage.
}
