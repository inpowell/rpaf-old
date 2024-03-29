% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/disease_est.R
\name{dpaf_est_diff}
\alias{dpaf_est_diff}
\title{Estimate significance of groupwise differences in disease PAFs}
\usage{
dpaf_est_diff(dpaf1, dpaf2, vv_l)
}
\arguments{
\item{dpaf1, dpaf2}{objects from \code{\link{dpaf_est_paf}} to be compared}

\item{vv_l}{named list of the covariance matrices of parameters for disease and
mortality}
}
\value{
a matrix containing the PAF differences, their standard errors, Z
  values and p values.
}
\description{
The elements of the output from this function should be able to be passed to
\code{\link[stats]{printCoefmat}}.
}
