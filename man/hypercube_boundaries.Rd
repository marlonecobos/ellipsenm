% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/short_helpers.R
\name{hypercube_boundaries}
\alias{hypercube_boundaries}
\title{Helper funtion to get data for Montecarlo simulations}
\usage{
hypercube_boundaries(ellipsoids, n_points = 1e+06)
}
\arguments{
\item{ellipsoids}{list of ellipsoid objects.}

\item{n_points}{(numeric) number of random points to be generated.}
}
\value{
A total of \code{n_points} created randomly considering the limits of the
ellipsoids considered.
}
\description{
Helper funtion to get data for Montecarlo simulations
}
