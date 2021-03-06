% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_classes.R
\docType{class}
\name{ellipsoid}
\alias{ellipsoid}
\alias{ellipsoid-class}
\alias{ellipsoid_model_sim-class}
\alias{ellipsoid_model_rep-class}
\alias{ellipsoid_model_sim}
\alias{ellipsoid_model_rep}
\title{An S4 class to organize data and results of ellipsoid* objects}
\description{
An S4 class to organize data and results of ellipsoid* objects
}
\section{Slots}{

\describe{
\item{\code{method}}{("character"), not in ellipsoid_model_rep.}

\item{\code{level}}{("numeric"), not in ellipsoid_model_rep.}

\item{\code{centroid}}{("numeric"), not in ellipsoid_model_rep.}

\item{\code{covariance_matrix}}{("matrix"), not in ellipsoid_model_rep.}

\item{\code{niche_volume}}{("numeric"), not in ellipsoid_model_rep.}

\item{\code{semi_axes_length}}{("numeric"), not in ellipsoid_model_rep.}

\item{\code{axes_coordinates}}{list of matrices}

\item{\code{ellipsoids}}{list of objects of class ellipsoid, only ellipsoid_model_rep.}

\item{\code{mahalanobis}}{if ellipsoid_model_sim, ("numeric"); if ellipsoid_model_rep,
numeric matrix.}

\item{\code{suitability}}{if ellipsoid_model_sim, ("numeric"); if ellipsoid_model_rep,
numeric matrix.}

\item{\code{prevalence}}{if ellipsoid_model_sim, ("numeric"); if ellipsoid_model_rep,
numeric matrix.}

\item{\code{prediction_maha}}{object of class Raster*.}

\item{\code{prediction_suit}}{object of class Raster*.}

\item{\code{mahalanobis_proj}}{if ellipsoid_model_sim, numeric matrix; if
ellipsoid_model_rep, list of numeric matrices.}

\item{\code{suitability_proj}}{if ellipsoid_model_sim, numeric matrix; if
ellipsoid_model_rep, list of numeric matrices.}

\item{\code{projections_maha}}{if ellipsoid_model_sim, object of class Raster*, if
ellipsoid_model_rep, list of Raster* objects.}

\item{\code{projections_suit}}{if ellipsoid_model_sim, object of class Raster*, if
ellipsoid_model_rep, list of Raster* objects.}
}}

\examples{
showClass("ellipsoid")
#' showClass("ellipsoid_model_rep")
}
