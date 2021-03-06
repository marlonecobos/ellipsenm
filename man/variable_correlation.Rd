% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variable_correlation.R
\name{variable_correlation}
\alias{variable_correlation}
\title{Evaluates correlation among variables}
\usage{
variable_correlation(variables, sample_size = 10000, correlation_limit = 0.8,
                     save = FALSE, name = "correlation", corrplot = FALSE,
                     magnify_to = NULL, ...)
}
\arguments{
\item{variables}{RasterStack, RasterBrick, or matrix. If matrix, columns
represent distinct variables for analysis, otherwise, group of raster layers.}

\item{sample_size}{(numeric) sample size to be taken from all variables;
default = 10000.}

\item{correlation_limit}{(numeric) absolute value of correlation limit;
default = 0.8.}

\item{save}{(logical) whether or not to save the results; default = FALSE.}

\item{name}{(character) name of the csv files to be writen;
default = "correlation".}

\item{corrplot}{(logical) whether or not to plot the results; default = FALSE.}

\item{magnify_to}{(numeric) optional value to be used to magnify all values
with absolute correlations above \code{correlation_limit}. Default = NULL.}

\item{...}{other arguments to be passed to \code{\link[corrplot]{corrplot}}.
Arguments "type", "tl.col", and "tl.srt" are fixed.}
}
\value{
A correlation matrix. If argument \code{corrplot} = TRUE correlation values
are shown in a plot.
}
\description{
variable_correlation helps in evaluating correlation among
distinct variables.
}
\details{
If \code{magnify_to} is defined and \code{save} = TRUE, an additional csv
file named as "\code{name}_magnified.csv" will be written.
}
\examples{
# raster layers of environmental data
vars <- raster::stack(list.files(system.file("extdata", package = "ellipsenm"),
                                 pattern = "bio", full.names = TRUE))

# simple correlation matrix
cors <- variable_correlation(variables, sample_size = 5000)

# correlation matrix and plot (values correlated above |0.8| are magnified)
cors <- variable_correlation(variables, sample_size = 5000, corrplot = TRUE,
                             magnified = 2)

# to save results check arguments "save" and "name"
}
