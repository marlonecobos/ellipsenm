#' Fit ellipsoids based on distinct methods
#'
#' @description
#'
#' @param data data.frame of occurrence records. Columns must be species,
#' longitude, and latitude. Other columns are values of environmental variables
#' to be used as dimensions of the ellipsoid.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @param method (character) method to construct the ellipsoid that characterizes
#' the species ecological niche. Available methods are: "covmat" and "mve".
#' See details. Default = "mve".
#' @param level (numeric) percentage of data to be considered when creating the
#' ellipsoid that characterizes the species ecological niche.
#' @param raster_layers optional RasterStack of environmental variables to be
#' extracted using geographic coordinates present in \code{data}.
#'
#' @return
#'
#'
#' @export
#'
#' @examples
#' occurrences <- read.csv(system.file("extdata", "occurrences_comp.csv",
#'                                     package = "ellipsenm"))
#'
#' vars <- raster::stack(list.files(system.file("extdata", package = "ellipsenm"),
#'                                  pattern = "m_bio", full.names = TRUE))

ellipsoid_fit <- function (data, longitude, latitude, method = "mve", level,
                           raster_layers = NULL) {

  if (!is.null(raster_layers)) {
    data <- raster::extract(raster_layers, data[, c(longitude, latitude)])
  } else {
    data <- data[, -c(longitude, latitude)]
  }

  if (method = "mve") {
    n <- ndata_quantile(nrow(data), level)
    cent_var <- MASS::cov.mve(data, quantile.used = n)
    centroid <- cent_var$center
    covari <- cent_var$cov
  } else {
    centroid <- colMeans(data)
    covari <- stats::cov(data)
  }

  sigma_i <- solve(covari) / stats::qchisq(level, df = ncol(data))
  s_eigen <- eigen(sigma_i)
  s_eigenval <- s_eigen$values
  s_eigenvec <- s_eigen$vectors
  stds <- 1 / sqrt(s_eigenval)
  axis_length <- NULL

  for (i in 1:dim(sigma_i)[1]) {
    axis_length[i] <- stds[i] * 2
  }

  names(axis_length) <- letters[1:dim(covari)[1]]
  n <- dim(covari)[1]
  vol2 <- ellipsoid_volume(n, axis_length / 2)
  axis_coordinates <- list()

  for (i in 1:dim(covari)[1]) {
    assign(paste0("l", i, "_inf"), centroid - s_eigenvec[, i] * stds[i])
    assign(paste0("l", i, "_sup"), centroid + s_eigenvec[, i] * stds[i])
    coord_matrix <- matrix(c(eval(parse(text = paste0("l", i, "_sup"))),
                             eval(parse(text = paste0("l", i, "_inf")))),
                           byrow = TRUE, nrow = 2)
    colnames(coord_matrix) <- names(centroid)
    rownames(coord_matrix) <- paste0("vec_", 1:2)
    axis_coordinates[[i]] <- coord_matrix
  }

  results <- ellipsoid(method = method,
                       centroid = centroid,
                       covariance_matrix = covari,
                       niche_volume = vol2,
                       semi_axis_length = axis_length / 2,
                       axis_coordinates = axis_coordinates)
  return(results)
}


#' Helper function to calculate niche volume
#'

ellipsoid_volume <- function (n, axis_length) {
  term1 <- 2 * pi^(n / 2)
  term2 <- n * gamma(n / 2)
  term3 <- prod(axis_length)
  term4 <- (term1 / term2) * term3

  return(term4)
}


#' Helper function to calculate quantiles

ndata_quantile <- function(n_data, level) {
  n <- floor(n_data * level)
  if (n > n_data) {n <- n_data}
  return(n)
}
