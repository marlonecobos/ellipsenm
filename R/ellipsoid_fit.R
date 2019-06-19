#' Fit ellipsoids based on distinct methods
#'
#' @description ellipsoid_fit helps in finding the centroid and matrix that
#' define an ellipsoid. It uses distinct methods with asumptions that differ
#' from each other.
#'
#' @param data data.frame of occurrence records. Columns must be longitude and
#' latitude. Other columns are optional and could values of environmental
#' variables to be used as dimensions for fitting the ellipsoid.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @param method (character) method to construct the ellipsoid that characterizes
#' the species ecological niche. Available methods are: "covmat", "mve1", and
#' "mve2". See details. Default = "mve1".
#' @param level (numeric) percentage of data to be considered when creating the
#' ellipsoid that characterizes the species ecological niche. Default = 95.
#' @param raster_layers optional RasterStack of environmental variables to be
#' extracted using geographic coordinates present in \code{data}.
#'
#' @return
#'
#' @details
#' Methods details are as follows:
#'
#' "covmat"
#'
#' "mve1"
#'
#' "mve2"
#'
#' @export
#'
#' @examples
#' occurrences <- read.csv(system.file("extdata", "occurrences_comp.csv",
#'                                     package = "ellipsenm"))[, -1]
#'
#' vars <- raster::stack(list.files(system.file("extdata", package = "ellipsenm"),
#'                                  pattern = "m_bio", full.names = TRUE))
#'
#' ellips <- ellipsoid_fit(data = occurrences, longitude = "longitude",
#'                         latitude = "latitude", method = "mve1",
#'                         level = 95, raster_layers = vars)

ellipsoid_fit <- function (data, longitude, latitude, method = "mve1",
                           level = 95, raster_layers = NULL) {

  # -----------
  # detecting potential errors
  if (missing(data)) {
    stop("Argument data is necessary to perform the analysis")
  }
  if (missing(longitude)) {
    stop("Argument longitude is not defined.")
  }
  if (missing(latitude)) {
    stop("Argument latitude is not defined.")
  }

  # -----------
  # preparing data
  if (!is.null(raster_layers)) {
    data <- raster::extract(raster_layers, data[, c(longitude, latitude)])
  } else {
    data <- data[, -c(longitude, latitude)]
  }

  # -----------
  # finding centroid and covariance matrix
  if (method == "covmat" | method == "mve1" | method == "mve2") {
    if (method == "covmat") {
      centroid <- colMeans(data)
      covari <- stats::cov(data)
    }

    if (method == "mve1") {
      n <- ndata_quantile(nrow(data), level)
      cent_var <- MASS::cov.mve(data, quantile.used = n)
      centroid <- cent_var$center
      covari <- cent_var$cov
    }

    if (method == "mve2") {
      mvee <- mbased_mve(data, tolerance = 0.001)
      centroid <- mvee[[1]]
      covari <- mvee[[2]] # check if it works
    }
  } else {
    stop("Argument method is not valid, please see function's help.")
  }

  # -----------
  # calculating ellipsoid characteristics
  sigma_i <- solve(covari) / stats::qchisq(level, df = ncol(data))
  s_eigen <- eigen(sigma_i)
  s_eigenval <- s_eigen$values
  s_eigenvec <- s_eigen$vectors
  stds <- 1 / sqrt(s_eigenval)
  axes_length <- NULL

  for (i in 1:dim(sigma_i)[1]) {
    axes_length[i] <- stds[i] * 2
  }

  names(axes_length) <- letters[1:dim(covari)[1]]
  n_dimensions <- dim(covari)[1]
  volume <- ellipsoid_volume(n_dimensions, axes_length / 2)
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

  # -----------
  # preparing results
  results <- ellipsoid(method = method,
                       centroid = centroid,
                       covariance_matrix = covari,
                       niche_volume = volume,
                       semi_axes_length = axes_length / 2,
                       axis_coordinates = axis_coordinates)
  return(results)
}


#' Helper function to calculate niche volume
#'
#' @param n_dimensions (numeric) number of dimensions to be considered.
#' @param axes_length (numeric) length of ellipsoid axes.
#'
#' @export

ellipsoid_volume <- function (n_dimensions, axes_length) {
  term1 <- 2 * pi^(n_dimensions / 2)
  term2 <- n_dimensions * gamma(n_dimensions / 2)
  term3 <- prod(axes_length)
  term4 <- (term1 / term2) * term3

  return(term4)
}


#' Helper function to calculate quantiles
#'
#' @param n_data (numeric) total number of data to be considered.
#' @param level (numeric) percentage of data to be considered when creating the
#' ellipsoid that characterizes the species ecological niche. Default = 95.
#'
#' @export

ndata_quantile <- function(n_data, level) {
  n <- floor(n_data * level)
  if (n > n_data) {n <- n_data}
  return(n)
}
