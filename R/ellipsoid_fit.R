#' Fit ellipsoids based on distinct methods
#'
#' @description ellipsoid_fit helps in finding the centroid and matrix that
#' define an ellipsoid. It uses distinct methods with asumptions that differ
#' from each other.
#'
#' @param data data.frame or matrix of occurrence records. Columns must include
#' longitude and latitude. Other columns are optional and could values of
#' environmental variables to be used as dimensions for fitting the ellipsoid.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @param method (character) method to construct the ellipsoid that characterizes
#' the species ecological niche. Available methods are: "covmat", "mve1", and
#' "mve2". See details. Default = "mve1".
#' @param level (numeric) the confidence level of a pairwise confidence region
#' for the ellipsoid, expresed as percentage. Default = 95.
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
#'
#' # using only a matrix of data and no raster layers
#' data1 <- cbind(occurrences, raster::extract(ras, occurrences))
#'
#' ell3 <- ellipsoid_fit(data1, longitude = "longitude", latitude = "latitude",
#'                       method = "mve1", level = 95)

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
    data <- data[, -which(colnames(data) %in% c(longitude, latitude))]
  }
  level <- level / 100

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
      covari <- mvee[[2]] # PROBLEM, see if level needs to be changed
    }
  } else {
    stop("Argument method is not valid, please see function's help.")
  }

  # -----------
  # calculating ellipsoid characteristics
  ndim <- ncol(data)
  sigma_i <- solve(covari) / stats::qchisq(level, df = ndim)
  s_eigen <- eigen(sigma_i)
  s_eigenval <- s_eigen$values
  s_eigenvec <- s_eigen$vectors

  ## semi axes length
  stds <- 1 / sqrt(s_eigenval)
  names(stds) <- letters[1:ndim]

  ## volume
  volume <- ellipsoid_volume(ndim, stds)

  ## axes coordinates
  axes_coordinates <- lapply(1:ndim, function(x) {
    coor <- rbind((centroid + s_eigenvec[, x] * stds[x]),
                  (centroid - s_eigenvec[, x] * stds[x]))
    rownames(coor) <- paste0("vec_", 1:2)
    return(coor)
  })
  names(axes_coordinates) <- names(stds)

  # -----------
  # preparing results
  results <- ellipsoid_basic(method = method,
                             centroid = centroid,
                             covariance_matrix = covari,
                             level = level * 100,
                             niche_volume = volume,
                             semi_axes_length = stds,
                             axes_coordinates = axes_coordinates)
  return(results)
}


#' Helper function to calculate niche volume
#'
#' @param n_dimensions (numeric) number of dimensions to be considered.
#' @param semi_axes_length (numeric) length of ellipsoid axes.
#'
#' @export

ellipsoid_volume <- function (n_dimensions, semi_axes_length) {
  term1 <- 2 * pi^(n_dimensions / 2)
  term2 <- n_dimensions * gamma(n_dimensions / 2)
  term3 <- prod(semi_axes_length)
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
