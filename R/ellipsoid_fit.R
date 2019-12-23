#' Fit ellipsoids based on distinct methods
#'
#' @description ellipsoid_fit helps in finding the centroid and matrix that
#' define an ellipsoid. It uses distinct methods with asumptions that differ
#' from each other.
#'
#' @param data data.frame or matrix of occurrence records. Columns must be
#' longitude and latitude. Other columns are optional and would represent values
#' of environmental variables to be used as dimensions for fitting the ellipsoid.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @param method (character) method to construct the ellipsoid that characterizes
#' the species ecological niche. Available methods are: "covmat", "mve1", and
#' "mve2". See details. Default = "covmat".
#' @param level (numeric) the confidence level of a pairwise confidence region
#' for the ellipsoid, expresed as percentage. Default = 95.
#' @param raster_layers optional RasterStack of environmental variables to be
#' extracted using geographic coordinates present in \code{data}. If not defined
#' \code{data} must include at least two other columns with values of the
#' variables (dimesnions) used to fit the ellipsoid.
#'
#' @return
#' An object of class \code{\link{ellipsoid}}.
#'
#' @usage
#' ellipsoid_fit(data, longitude, latitude, method = "covmat",
#'               level = 95, raster_layers = NULL)
#'
#' @details
#' The number of variables that can be used to created "ellipsoids" should be >=
#' 2 for the purposes of this package. When two variables are used "ellipsoids"
#' are ellipses, with three variables "ellipsoids" are ellipsoids, and with more
#' than three variables "ellipsoids" are hyper-ellipsoids.
#'
#' Method details are as follows:
#'
#' "covmat" creates ellipsoids based in the centriod and a matrix of covariances
#' of the variables used. A particularity of this method is that the centroid
#' will always be located in the center of the destribution of the entire data.
#' This is, the density of points in the plane of analyses matters. Analyses are
#' performed with base functions from R.
#'
#' "mve1" generates an ellipsoid that reduces the volume contained it without
#' loosing the data contained (i.e., minimum volume ellipsoid). This method may
#' modify the position of the centroid and the values in the covariance matrix
#' that can be obtained using the "covmat" method.
#'
#' "mve2" as with the previous method, this one also creates a minimum volume
#' ellipsoid. However, the algorithm for creating ellipsoids used here is
#' different. The ellipsoids created in this method are called moment based
#' minimum volume ellipsoids and, as with the previous method, they may suffer
#' changes in the position of the centroid and the values of the covariance
#' matrix if compared to the "covmat" method. In general ellipsoids created with
#' this method have smaller volumes than the ones created with previous methods.
#'
#' @export
#'
#' @examples
#' # reading data
#' occurrences <- read.csv(system.file("extdata", "occurrences.csv",
#'                                     package = "ellipsenm"))
#'
#' # raster layers of environmental data
#' vars <- raster::stack(list.files(system.file("extdata", package = "ellipsenm"),
#'                                  pattern = "bio", full.names = TRUE))
#'
#' # fitting an ellipsoid using normal covariance matrix
#' ellips <- ellipsoid_fit(data = occurrences, longitude = "longitude",
#'                         latitude = "latitude", method = "covmat",
#'                         level = 99, raster_layers = vars)
#'
#' class(ellips) # the class created in this package for this type of object
#' str(ellips)
#'
#' # using only a matrix of data and no raster layers, also another method
#' occurrences1 <- cbind(occurrences[, 2:3], raster::extract(vars, occurrences[, 2:3]))
#'
#' ellips1 <- ellipsoid_fit(occurrences1, longitude = "longitude", latitude = "latitude",
#'                          method = "mve1", level = 99)

ellipsoid_fit <- function (data, longitude, latitude, method = "covmat",
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
    data <- na.omit(raster::extract(raster_layers, data[, c(longitude, latitude)]))
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
      mvee <- mbased_mve(data, fitting_tolerance = 0.001)
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
  results <- ellipsoid(method = method,
                       centroid = centroid,
                       covariance_matrix = covari,
                       level = level * 100,
                       niche_volume = volume,
                       semi_axes_length = stds,
                       axes_coordinates = axes_coordinates)
  return(results)
}
