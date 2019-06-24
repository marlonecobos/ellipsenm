#' Mean, min, and max ellipsoids
#'
#' @description mmm_ellipsoid helps in finding the mean, minimum, and maximum
#' ellipsoids of a list of ellipsoid objects.
#'
#' @param ellipsoid_list list of ellipsoid objects created with the function
#' \code{\link{ellipsoid_fit}}.
#'
#' @details
#' Minimum and maximum ellipsoids are calculated based on the volumes of all
#' elements in \code{ellipsoid_list}.
#'
#' The mean ellipsoid is calculated by finsing the mean of all centroids and all
#' covariance matrices of all elements in \code{ellipsoid_list}.
#'
#' @export


mmm_ellipsoid <- function(ellipsoid_list) {
  # -----------
  # detecting potential errors
  if (missing(ellipsoid_list)) {
    stop("Argument ellipsoid_list is necessary to perform the analysis.")
  }

  # -----------
  # getting max, min, and preparing mean
  ## volumes
  volumes <- sapply(ellipsoid_list, function(x) {x@niche_volume})
  min_vol <- min(volumes)
  max_vol <- max(volumes)

  ## min
  min_ellipsoid <- ellipsoid_list[sapply(ellipsoid_list, function(x) {
    x@niche_volume == min_vol
  })]

  ## max
  max_ellipsoid <- ellipsoid_list[sapply(ellipsoid_list, function(x) {
    x@niche_volume == max_vol
  })]

  ## mean
  centroids <- lapply(ellipsoid_list, function(x) {x@centroid})
  covariances <- lapply(ellipsoid_list, function(x) {x@covariance_matrix})

  mean_cen <- Reduce('+', centroids) / length(centroids)
  mean_cov <- Reduce('+', covariances) / length(covariances)

  ### calculating mean ellipsoid characteristics
  level <- ellipsoid_list[[1]]@level / 100
  ndim <- length(mean_cen)
  sigma_i <- solve(mean_cov) / stats::qchisq(level, df = ndim)
  s_eigen <- eigen(sigma_i)
  s_eigenval <- s_eigen$values
  s_eigenvec <- s_eigen$vectors

  ### semi axes length for mean ellipsoid
  stds <- 1 / sqrt(s_eigenval)
  names(stds) <- letters[1:ndim]

  ### volume
  volume <- ellipsoid_volume(ndim, stds)

  ### axes coordinates for mean ellipsoid
  axes_coordinates <- lapply(1:ndim, function(x) {
    coor <- rbind((mean_cen + s_eigenvec[, x] * stds[x]),
                  (mean_cen - s_eigenvec[, x] * stds[x]))
    rownames(coor) <- paste0("vec_", 1:2)
    return(coor)
  })
  names(axes_coordinates) <- names(stds)

  ### mean ellipsoid
  mean_ellipsoid <- ellipsoid(method = ellipsoid_list[[1]]@method,
                              centroid = mean_cen,
                              covariance_matrix = mean_cov,
                              level = level * 100,
                              niche_volume = volume,
                              semi_axes_length = stds,
                              axes_coordinates = axes_coordinates)

  # -----------
  # returning results
  return(list(mean_ellipsoid = mean_ellipsoid, min_ellipsoid = min_ellipsoid[[1]],
              max_ellipsoid = max_ellipsoid[[1]]))
}
