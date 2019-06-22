#setMethod("predict", signature(object = "Raster", model = "ellipsoid"),
#          function(object, model, filename="", fun=predict, ext=NULL, const=NULL,
#                   index=1, na.rm=TRUE, inf.rm=FALSE, factors=NULL, format,
#                   datatype, overwrite=FALSE, progress='', ...) {})


#' Predict suitability derived from ellipsoid envelope models
#'
#' @description predict_ncdsuit predicts of suitability vlues based on
#' mahalanobis distances based on a centroid and a covariance matrix.
#'
#' @param ellipsoid an ellipsoid* object. If defined, arguments \code{centroid}
#' and \code{covariance_matrix} are not required.
#' @param centroid (optional) centroid to wich distance will be measured in
#' \code{projection_layers}. Length must correspond with the number of layers
#' in \code{projection_layers}. Ignored if \code{ellipsoid} is defined.
#' @param covariance_matrix (optional) covariance matrix to be used in measuring
#' the mahalanobis distance to each point in \code{projection_layers}. Ignored
#' if \code{ellipsoid} is defined.
#' @param projection_layers RasterStack of variables representing environmental
#' conditions of the scenario to which the \code{ellipsoid} model will be projected.
#' @param level (numeric) the confidence level of a pairwise confidence region
#' for the ellipsoid, expresed as percentage. Default = 95.
#' @param tolerance the tolerance for detecting linear dependencies. Default = 1e-60.

predict_ncdsuit <- function(ellipsoid, centroid = NULL, covariance_matrix = NULL,
                            projection_layers, level = 95, tolerance = 1e-60) {
  # -----------
  # detecting potential errors and preparing the data
  if (!missing(ellipsoid)) {
    centroid <- ellipsoid@centroid
    covariance_matrix <- ellipsoid@covariance_matrix
  }else {
    if (missing(centroid) | missing(covariance_matrix)) {
      stop("If ellipsoid is missing, centroid and covariance_matrix must be defined.")
    }
  }

  # raster data
  back <- na.omit(raster::values(variables))

  # -----------
  # analyses
  ## Mahalanobis distance
  maha <- mahalanobis(x = back, center = centroid, cov = covariance_matrix,
                      tol = tolerance)

  ## a point is inside the confidence region (1-alpha=confidence%) if the distance
  ## divided by the quantile of a Chi-square variable with k d.f. is less than 1
  alpha <- level / 100
  chi_sq <- qchisq(alpha, ncol(back))

  ## distances to suitabilities considering a multivariate normal distribution
  suitability <- exp(-0.5 * maha)
  suitability <- ifelse(maha / chi_sq <= 1, suitability, 0) # inside only

  # -----------
  # preparing further results
  ## prevalences
  p_no_suit_g <- length(back[suitability != 0, 1]) / length(back[, 1])
  p_no_suit_e <- length(unique(back[suitability != 0, ])[, 1]) / length(back[, 1])

  not_suitable <- data.frame(c("Geographic_space", "Environmental_space"),
                             c(p_no_suit_g, p_no_suit_e))
  names(not_suitable) <- c("Space", "Proportion_not_suitable")

  # raster generation
  suit_layer <- variables[[1]]
  suit_layer[!is.na(raster::values(suit_layer))] <- suitability

  results <- list(centroid = centroid, covariance_matrix = covariance_matrix,
                  non_suitable_area = not_suitable, suitability_layer = suit_layer)

  return(results)
}
