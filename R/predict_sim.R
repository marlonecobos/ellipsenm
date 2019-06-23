#' Predict suitability derived from simple ellipsoid envelope models
#'
#' @description predict_sim predicts suitability values based on mahalanobis
#' distances based on a centroid and a covariance matrix.
#'
#' @param object a fitted object of class ellipsoid or ellipsoid_model_sim.
#' @param projection_layers RasterStack or matrix of variables representing
#' environmental conditions of the scenario to which \code{object} will be
#' projected. See details.
#' @param prediction (character) type of prediction to be made, options are:
#' "suitability", "mahalanobis", and "both". Default = "suitability".
#' @param return_numeric (logical) whether or not to return values mahalanobis
#' distance and suitability as part of the results (it depends on the type of
#' \code{prediction} selected). Default = FALSE.
#' @param tolerance the tolerance for detecting linear dependencies.
#' Default = 1e-60.
#' @param name (character) optional, name of the file to be writen. Must include
#' format extension (e.g., ".tif"). When defined, raster predictions are not
#' returned to the environment. Default = NULL. See detals.
#' @param format (charater) if \code{name} is defined, file type to be written.
#' Must correspond with format extension in \code{name}.
#' See \code{\link[raster]{writeFormats}}.
#' @param overwrite (logical) if \code{name} is defined, whether or not to
#' overwrite an exitent file with the exact same name. Default = FALSE.
#'
#' @return
#' An ellipsoid_model_sim with new predictions.
#'
#' @details
#' Argument \code{object} must be of one of the following classes: "ellipsoid"
#' or "ellipsoid_model_sim". The prefix "suitability" or "mahalanobis" will be
#' added to \code{name} depending on the type of prediction defined in
#' \code{prediction}.
#'
#' For \code{projection_layers} variables can be given either as a RasterStack
#' or as a matrix. If a matrix is given each column represents a variable and
#' predictions are returned only as numeric vectors. In both cases, variable
#' names must match exactly the order and name of variables used to create
#' \code{object}.
#'
#' @export
#'

predict_sim <- function(object, projection_layers, prediction = "suitability",
                        return_numeric = FALSE, tolerance = 1e-60,
                        name = NULL, format, overwrite = FALSE) {
  # -----------
  # detecting potential errors
  if (!missing(object)) {
    if (!class(object)[1] %in%
        c("ellipsoid", "ellipsoid_model_sim", "ellipsoid_model_rep")) {
      stop("Argument object must be of class ellipsoid*.")
    }
  }else {
    stop("Argument object is necessary to perform the analysis.")
  }
  if (!is.null(name)) {
    if (missing(format)) {
      stop("Argument format needs to be defined if argument name is given.")
    }
  }
  if (!class(projection_layers)[1] %in% c("RasterStack", "matrix", "data.frame")) {
    stop("Argument projection_layers needs to be either a RasterStack or a matrix.")
  }
  if (!is.null(name) & class(projection_layers)[1] != "RasterStack") {
    warning("Argument projection_layers is matrix, no raster predictions will be written.")
  }

  ## solving potential problems with name
  name <- gsub("\\\\", "/", name)
  name <- unlist(strsplit(name, "/"))
  ndir <- paste0(paste(name[-length(name)], collapse = "/"), "/")
  name <- name[length(name)]

  # ellipsoid data
  centroid <- object@centroid
  covariance_matrix <- object@covariance_matrix
  level <- object@level

  # raster data
  if (class(projection_layers)[1] == "RasterStack") {
    back <- na.omit(raster::values(projection_layers))
  } else {
    back <- as.matrix(projection_layers)
  }
  db <- !duplicated.matrix(back)

  # -----------
  # analyses and preparation of results
  ## Mahalanobis distance
  maha <- mahalanobis(x = back, center = centroid, cov = covariance_matrix,
                      tol = tolerance)

  # -----------
  # conditioned results
  if (prediction == "suitability" | prediction == "mahalanobis" | prediction == "both") {
    if (prediction != "mahalanobis") {
      ## a point is inside the confidence region (1-alpha=confidence%) if the distance
      ## divided by the quantile of a Chi-square variable with k d.f. is less than 1
      alpha <- level / 100
      chi_sq <- qchisq(alpha, ncol(back))

      ## distances to suitabilities considering a multivariate normal distribution
      suitability <- exp(-0.5 * maha)
      suitability <- ifelse(maha / chi_sq <= 1, suitability, 0) # inside only

      ## prevalence
      p_suit_g <- sum(suitability != 0) / length(suitability)
      u_suit <- suitability[db]
      p_suit_e <- sum(u_suit != 0) / length(u_suit)

      prevalence <- c(prevalence_E_space = p_suit_e, prevalence_G_space = p_suit_g)

      ## preparing RasterLayers suit
      if (class(projection_layers)[1] == "RasterStack") {
        suit_layer <- projection_layers[[1]]
        suit_layer[!is.na(suit_layer[])] <- suitability
      } else {
        suit_layer <- vector()
      }

      if (prediction == "both") {
        ## preparing RasterLayers maha
        if (class(projection_layers)[1] == "RasterStack") {
          maha_layer <- projection_layers[[1]]
          maha_layer[!is.na(maha_layer[])] <- maha
        } else {
          maha_layer <- vector()
        }

        ## returning results for both type of predictions
        if (return_numeric == TRUE) {
          results <- new("ellipsoid_model_sim", method = object@method,
                                         centroid = centroid,
                                         covariance_matrix = covariance_matrix,
                                         level = level,
                                         niche_volume = object@niche_volume,
                                         semi_axes_length = object@semi_axes_length,
                                         axes_coordinates = object@axes_coordinates,
                                         mahalanobis = maha,
                                         suitability = suitability,
                                         prevalence = prevalence)
          slot(results, "prediction_maha", check = FALSE) <- maha_layer
          slot(results, "prediction_suit", check = FALSE) <- suit_layer
        } else {
          results <- new("ellipsoid_model_sim", method = object@method,
                                         centroid = centroid,
                                         covariance_matrix = covariance_matrix,
                                         level = level,
                                         niche_volume = object@niche_volume,
                                         semi_axes_length = object@semi_axes_length,
                                         axes_coordinates = object@axes_coordinates,
                                         prevalence = prevalence)
          slot(results, "prediction_maha", check = FALSE) <- maha_layer
          slot(results, "prediction_suit", check = FALSE) <- suit_layer
        }
      } else {
        ## returning results for suitability type of predictions
        if (return_numeric == TRUE) {
          results <- new("ellipsoid_model_sim", method = object@method,
                                         centroid = centroid,
                                         covariance_matrix = covariance_matrix,
                                         level = level,
                                         niche_volume = object@niche_volume,
                                         semi_axes_length = object@semi_axes_length,
                                         axes_coordinates = object@axes_coordinates,
                                         suitability = suitability,
                                         prevalence = prevalence)
          slot(results, "prediction_suit", check = FALSE) <- suit_layer
        } else {
          results <- new("ellipsoid_model_sim", method = object@method,
                                         centroid = centroid,
                                         covariance_matrix = covariance_matrix,
                                         level = level,
                                         niche_volume = object@niche_volume,
                                         semi_axes_length = object@semi_axes_length,
                                         axes_coordinates = object@axes_coordinates,
                                         prevalence = prevalence)
          slot(results, "prediction_suit", check = FALSE) <- suit_layer
        }
      }

    } else {
      ## preparing RasterLayers maha
      if (class(projection_layers)[1] == "RasterStack") {
        maha_layer <- projection_layers[[1]]
        maha_layer[!is.na(maha_layer[])] <- maha
      } else {
        maha_layer <- vector()
      }

      ## returning results for mahalanobis type of predictions
      if (return_numeric == TRUE) {
        results <- new("ellipsoid_model_sim", method = object@method,
                                       centroid = centroid,
                                       covariance_matrix = covariance_matrix,
                                       level = level,
                                       niche_volume = object@niche_volume,
                                       semi_axes_length = object@semi_axes_length,
                                       axes_coordinates = object@axes_coordinates,
                                       mahalanobis = maha)
        slot(results, "prediction_maha", check = FALSE) <- maha_layer
      } else {
        results <- new("ellipsoid_model_sim", method = object@method,
                                       centroid = centroid,
                                       covariance_matrix = covariance_matrix,
                                       level = level,
                                       niche_volume = object@niche_volume,
                                       semi_axes_length = object@semi_axes_length,
                                       axes_coordinates = object@axes_coordinates)
        slot(results, "prediction_maha", check = FALSE) <- maha_layer
      }
    }
  } else {
    stop("Argument prediction is not valid, see function's help.")
  }

  if (!is.null(name) & class(projection_layers)[1] == "RasterStack") {
    ## excluding raster predictions if name exist and layers were rasters
    if (prediction != "mahalanobis") {
      if (prediction == "both") {
        ## writing raster layers
        mname <- paste0(ndir, "mahalanobis_", name)
        sname <- paste0(ndir, "suitability_", name)
        raster::writeRaster(maha_layer, filename = mname, format = format, overwrite = overwrite)
        raster::writeRaster(suit_layer, filename = sname, format = format, overwrite = overwrite)

        ## erasing slots
        slot(results, "prediction_maha", check = FALSE) <- vector()
        slot(results, "prediction_suit", check = FALSE) <- vector()
      } else {
        sname <- paste0(ndir, "suitability_", name)
        raster::writeRaster(suit_layer, filename = sname, format = format, overwrite = overwrite)

        slot(results, "prediction_suit", check = FALSE) <- vector()
      }
    } else {
      mname <- paste0(ndir, "mahalanobis_", name)
      raster::writeRaster(maha_layer, filename = mname, format = format, overwrite = overwrite)

      slot(results, "prediction_maha", check = FALSE) <- vector()
    }
  }

  # -----------
  # returning results
  return(results)
}
