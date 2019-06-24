#' Predict suitability derived from ellipsoid envelope models
#'
#' @description prediction of suitability values based on mahalanobis
#' distances based on a centroid and a covariance matrix.
#'
#' @param object a fitted object of class ellipsoid*. See details.
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
#' An ellipsoid_model* with new predictions.
#'
#' @details
#' Argument \code{object} must be of one of the following classes: "ellipsoid",
#' "ellipsoid_model_sim", or "ellipsoid_model_rep". A prefix depending on the
#' type of prediction "suitability" or "mahalanobis" will be added to \code{name},
#' if defined.
#'
#' If \code{object} is a replicated model (i.e., "ellipsoid_model_rep"),
#' predictions for all replicates will be performed. For this class of ellipsoid,
#' if \code{name} is defined, a prefix starting in 1 will be added to each
#' replicate.
#'
#' For \code{projection_layers} variables can be given either as a RasterStack
#' or as a matrix. If a matrix is given each column represents a variable and
#' predictions are returned only as numeric vectors. In both cases, variable
#' names must match exactly the order and name of variables used to create
#' \code{object}.
#'
#' @export
#'

setMethod("predict", signature(object = "ellipsoid"),
          function(object, projection_layers, prediction = "suitability",
                   return_numeric = FALSE, tolerance = 1e-60,
                   name = NULL, format, overwrite = FALSE) {
            # -----------
            # detecting potential errors is mostly done in internal functions
            if (!missing(object)) {
              if (!class(object)[1] %in%
                  c("ellipsoid_basic", "ellipsoid_model_sim", "ellipsoid_model_rep")) {
                stop("Argument object must be one of class ellipsoid* that\ncontains centroid and coariance_matrix.")
              }
            }else {
              stop("Argument object is necessary to perform the analysis.")
            }

            # -----------
            # prediction depending on ellipsoid type
            if (class(object)[1] == "ellipsoid" | class(object)[1] == "ellipsoid_model_sim") {

              results <- predict_sim(object, projection_layers, prediction,
                                     return_numeric, tolerance,
                                     name, format, overwrite)

            } else {

              results <- predict_rep(object, projection_layers, prediction,
                                     return_numeric, tolerance,
                                     name, format, overwrite)

            }

            # -----------
            # returning results
            return(results)
          }

)
