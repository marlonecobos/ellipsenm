#' Predict suitability derived from simple ellipsoid envelope models
#'
#' @description predicts suitability values based on mahalanobis distances
#' based on a centroid and a covariance matrix.
#'
#' @param object a fitted object of class ellipsoid or ellipsoid_model_sim.
#' @param projection_variables RasterStack or matrix of variables representing
#' environmental conditions of the scenario to which \code{object} will be
#' projected. See details.
#' @param prediction (character) type of prediction to be made, options are:
#' "suitability", "mahalanobis", and "both". Default = "suitability".
#' @param truncate (logical) whether or not to truncate values of suitability
#' based on ellipsoid limits. All values outside the ellipsoid will be zero.
#' Default = TRUE.
#' @param return_numeric (logical) whether or not to return values of mahalanobis
#' distance and suitability as part of the results (it depends on the type of
#' \code{prediction} selected). If \code{projection_variables} is a RasterStack,
#' default = FALSE, but can be changed to TRUE; if \code{projection_variables} is
#' a matrix, default = TRUE and cannot be changed. See details.
#' @param tolerance the tolerance for detecting linear dependencies.
#' Default = 1e-60.
#' @param name (character) optional, a name for the files to be writen. When
#' defined, raster predictions and numeric results are not returned as part of
#' the ellipsoid* object unless \code{force_return} = TRUE. File extensions will
#' be added as needed for writing raster and numeric results. Default = NULL.
#' See details.
#' @param format (charater) if \code{name} is defined, raster type to be written.
#' See \code{\link[raster]{writeFormats}} for details and options.
#' @param overwrite (logical) if \code{name} is defined, whether or not to
#' overwrite an exitent file with the exact same name. Default = FALSE.
#' @param force_return (logical) whether or not to force returning numeric and
#' raster results as part of the ellipsoid* object when \code{name} is defined.
#'
#' @return
#' An ellipsoid_model_sim with new predictions. If \code{name} is defined, csv
#' files with numeric results and raster files with the geographic predictions
#' will be written.
#'
#' @details
#' Argument \code{object} must be of one of the following classes: "ellipsoid"
#' or "ellipsoid_model_sim". The prefix "suitability" or "mahalanobis" will be
#' added to \code{name} depending on the type of prediction defined in
#' \code{prediction}. File type (extention) will be added to \code{name}, if
#' defined, .csv for numeric results and any of the ones described in
#' \code{\link[raster]{writeFormats}} depending on \code{format}.
#'
#' Argument \code{projection_variables} variables can be defined either as a
#' RasterStack or as a matrix. If a matrix is given each column represents a
#' variable and predictions are returned only as numeric vectors. In both cases,
#' variable names must match exactly the order and name of variables used to
#' create \code{object}.
#'
#' The only scenarios in which none of the numeric results will be returned are:
#' if \code{projection_variables} is a RasterStack and \code{return numeric} is
#' set as FALSE, and if \code{name} is defined and \code{force_return} is set as
#' FALSE, even if \code{return numeric} = TRUE.
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
#' # fitting a minimum volume ellipsoid
#' ellips1 <- ellipsoid_fit(data = occurrences, longitude = "longitude",
#'                          latitude = "latitude", method = "mve1",
#'                          level = 99, raster_layers = vars)
#'
#' # predicting suitability (some slots will be empty if not required)
#' prediction <- predict(object = ellips1, projection_variables = vars,
#'                       prediction = "suitability")
#'
#' class(prediction)
#'
#' # predicting mahalanobis distance
#' prediction1 <- predict(object = ellips1, projection_variables = vars,
#'                        prediction = "mahalanobis")
#'
#'
#' # predicting both things
#' prediction2 <- predict(object = ellips1, projection_variables = vars,
#'                        prediction = "both")

setMethod("predict", signature(object = "ellipsoid"),
          function(object, projection_variables, prediction = "suitability",
                   truncate = TRUE, return_numeric, tolerance = 1e-60, name = NULL,
                   format, overwrite = FALSE, force_return = FALSE) {
            # -----------
            # detecting potential errors
            if (!missing(object)) {
              if (!class(object)[1] %in%
                  c("ellipsoid", "ellipsoid_model_sim")) {
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
            if (!class(projection_variables)[1] %in% c("RasterStack", "RasterBrick", "matrix", "data.frame")) {
              stop("Argument projection_variables needs to be either a RasterStack or a matrix.")
            } else {
              if (class(projection_variables)[1] == "RasterBrick") {
                projection_variables <- raster::stack(projection_variables)
              }
            }
            if (!is.null(name)) {
              if (class(projection_variables)[1] != "RasterStack") {
                message("Argument projection_variables is a matrix, no raster predictions will be returned.")
              }
            } else {
              force_return <- FALSE
            }

            ## solving potential problems with name
            if (!is.null(name)) {
              name <- gsub("\\\\", "/", name)
              name <- unlist(strsplit(name, "/"))
              ndir <- paste0(paste(name[-length(name)], collapse = "/"), "/")
              name <- name[length(name)]
              num_format <- ".csv"
              ras_format <- rformat_type(format)
            }

            # ellipsoid data
            centroid <- object@centroid
            covariance_matrix <- object@covariance_matrix
            level <- object@level

            # raster data
            if (class(projection_variables)[1] == "RasterStack") {
              return_numeric <- ifelse(missing(return_numeric), FALSE, return_numeric)
              back <- na.omit(raster::values(projection_variables))
            } else {
              return_numeric <- TRUE
              back <- as.matrix(projection_variables)
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
                suitab <- exp(-0.5 * maha)
                suitability <- ifelse(maha / chi_sq <= 1, suitab, 0) # inside only

                ## prevalence
                p_suit_g <- sum(suitability != 0) / length(suitability)
                u_suit <- suitability[db]
                p_suit_e <- sum(u_suit != 0) / length(u_suit)

                prevalence <- c(prevalence_E_space = p_suit_e, prevalence_G_space = p_suit_g)

                ## redefining suitability by trucate option
                if (truncate == FALSE) {suitability <- suitab}

                ## preparing RasterLayers suit
                if (class(projection_variables)[1] == "RasterStack") {
                  suit_layer <- projection_variables[[1]]
                  suit_layer[!is.na(suit_layer[])] <- suitability
                  names(suit_layer) <- "suitability"
                } else {
                  suit_layer <- vector()
                }

                if (prediction == "both") {
                  ## preparing RasterLayers maha
                  if (class(projection_variables)[1] == "RasterStack") {
                    maha_layer <- projection_variables[[1]]
                    maha_layer[!is.na(maha_layer[])] <- maha
                    names(maha_layer) <- "mahalanobis"
                  } else {
                    maha_layer <- vector()
                  }

                  ## returning results for both type of predictions
                  if (return_numeric == TRUE) {
                    results <- ellipsoid_model_sim(method = object@method,
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
                    results <- ellipsoid_model_sim(method = object@method,
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
                    results <- ellipsoid_model_sim(method = object@method,
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
                    results <- ellipsoid_model_sim(method = object@method,
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
                if (class(projection_variables)[1] == "RasterStack") {
                  maha_layer <- projection_variables[[1]]
                  maha_layer[!is.na(maha_layer[])] <- maha
                  names(maha_layer) <- "mahalanobis"
                } else {
                  maha_layer <- vector()
                }

                ## returning results for mahalanobis type of predictions
                if (return_numeric == TRUE) {
                  results <- ellipsoid_model_sim(method = object@method,
                                                 centroid = centroid,
                                                 covariance_matrix = covariance_matrix,
                                                 level = level,
                                                 niche_volume = object@niche_volume,
                                                 semi_axes_length = object@semi_axes_length,
                                                 axes_coordinates = object@axes_coordinates,
                                                 mahalanobis = maha)
                  slot(results, "prediction_maha", check = FALSE) <- maha_layer
                } else {
                  results <- ellipsoid_model_sim(method = object@method,
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

            if (!is.null(name)) {
              ## excluding raster predictions if name exist and layers were rasters
              if (prediction != "mahalanobis") {
                if (prediction == "both") {
                  mnamec <- paste0(ndir, "mahalanobis_", name, num_format)
                  snamec <- paste0(ndir, "suitability_", name, num_format)
                  suppressMessages(data.table::fwrite(data.frame(mahalanobis = maha),
                                                      file = mnamec))
                  suppressMessages(data.table::fwrite(data.frame(suitability = suitability),
                                                      file = snamec))
                  if (force_return == FALSE) {
                    slot(results, "mahalanobis", check = FALSE) <- vector()
                    slot(results, "suitability", check = FALSE) <- vector()
                  }
                } else {
                  snamec <- paste0(ndir, "suitability_", name, num_format)
                  suppressMessages(data.table::fwrite(data.frame(suitability = suitability),
                                                      file = snamec))
                  if (force_return == FALSE) {
                    slot(results, "suitability", check = FALSE) <- vector()
                  }
                }
              } else {
                mnamec <- paste0(ndir, "mahalanobis_", name, num_format)
                suppressMessages(data.table::fwrite(data.frame(mahalanobis = maha),
                                                    file = mnamec))
                if (force_return == FALSE) {
                  slot(results, "mahalanobis", check = FALSE) <- vector()
                }
              }

              if (class(projection_variables)[1] == "RasterStack") {
                if (prediction != "mahalanobis") {
                  if (prediction == "both") {
                    ## writing raster layers
                    mname <- paste0(ndir, "mahalanobis_", name, ras_format)
                    sname <- paste0(ndir, "suitability_", name, ras_format)
                    raster::writeRaster(maha_layer, filename = mname, format = format,
                                        overwrite = overwrite)
                    raster::writeRaster(suit_layer, filename = sname, format = format,
                                        overwrite = overwrite)

                    ## erasing slots
                    if (force_return == FALSE) {
                      slot(results, "prediction_maha", check = FALSE) <- vector()
                      slot(results, "prediction_suit", check = FALSE) <- vector()
                    }
                  } else {
                    sname <- paste0(ndir, "suitability_", name, ras_format)
                    raster::writeRaster(suit_layer, filename = sname, format = format,
                                        overwrite = overwrite)

                    if (force_return == FALSE) {
                      slot(results, "prediction_suit", check = FALSE) <- vector()
                    }
                  }
                } else {
                  mname <- paste0(ndir, "mahalanobis_", name, ras_format)
                  raster::writeRaster(maha_layer, filename = mname, format = format,
                                      overwrite = overwrite)

                  if (force_return == FALSE) {
                    slot(results, "prediction_maha", check = FALSE) <- vector()
                  }
                }
              }
            }

            # -----------
            # returning results
            return(results)
          }
)
