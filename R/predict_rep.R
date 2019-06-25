#' Predict suitability derived from replicated ellipsoid envelope models
#'
#' @description predicts suitability values based on mahalanobis distances
#' based on a centroid and a covariance matrix.
#'
#' @param object a fitted object of class ellipsoid_model_rep.
#' @param projection_layers RasterStack or matrix of variables representing
#' environmental conditions of the scenario to which \code{object} will be
#' projected. See details.
#' @param prediction (character) type of prediction to be made, options are:
#' "suitability", "mahalanobis", and "both". Default = "suitability".
#' @param return_numeric (logical) whether or not to return values of mahalanobis
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
#' An ellipsoid_model_rep with new predictions.
#'
#' @details
#' Predictions for all replicates will be performed. If \code{name} is defined,
#' a prefix starting in 1 will be added to each replicate. After replicate number
#' the word "suitability" or "mahalanobis" will be added depending on the type
#' of prediction defined in \code{prediction}.
#'
#' For \code{projection_layers} variables can be given either as a RasterStack
#' or as a matrix. If a matrix is given each column represents a variable and
#' predictions are returned only as numeric vectors. In both cases, variable
#' names must match exactly the order and name of variables used to create
#' \code{object}.
#'
#' @export

setMethod("predict", signature(object = "ellipsoid_model_rep"),
          function(object, projection_layers, prediction = "suitability",
                   return_numeric = FALSE, tolerance = 1e-60,
                   name = NULL, format, overwrite = FALSE) {
            # -----------
            # detecting potential errors
            if (!missing(object)) {
              if (!class(object)[1] %in% c("ellipsoid_model_rep")) {
                stop("Argument object must be of class ellipsoid* that\ncontains centroid and coariance_matrix.")
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

            # ellipsoids data
            centroid <- sapply(object@ellipsoids, function(x) {x@centroid})
            covariance_matrix <- lapply(object@ellipsoids, function(x) {x@covariance_matrix})
            level <- sapply(object@ellipsoids, function(x) {x@level})
            n_ell <- ncol(centroid)
            nam_ell <- names(object@ellipsoids)

            # raster data
            if (class(projection_layers)[1] == "RasterStack") {
              back <- na.omit(raster::values(projection_layers))
            } else {
              back <- as.matrix(projection_layers)
            }
            db <- !duplicated.matrix(back)

            # -----------
            # analyses and preparation of results
            ## layers to be filled
            if (prediction == "suitability" | prediction == "mahalanobis" | prediction == "both") {
              if (class(projection_layers)[1] == "RasterStack") {
                if (prediction != "mahalanobis") {
                  suit_layer <- projection_layers[[1]]
                  if (prediction == "both") {
                    maha_layer <- projection_layers[[1]]
                  }
                } else {
                  maha_layer <- projection_layers[[1]]
                }
              }
            } else {
              stop("Argument prediction is not valid, see function's help.")
            }

            # -----------
            # conditioned results
            if (prediction != "mahalanobis") {
              alpha <- level / 100
              chi_sq <- qchisq(alpha, ncol(back))

              maha_suit <- lapply(1:n_ell, function(x) {
                cat("\tPrediction", x, "of", n_ell, "\n")
                mah <-  mahalanobis(x = back, center = centroid[, x],
                                    cov = covariance_matrix[[x]], tol = tolerance)
                if (prediction == "both") {
                  if (class(projection_layers)[1] == "RasterStack") {
                    maha_layer[!is.na(maha_layer[])] <- mah
                  } else {
                    maha_layer <- vector()
                  }
                }

                suitability <- exp(-0.5 * mah)
                suitability <- ifelse(mah / chi_sq[x] <= 1, suitability, 0)

                if (class(projection_layers)[1] == "RasterStack") {
                  suit_layer[!is.na(suit_layer[])] <- suitability
                } else {
                  suit_layer <- vector()
                }

                p_suit_g <- sum(suitability != 0) / length(suitability)
                u_suit <- suitability[db]
                p_suit_e <- sum(u_suit != 0) / length(u_suit)
                prevalence <- c(prevalence_E_space = p_suit_e, prevalence_G_space = p_suit_g)

                if (!is.null(name)) {
                  if (class(projection_layers)[1] == "RasterStack") {
                    if (prediction == "both") {
                      mname <- paste0(ndir, x, "_mahalanobis_", name)
                      sname <- paste0(ndir, x, "_suitability_", name)
                      raster::writeRaster(maha_layer, filename = mname, format = format, overwrite = overwrite)
                      raster::writeRaster(suit_layer, filename = sname, format = format)
                    } else {
                      sname <- paste0(ndir, x, "_suitability_", name)
                      raster::writeRaster(suit_layer, filename = sname, format = format, overwrite = overwrite)
                    }
                  }

                  if (return_numeric == TRUE) {
                    if (prediction == "both") {
                      return(list(mah, suitability, prevalence))
                    } else {
                      return(list(suitability, prevalence))
                    }
                  } else {
                    return(list(prevalence))
                  }

                } else {
                  if (return_numeric == TRUE) {
                    if (prediction == "both") {
                      return(list(mah, suitability, maha_layer, suit_layer, prevalence))
                    } else {
                      return(list(suitability, suit_layer, prevalence))
                    }
                  } else {
                    if (prediction == "both") {
                      return(list(maha_layer, suit_layer, prevalence))
                    } else {
                      return(list(suit_layer, prevalence))
                    }
                  }
                }
              })

              if (return_numeric == TRUE) {
                ids <- ifelse(prediction == "both", 2, 1)
                suitability <- do.call(cbind, lapply(maha_suit, function(x) {x[[ids]]}))
                colnames(suitability) <- nam_ell
              }

              idp <- length(maha_suit[[1]])
              prevalence <- do.call(cbind, lapply(maha_suit, function(x) {x[[idp]]}))
              colnames(prevalence) <- nam_ell

              if (is.null(name)) {
                if (class(projection_layers)[1] == "RasterStack") {
                  if (prediction == "both") {
                    inds <- ifelse(return_numeric == TRUE, 4, 2)
                  } else {
                    inds <- ifelse(return_numeric == TRUE, 2, 1)
                  }
                  suit_layers <- do.call(raster::stack, lapply(maha_suit, function(x) {x[[inds]]}))
                  names(suit_layers) <- nam_ell
                } else {
                  suit_layers <- vector()
                }
              } else {
                suit_layers <- vector()
              }

              if (prediction == "both") {
                if (is.null(name)) {
                  if (class(projection_layers)[1] == "RasterStack") {
                    indm <- ifelse(return_numeric == TRUE, 3, 1)
                    maha_layers <- do.call(raster::stack, lapply(maha_suit, function(x) {x[[indm]]}))
                    names(maha_layers) <- nam_ell
                  } else {
                    maha_layers <- vector()
                  }
                } else {
                  maha_layers <- vector()
                }

                ## returning results for both type of predictions
                if (return_numeric == TRUE) {
                  maha <- do.call(cbind, lapply(maha_suit, function(x) {x[[1]]}))
                  colnames(maha) <- nam_ell

                  results <- ellipsoid_model_rep(ellipsoids = object@ellipsoids,
                                                 mahalanobis = maha,
                                                 suitability = suitability,
                                                 prevalence = prevalence)
                  slot(results, "prediction_maha", check = FALSE) <- maha_layers
                  slot(results, "prediction_suit", check = FALSE) <- suit_layers
                } else {
                  results <- ellipsoid_model_rep(ellipsoids = object@ellipsoids,
                                                 prevalence = prevalence)
                  slot(results, "prediction_maha", check = FALSE) <- maha_layers
                  slot(results, "prediction_suit", check = FALSE) <- suit_layers
                }
              } else {
                ## returning results for suitability type of predictions
                if (return_numeric == TRUE) {
                  results <- ellipsoid_model_rep(ellipsoids = object@ellipsoids,
                                                 suitability = suitability,
                                                 prevalence = prevalence)
                  slot(results, "prediction_suit", check = FALSE) <- suit_layers
                } else {
                  results <- ellipsoid_model_rep(ellipsoids = object@ellipsoids,
                                                 prevalence = prevalence)
                  slot(results, "prediction_suit", check = FALSE) <- suit_layers
                }
              }

            } else {
              maha <- lapply(1:n_ell, function(x) {
                cat("\tPrediction", x, "of", n_ell, "\n")
                mah <-  mahalanobis(x = back, center = centroid[, x],
                                    cov = covariance_matrix[[x]], tol = tolerance)

                if (class(projection_layers)[1] == "RasterStack") {
                  maha_layer[!is.na(maha_layer[])] <- mah
                } else {
                  maha_layer <- vector()
                }

                if (!is.null(name)) {
                  if (class(projection_layers)[1] == "RasterStack") {
                    mname <- paste0(ndir, x, "_mahalanobis_", name)
                    raster::writeRaster(maha_layer, filename = mname, format = format, overwrite = overwrite)
                  }

                  if (return_numeric == TRUE) {
                    return(list(mah))
                  } else {
                    return(list())
                  }

                } else {
                  if (return_numeric == TRUE) {
                    return(list(mah, maha_layer))
                  } else {
                    return(list(maha_layer))
                  }
                }
              })

              if (is.null(name)) {
                if (class(projection_layers)[1] == "RasterStack") {
                  indm <- ifelse(return_numeric == TRUE, 2, 1)
                  maha_layers <- do.call(raster::stack, lapply(maha, function(x) {x[[indm]]}))
                  names(maha_layers) <- nam_ell
                } else {
                  maha_layers <- vector()
                }
              } else {
                maha_layers <- vector()
              }

              ## returning results for mahalanobis type of predictions
              if (return_numeric == TRUE) {
                maha <- do.call(cbind, lapply(maha, function(x) {x[[1]]}))
                colnames(maha) <- nam_ell

                results <- ellipsoid_model_rep(ellipsoids = object@ellipsoids,
                                               mahalanobis = maha)
                slot(results, "prediction_maha", check = FALSE) <- maha_layers
              } else {
                results <- ellipsoid_model_rep(ellipsoids = object@ellipsoids)
                slot(results, "prediction_maha", check = FALSE) <- maha_layers
              }
            }

            # -----------
            # returning results
            return(results)
          }
)
