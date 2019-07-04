#' Predict suitability derived from replicated ellipsoid envelope models
#'
#' @description predicts suitability values based on mahalanobis distances
#' based on a centroid and a covariance matrix.
#'
#' @param object a fitted object of class ellipsoid_model_rep.
#' @param projection_variables RasterStack or matrix of variables representing
#' environmental conditions of the scenario to which \code{object} will be
#' projected. See details.
#' @param prediction (character) type of prediction to be made, options are:
#' "suitability", "mahalanobis", and "both". Default = "suitability".
#' @param return_numeric (logical) whether or not to return values of mahalanobis
#' distance and suitability as part of the results (it depends on the type of
#' \code{prediction} selected). If \code{projection_variables} is a RasterStack,
#' default = FALSE; if \code{projection_variables} is a matrix, default = TRUE. For
#' both options the default can be changed. See details.
#' @param tolerance the tolerance for detecting linear dependencies.
#' Default = 1e-60.
#' @param name (character) optional, a name for the files to be writen. When
#' defined, raster predictions and numeric results (with exceptions) are not
#' returned as part of the ellipsoid_model_rep object unless
#' \code{force_return} = TRUE. File extensions will be added as needed for
#' writing raster and numeric results. Default = NULL. See details.
#' @param format (charater) if \code{name} is defined, raster type to be written.
#' See \code{\link[raster]{writeFormats}} for details and options.
#' @param overwrite (logical) if \code{name} is defined, whether or not to
#' overwrite an exitent file with the exact same name. Default = FALSE.
#' @param force_return (logical) whether or not to force returning numeric and
#' raster results for one of the ellipsoids defined by name in \code{return_name},
#' as part of the ellipsoid_model_rep object when \code{name} is defined. See
#' details.
#' @param return_name (character) names of the ellipsoid (part of \code{object})
#' for which numeric and raster results will be forced to return if \code{name}
#' is defined. Default = NULL.
#'
#' @return
#' An ellipsoid_model_rep with new predictions. If \code{name} is defined, csv
#' files with numeric results and raster files with the geographic predictions
#' will be written.
#'
#' @details
#' Predictions for all replicates will be performed. If \code{name} is defined,
#' a prefix starting in 1 will be added to each replicate. After replicate number
#' the word "suitability" or "mahalanobis" will be added depending on the type
#' of prediction defined in \code{prediction}. File type (extention) will be
#' added to \code{name}, if defined, .csv for numeric results and any of the ones
#' described in \code{\link[raster]{writeFormats}} depending on \code{format}.
#'
#' Argument \code{projection_variables} variables can be defined either as a
#' RasterStack or as a matrix. If a matrix is given each column represents a
#' variable and predictions are returned only as numeric matrices. In both cases,
#' variable names must match exactly the order and name of variables used to
#' create \code{object}.
#'
#' If \code{projection_variables} is a matrix at least one numeric result will be
#' returned even if \code{return numeric} is set as FALSE; if \code{return_name}
#' is defined this indicates the ellipsoid for which the numeric result will
#' return, if not defined, the results for the first ellipsoid will return.
#'
#' The only scenarios in which none of the numeric results will be returned are:
#' if \code{projection_variables} is a RasterStack and \code{return numeric} is
#' set as FALSE, and if \code{name} is defined and \code{force_return} is set as
#' FALSE. However, for the latter, if \code{force_return} is set as TRUE, raster
#' and numeric results will be returned for the ellipsoid defined in
#' \code{return_name}.
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
#' # getting values of variables in data
#' occurrences1 <- cbind(occurrences[, 2:3], raster::extract(vars, occurrences[, 2:3]))
#'
#' # subsampling data for construction of multiple ellipsoids
#' subsamples <- data_subsample(data = occurrences1, replicates = 10,
#'                              replicate_type = "bootstrap")
#'
#' # fitting ellipsoids for the 10 subsamples
#' ellipsoids <- lapply(subsamples, function (x) {
#'   ellipsoid_fit(data = x, longitude = "longitude",
#'                 latitude = "latitude", method = "mve1",
#'                 level = 99)
#' })
#' length(ellipsoids)
#' names(ellipsoids) <- paste0("replicate_", 1:10)
#'
#' # creating an ellipsoid_model_rep object with replicates and mean, min, max
#' ell_rep <- new("ellipsoid_model_rep",
#'                ellipsoids = ellipsoids) # again some slots empty here
#'
#' # predicting suitability
#' prediction_rep <- predict(object = ell_rep, projection_variables = vars,
#'                           prediction = "suitability")
#'
#' # predicting mahalanobis distance
#' prediction_rep <- predict(object = ell_rep, projection_variables = vars,
#'                           prediction = "mahalanobis")
#'
#' # predicting both things
#' prediction_rep <- predict(object = ell_rep, projection_variables = vars,
#'                           prediction = "both")

setMethod("predict", signature(object = "ellipsoid_model_rep"),
          function(object, projection_variables, prediction = "suitability",
                   return_numeric, tolerance = 1e-60, name = NULL, format,
                   overwrite = FALSE, force_return = FALSE, return_name = NULL) {
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
            if (!class(projection_variables)[1] %in% c("RasterStack", "RasterBrick", "matrix", "data.frame")) {
              stop("Argument projection_variables needs to be either a RasterStack or a matrix.")
            } else {
              if (class(projection_variables)[1] == "RasterBrick") {
                projection_variables <- raster::stack(projection_variables)
              }
              if (class(projection_variables)[1] == "RasterStack" & missing(return_numeric)) {
                return_numeric <- FALSE
              }
              if (class(projection_variables)[1] != "RasterStack" & missing(return_numeric)){
                return_numeric <- TRUE
              }
            }
            if (!is.null(name)) {
              if (class(projection_variables)[1] != "RasterStack") {
                message("\nArgument projection_variables is a matrix, no raster predictions will be returned.\n")
              }
            } else {
              force_return <- FALSE
            }

            # -----------
            # preparing data
            ## solving potential problems with name and creting prefixes and colnames
            if (!is.null(name)) {
              name <- gsub("\\\\", "/", name)
              name <- unlist(strsplit(name, "/"))
              ndir <- paste0(paste(name[-length(name)], collapse = "/"), "/")
              name <- name[length(name)]
              num_format <- ".csv"
              ras_format <- rformat_type(format)
            }

            nam <- names(object@ellipsoids)
            if (is.null(nam)) {
              enames <- as.character(1:length(object@ellipsoids))
              nam_ell <- paste0("ellipsoid", enames)
            } else {
              if (length(grep("replicate", nam)) > 0 & length(grep("mean", nam)) > 0) {
                enames <- c(1:(length(nam) - 3), "mean", "min", "max")
                nam_ell <- nam
              } else {
                enames <- nam
                nam_ell <- nam
              }

              if (force_return == TRUE & is.null(return_name)) {
                return_name <- nam_ell[1]
              }
            }

            ## preparing variants
            if (return_numeric == FALSE & is.null(name) & is.null(return_name) &
                class(projection_variables)[1] != "RasterStack") {
              return_name <- nam_ell[1]
            }
            if (return_numeric == FALSE & is.null(name) & is.null(return_name) &
                class(projection_variables)[1] == "RasterStack") {
              return_name <- "ipmoisbsel392"
            }
            if (return_numeric == FALSE & !is.null(name) & is.null(return_name) &
                class(projection_variables)[1] == "RasterStack") {
              return_name <- "ipmoisbsel392"
            }
            if (return_numeric == TRUE & !is.null(name) & is.null(return_name) &
                class(projection_variables)[1] == "RasterStack") {
              return_name <- "ipmoisbsel392"
            }

            # ellipsoids data
            centroid <- sapply(object@ellipsoids, function(x) {x@centroid})
            covariance_matrix <- lapply(object@ellipsoids, function(x) {x@covariance_matrix})
            level <- sapply(object@ellipsoids, function(x) {x@level})
            n_ell <- ncol(centroid)

            # raster data
            if (class(projection_variables)[1] == "RasterStack") {
              back <- na.omit(raster::values(projection_variables))
            } else {
              back <- as.matrix(projection_variables)
            }
            db <- !duplicated.matrix(back)

            # -----------
            # analyses and preparation of results
            ## layers to be filled
            if (prediction == "suitability" | prediction == "mahalanobis" | prediction == "both") {
              if (class(projection_variables)[1] == "RasterStack") {
                if (prediction != "mahalanobis") {
                  suit_layer <- projection_variables[[1]]
                  if (prediction == "both") {
                    maha_layer <- projection_variables[[1]]
                  }
                } else {
                  maha_layer <- projection_variables[[1]]
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
                  if (class(projection_variables)[1] == "RasterStack") {
                    maha_layer[!is.na(maha_layer[])] <- mah
                  } else {
                    maha_layer <- vector()
                  }
                }

                suitability <- exp(-0.5 * mah)
                suitability <- ifelse(mah / chi_sq[x] <= 1, suitability, 0)

                if (class(projection_variables)[1] == "RasterStack") {
                  suit_layer[!is.na(suit_layer[])] <- suitability
                } else {
                  suit_layer <- vector()
                }

                p_suit_g <- sum(suitability != 0) / length(suitability)
                u_suit <- suitability[db]
                p_suit_e <- sum(u_suit != 0) / length(u_suit)
                prevalence <- c(prevalence_E_space = p_suit_e, prevalence_G_space = p_suit_g)

                if (!is.null(name)) {
                  if (prediction == "both") {
                    mnamec <- paste0(ndir, enames[x], "_mahalanobis_", name, num_format)
                    snamec <- paste0(ndir, enames[x], "_suitability_", name, num_format)
                    suppressMessages(data.table::fwrite(data.frame(mahalanobis = mah),
                                                        file = mnamec))
                    suppressMessages(data.table::fwrite(data.frame(suitability = suitability),
                                                        file = snamec))
                  } else {
                    snamec <- paste0(ndir, enames[x], "_suitability_", name, num_format)
                    suppressMessages(data.table::fwrite(data.frame(suitability = suitability),
                                                        file = snamec))
                  }
                  if (class(projection_variables)[1] == "RasterStack") {
                    if (prediction == "both") {
                      mname <- paste0(ndir, enames[x], "_mahalanobis_", name, ras_format)
                      sname <- paste0(ndir, enames[x], "_suitability_", name, ras_format)
                      raster::writeRaster(maha_layer, filename = mname, format = format,
                                          overwrite = overwrite)
                      raster::writeRaster(suit_layer, filename = sname, format = format,
                                          overwrite = overwrite)
                    } else {
                      sname <- paste0(ndir, enames[x], "_suitability_", name, ras_format)
                      raster::writeRaster(suit_layer, filename = sname, format = format,
                                          overwrite = overwrite)
                    }
                  }

                  if (force_return == TRUE) {
                    if (nam_ell[x] == return_name) {
                      if (return_numeric == TRUE) {
                        if (prediction == "both") {
                          return(list(maha = mah, suit = suitability, m_layer = maha_layer,
                                      s_layer = suit_layer, prev = prevalence))
                        } else {
                          return(list(suit = suitability, s_layer = suit_layer,
                                      prev = prevalence))
                        }
                      } else {
                        if (prediction == "both") {
                          return(list(m_layer = maha_layer, s_layer = suit_layer,
                                      prev = prevalence))
                        } else {
                          return(list(s_layer = suit_layer, prev = prevalence))
                        }
                      }
                    } else {
                      return(list(prev = prevalence))
                    }
                  } else {
                    return(list(prev = prevalence))
                  }
                } else {
                  if (return_numeric == TRUE) {
                    if (prediction == "both") {
                      return(list(maha = mah, suit = suitability, m_layer = maha_layer,
                                  s_layer = suit_layer, prev = prevalence))
                    } else {
                      return(list(suit = suitability, s_layer = suit_layer,
                                  prev = prevalence))
                    }
                  } else {
                    if (nam_ell[x] == return_name) {
                      if (prediction == "both") {
                        return(list(maha = mah, suit = suitability, m_layer = maha_layer,
                                    s_layer = suit_layer, prev = prevalence))
                      } else {
                        return(list(suit = suitability, s_layer = suit_layer,
                                    prev = prevalence))
                      }
                    } else {
                      if (prediction == "both") {
                        return(list(m_layer = maha_layer, s_layer = suit_layer,
                                    prev = prevalence))
                      } else {
                        return(list(s_layer = suit_layer, prev = prevalence))
                      }
                    }
                  }
                }
              })
              names(maha_suit) <- nam_ell

              prevalence <- do.call(cbind, lapply(maha_suit, function(x) {x[["prev"]]}))
              colnames(prevalence) <- nam_ell

              if (is.null(name)) {
                if (return_numeric == TRUE) {
                  suitability <- do.call(cbind, lapply(maha_suit, function(x) {x[["suit"]]}))
                  colnames(suitability) <- nam_ell
                } else {
                  suitability <- data.frame(maha_suit[[return_name]][["suit"]])
                  colnames(suitability) <- return_name
                }

                if (class(projection_variables)[1] == "RasterStack") {
                  suit_layers <- do.call(raster::stack,
                                         lapply(maha_suit, function(x) {x[["s_layer"]]}))
                  names(suit_layers) <- nam_ell
                } else {
                  suit_layers <- vector()
                }
              } else {
                if (force_return == TRUE) {
                  suitability <- data.frame(maha_suit[[return_name]][["suit"]])
                  colnames(suitability) <- return_name
                  suit_layers <- maha_suit[[return_name]][["s_layer"]]
                  names(suit_layers) <- return_name
                } else {
                  suitability <- vector()
                  suit_layers <- vector()
                }
              }

              if (prediction == "both") {
                if (is.null(name)) {
                  if (return_numeric == TRUE) {
                    maha <- do.call(cbind, lapply(maha_suit, function(x) {x[["maha"]]}))
                    colnames(maha) <- nam_ell
                  } else {
                    maha <- data.frame(maha_suit[[return_name]][["maha"]])
                    colnames(maha) <- return_name
                  }

                  if (class(projection_variables)[1] == "RasterStack") {
                    maha_layers <- do.call(raster::stack,
                                           lapply(maha_suit, function(x) {x[["m_layer"]]}))
                    names(maha_layers) <- nam_ell
                  } else {
                    maha_layers <- vector()
                  }
                } else {
                  if (force_return == TRUE) {
                    maha <- data.frame(maha_suit[[return_name]][["maha"]])
                    colnames(maha) <- return_name
                    maha_layers <- maha_suit[[return_name]][["m_layer"]]
                    names(maha_layers) <- return_name
                  } else {
                    maha <- vector()
                    maha_layers <- vector()
                  }
                }

                ## returning results for both type of predictions
                results <- ellipsoid_model_rep(ellipsoids = object@ellipsoids,
                                               prevalence = prevalence)
                slot(results, "mahalanobis", check = FALSE) <- maha
                slot(results, "suitability", check = FALSE) <- suitability
                slot(results, "prediction_maha", check = FALSE) <- maha_layers
                slot(results, "prediction_suit", check = FALSE) <- suit_layers
              } else {
                ## returning results for suitability type of predictions
                results <- ellipsoid_model_rep(ellipsoids = object@ellipsoids,
                                               prevalence = prevalence)
                slot(results, "suitability", check = FALSE) <- suitability
                slot(results, "prediction_suit", check = FALSE) <- suit_layers
              }

            } else {
              maha <- lapply(1:n_ell, function(x) {
                cat("\tPrediction", x, "of", n_ell, "\n")
                mah <-  mahalanobis(x = back, center = centroid[, x],
                                    cov = covariance_matrix[[x]], tol = tolerance)

                if (class(projection_variables)[1] == "RasterStack") {
                  maha_layer[!is.na(maha_layer[])] <- mah
                } else {
                  maha_layer <- vector()
                }

                if (!is.null(name)) {
                  mnamec <- paste0(ndir, enames[x], "_mahalanobis_", name, num_format)
                  suppressMessages(data.table::fwrite(data.frame(mahalanobis = mah),
                                                      file = mnamec))

                  if (class(projection_variables)[1] == "RasterStack") {
                    mname <- paste0(ndir, enames[x], "_mahalanobis_", name, ras_format)
                    raster::writeRaster(maha_layer, filename = mname, format = format,
                                        overwrite = overwrite)
                  }

                  if (force_return == TRUE) {
                    if (nam_ell[x] == return_name) {
                      if (return_numeric == TRUE) {
                        return(list(maha = mah, m_layer = maha_layer))
                      } else {
                        return(list(m_layer = maha_layer))
                      }
                    } else {
                      return(list())
                    }
                  } else {
                    return(list())
                  }
                } else {
                  if (return_numeric == TRUE) {
                    return(list(maha = mah, m_layer = maha_layer))
                  } else {
                    if (nam_ell[x] == return_name) {
                      return(list(maha = mah, m_layer = maha_layer))
                    } else {
                      return(list(m_layer = maha_layer))
                    }
                  }
                }
              })
              names(maha) <- nam_ell

              if (is.null(name)) {
                if (class(projection_variables)[1] == "RasterStack") {
                  maha_layers <- do.call(raster::stack,
                                         lapply(maha, function(x) {x[["m_layer"]]}))
                  names(maha_layers) <- nam_ell
                } else {
                  maha_layers <- vector()
                }

                if (return_numeric == TRUE) {
                  maha <- do.call(cbind, lapply(maha, function(x) {x[["maha"]]}))
                  colnames(maha) <- nam_ell
                } else {
                  maha <- data.frame(maha[[return_name]][["maha"]])
                  colnames(maha) <- return_name
                }
              } else {
                if (force_return == TRUE) {
                  maha_layers <- maha[[return_name]][["m_layer"]]
                  names(maha_layers) <- return_name
                  maha <- data.frame(maha[[return_name]][["maha"]])
                  colnames(maha) <- return_name
                } else {
                  maha <- vector()
                  maha_layers <- vector()
                }
              }

              ## returning results for mahalanobis type of predictions
              results <- ellipsoid_model_rep(ellipsoids = object@ellipsoids)
              slot(results, "mahalanobis", check = FALSE) <- maha
              slot(results, "prediction_maha", check = FALSE) <- maha_layers
            }

            # -----------
            # returning results
            return(results)
          }
)
