#' Summary of attributes and results
#' @name summary
#' @aliases summary,ellipsoid-method summary,ellipsoid_model_rep-method
#' @aliases summary,calibration_ellipsoid-method summary,data_overlap-method
#' @aliases summary,overlap_ellipsoid-method
#' @param object object of class ellipsoid*, ellipsoid_model_rep,
#' calibration_ellipsoid, data_overlap, or overlap_ellipsoid.
#' @export
#' @return
#' A written summary.
#' @rdname summary
setMethod("summary", signature(object = "ellipsoid"),
          function(object) {
            # -----------
            # detecting potential errors
            if (!missing(object)) {
              clo <- class(object)[1]
              if (!clo %in% c("ellipsoid", "ellipsoid_model_sim")) {
                stop("Argument 'object' must be of class ellipsoid*.")
              }
            }else {
              stop("Argument 'object' is necessary.")
            }

            if (clo == "ellipsoid") {
              cat("\n                       Summary of ellipsoid object\n")
              cat("---------------------------------------------------------------------------\n")
              cat("\n----------------------------Ellipsoid metadata-----------------------------\n")
              cat("\nMethod:\t", object@method)
              cat("\n\nLevel:\t", object@level)
              cat("\n\nCentroid:\n", paste0(names(object@centroid), "\t",
                                            object@centroid, "\n"))
              cat("\nCovariance matrix:\n")
              print(object@covariance_matrix)
              cat("\nVolume:\t", object@niche_volume)
              #cat("\n\nSemi-axes length:\n", paste0(names(object@semi_axes_length),
              #                                      "\t", object@semi_axes_length,
              #                                      "\n"))
              #cat("\nAxes coordinates:")
              #a_cord <- object@axes_coordinates
              #cords <- lapply(1:length(a_cord), function(x) {
              #  cat("\n", letters[x], "\n"); print(a_cord[[x]])
              #})

            } else {
              cat("\n                  Summary of ellipsoid_model_sim object\n")
              cat("---------------------------------------------------------------------------\n")
              cat("\n----------------------------Ellipsoid metadata-----------------------------\n")
              cat("\nMethod:\t", object@method)
              cat("\n\nLevel:\t", object@level)
              cat("\n\nCentroid:\n", paste0(names(object@centroid), "\t",
                                            object@centroid, "\n"))
              cat("\nCovariance matrix:\n")
              print(object@covariance_matrix)
              cat("\nVolume:\t", object@niche_volume)
              #cat("\n\nSemi-axes length:\n",
              #    paste0(names(object@semi_axes_length), "\t",
              #           object@semi_axes_length, "\n"))
              #cat("\nAxes coordinates:")
              #a_cord <- object@axes_coordinates
              #cords <- lapply(1:length(a_cord), function(x) {
              #  cat("\n", letters[x], "\n"); print(a_cord[[x]])
              #})
              cat("\n\n--------------------------Summary of predictions---------------------------\n")
              if (length(object@mahalanobis) > 0) {
                cat("\nMahalanobis quantiles:\n")
                print(quantile(object@mahalanobis))
              }
              if (class(object@prediction_maha)[1] == "RasterLayer") {
                cat("\n\nRaster mahalanobis:\n")
                print(summary(object@prediction_maha))
              }
              if (length(object@suitability) > 0) {
                cat("\n\nSuitability  quantiles:\n")
                print(quantile(object@suitability))
              }
              if (class(object@prediction_suit)[1] == "RasterLayer") {
                cat("\n\nRaster suitability:\n")
                print(summary(object@prediction_suit))
              }
              if (length(object@prevalence) > 0) {
                cat("\n\nPrevalence:\t", object@prevalence)
              }
              #if (nrow(object@mahalanobis_proj) > 0) {
              #  cat("\n\nMahalanobis projection quantiles:")
              #  mq <- lapply(1:ncol(object@mahalanobis_proj), function(x) {
              #    cat("\n", colnames(object@mahalanobis_proj)[x], "\n")
              #    print(quantile(object@mahalanobis_proj[, x]))
              #  })
              #}
              #if (class(object@projections_maha)[1] %in%
              #    c("RasterLayer", "RasterStack")) {
              #  cat("\n\nRaster mahalanobis projection:\n")
              #  print(summary(object@projections_maha))
              #}
              #if (nrow(object@suitability_proj) > 0) {
              #  cat("\n\nSuitability projection quantiles:")
              #  mq <- lapply(1:ncol(object@suitability_proj), function(x) {
              #    cat("\n", colnames(object@suitability_proj)[x], "\n")
              #    print(quantile(object@suitability_proj[, x]))
              #  })
              #}
              #if (class(object@projections_suit)[1] %in%
              #    c("RasterLayer", "RasterStack")) {
              #  cat("\n\nRaster mahalanobis projection:\n")
              #  print(summary(object@projections_suit))
              #}
            }
          }
)

#' @rdname summary
setMethod("summary", signature(object = "ellipsoid_model_rep"),
          function(object) {
            # -----------
            # detecting potential errors
            if (!missing(object)) {
              clo <- class(object)[1]
              if (clo != "ellipsoid_model_rep") {
                stop("Argument 'object' must be of class ellipsoid*.")
              }
            }else {
              stop("Argument 'object' is necessary.")
            }

            nam <- names(object@ellipsoids)
            if (length(grep("mean", nam)) > 0 & length(grep("min", nam)) > 0 &
                length(grep("max", nam)) > 0 ) {
              nums <- c(grep("mean", nam), grep("min", nam), grep("max", nam))
              ellipsoids <- object@ellipsoids[nums]
            } else {
              ellipsoids <- mmm_ellipsoid(object@ellipsoids)
            }

            cat("\n                  Summary of ellipsoid_model_rep object\n")
            cat("---------------------------------------------------------------------------\n")
            cat("\n---------------------------Ellipsoids metadata-----------------------------\n")
            cat("\nMethod:\t", ellipsoids[[1]]@method)
            cat("\n\nLevel:\t", ellipsoids[[1]]@level)
            cat("\n\nCentroids:\n")
            print(sapply(ellipsoids, function(x) {x@centroid}))
            cat("\n\nCovariance matrices:\n")
            cov <- lapply(1:length(ellipsoids), function(x) {
              cat(names(ellipsoids)[x], "\n")
              print(ellipsoids[[x]]@covariance_matrix); cat("\n")
            })
            cat("\n\nVolumes:\n")
            print(sapply(ellipsoids, function(x) {x@niche_volume}))
            #cat("\n\nSemi-axes lengths:\n")
            #print(sapply(ellipsoids, function(x) {x@semi_axes_length}))
            #cat("\n\nAxes coordinates:\n")
            #acor <- lapply(1:length(ellipsoids), function(x) {
            #  cat(names(ellipsoids)[x])
            #  a_cord <- ellipsoids[[x]]@axes_coordinates
            #  cords <- lapply(1:length(a_cord), function(x) {
            #    cat("\n", letters[x], "\n"); print(a_cord[[x]])
            #  })
            #  cat("\n")
            #})
          }
)

#' @rdname summary
setMethod("summary", signature(object = "calibration_ellipsoid"),
          function(object) {
            # -----------
            # detecting potential errors
            if (!missing(object)) {
              clo <- class(object)[1]
              if (clo != "calibration_ellipsoid") {
                stop("Argument 'object' must be of class calibration_ellipsoid.")
              }
            }else {
              stop("Argument 'object' is necessary.")
            }

            cat("\n                 Summary of calibration_ellipsoid object\n")
            cat("---------------------------------------------------------------------------\n")
            cat("\nMethods tested:\t", object@methods)
            cat("\n\nLevel used:\t", object@level)
            cat("\n\nVariable sets tested:\n")
            print(object@variable_sets[[2]])
            cat("Selection criteria:\t")
            if (object@selection_criteria == "S_OR") {
              cat("Significance and omission rate")
            }
            if (object@selection_criteria == "S_OR_P") {
              cat("Significance, omission rate, and prevalence")
            }
            cat("\n\nSelected parameters:\n")
            print(object@selected_parameters, row.names = FALSE)
          }
)

#' @rdname summary
setMethod("summary", signature(object = "data_overlap"),
          function(object) {
            # -----------
            # detecting potential errors
            if (!missing(object)) {
              clo <- class(object)[1]
              if (clo != "data_overlap") {
                stop("Argument 'object' must be of class data_overlap.")
              }
            }else {
              stop("Argument 'object' is necessary to perform the analysis.")
            }

            cat("\n                      Summary of data_overlap object\n")
            cat("---------------------------------------------------------------------------\n")
            cat("\nMethod:\t", object@method)
            cat("\n\nLevel:\t", object@level)
            cat("\n\nData:\n")
            print(summary(object@data))
            cat("\nVariables:\n")
            print(summary(object@variables))
          }
)

#' @rdname summary
setMethod("summary", signature(object = "overlap_ellipsoid"),
          function(object) {
            # -----------
            # detecting potential errors
            if (!missing(object)) {
              clo <- class(object)[1]
              if (clo != "overlap_ellipsoid") {
                stop("Argument 'object' must be of class overlap_ellipsoid.")
              }
            }else {
              stop("Argument 'object' is necessary to perform the analysis.")
            }

            nam <- names(object@ellipsoids)
            ellipsoids <- object@ellipsoids

            cat("\n                   Summary of overlap_ellipsoid object\n")
            cat("---------------------------------------------------------------------------\n")
            cat("\n---------------------------Ellipsoids metadata-----------------------------\n")
            cat("\nMethod:\t", ellipsoids[[1]]@method)
            cat("\n\nLevel:\t", ellipsoids[[1]]@level)
            cat("\n\nCentroids:\n")
            print(sapply(ellipsoids, function(x) {x@centroid}))
            cat("\n\nCovariance matrices:\n")
            cov <- lapply(1:length(ellipsoids), function(x) {
              cat(names(ellipsoids)[x], "\n")
              print(ellipsoids[[x]]@covariance_matrix); cat("\n")
            })
            cat("\n\nVolumes:\n")
            print(sapply(ellipsoids, function(x) {x@niche_volume}))
            cat("\n\n----------------------------Summary of results--------------------------\n")
            if (nrow(object@full_overlap) > 0) {
              cat("\nFull overlap results:\n"); print(object@full_overlap)
            }
            if (nrow(object@union_overlap) > 0) {
              cat("\nUnion overlap results:\n"); print(object@union_overlap)
            }
          }
)
