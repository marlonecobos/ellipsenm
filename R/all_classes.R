# R classes for ellipsenm objects
# May 2019
# Version 0.1.5
# Licence GPL v3
#' @export

ellipsoid <- setClass("ellipsoid",
                      slots = c(method = "character",
                                centroid = "numeric",
                                covariance_matrix = "matrix",
                                level = "numeric",
                                niche_volume = "numeric",
                                semi_axes_length = "numeric",
                                axes_coordinates = "list"))

ellipsoid_model_sim <- setClass("ellipsoid_model_sim",
                              slots = c(mahalanobis = "numeric",
                                        suitability = "numeric",
                                        prevalence = "numeric",
                                        prediction_maha = "RasterLayer",
                                        prediction_suit = "RasterLayer",
                                        projections_maha = "character",
                                        projections_suit = "character"),
                              contains = "ellipsoid")

ellipsoid_model_rep <- setClass("ellipsoid_model_rep",
                              slots = c(ellipsoids = "list",
                                        mahalanobis = "matrix",
                                        suitability = "matrix",
                                        prevalence = "matrix",
                                        prediction_maha = "RasterStack",
                                        prediction_suit = "RasterStack",
                                        projections_maha = "list",
                                        projections_suit = "list"))

calibration_ellipsoid <- setClass("calibration_ellipsoid",
                                  slots = c(methods = "character",
                                            data = "data.frame",
                                            variable_sets = "list",
                                            level = "numeric",
                                            results = "data.frame",
                                            selection_criteria = "character",
                                            selected_parameters = "character"))

data_overlap <- setClass("data_overlap",
                         slots = c(method = "character",
                                   data = "data.frame",
                                   raster_layers = "RasterStack",
                                   level = "numeric"))

overlap_ellipsoid <- setClass("overlap_ellipsoid",
                              slots = c(method = "character",
                                        data = "list",
                                        models = "list",
                                        overlap = "data.frame"))
