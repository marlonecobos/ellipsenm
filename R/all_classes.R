# R classes for ellipsenm objects
# May 2019
# Version 0.1.5
# Licence GPL v3
#' @export

ellipsoid <- setClass("ellipsoid",
                      slots = c(method = "character",
                                level = "numeric"))

ellipsoid_basic <- setClass("ellipsoid_basic",
                            slots = c(centroid = "numeric",
                                      covariance_matrix = "matrix",
                                      niche_volume = "numeric",
                                      semi_axes_length = "numeric",
                                      axes_coordinates = "list"),
                            contains = c("ellipsoid", "VIRTUAL"))

ellipsoid_model_sim <- setClass("ellipsoid_model_sim",
                              slots = c(mahalanobis = "numeric",
                                        suitability = "numeric",
                                        prevalence = "numeric",
                                        prediction_maha = "S4",
                                        prediction_suit = "S4",
                                        projections_maha = "character",
                                        projections_suit = "character"),
                              contains = "ellipsoid_basic")

ellipsoid_model_rep <- setClass("ellipsoid_model_rep",
                                slots = c(ellipsoids = "list",
                                          mahalanobis = "matrix",
                                          suitability = "matrix",
                                          prevalence = "matrix",
                                          prediction_maha = "S4",
                                          prediction_suit = "S4",
                                          projections_maha = "list",
                                          projections_suit = "list"),
                                contains = c("ellipsoid", "VIRTUAL"))

calibration_ellipsoid <- setClass("calibration_ellipsoid",
                                  slots = c(methods = "character",
                                            data = "data.frame",
                                            variable_sets = "list",
                                            level = "numeric",
                                            results = "data.frame",
                                            selection_criteria = "character",
                                            selected_parameters = "character"))

data_overlap <- setClass("data_overlap",
                         slots = c(data = "data.frame",
                                   raster_layers = "S4"),
                         contains = "ellipsoid")

overlap_ellipsoid <- setClass("overlap_ellipsoid",
                              slots = c(spp_data = "list",
                                        ellipsoids = "list",
                                        overlap = "data.frame"))
