# R classes for ellipsenm objects
# May 2019
# Version 0.1.5
# Licence GPL v3
#' @export

ellipsoid <- setClass("ellipsoid",
                      slots = c(method = "character",
                                level = "numeric",
                                centroid = "numeric",
                                covariance_matrix = "matrix",
                                niche_volume = "numeric",
                                semi_axes_length = "numeric",
                                axes_coordinates = "list"))

ellipsoid_model_sim <- setClass("ellipsoid_model_sim",
                                slots = c(mahalanobis = "numeric",
                                          suitability = "numeric",
                                          prevalence = "numeric",
                                          prediction_maha = "S4",
                                          prediction_suit = "S4",
                                          mahalanobis_proj = "matrix",
                                          suitability_proj = "matrix",
                                          projections_maha = "S4",
                                          projections_suit = "S4"),
                                contains = "ellipsoid")

ellipsoid_model_rep <- setClass("ellipsoid_model_rep",
                                slots = c(ellipsoids = "list",
                                          mahalanobis = "matrix",
                                          suitability = "matrix",
                                          prevalence = "matrix",
                                          prediction_maha = "S4",
                                          prediction_suit = "S4",
                                          mahalanobis_proj = "list",
                                          suitability_proj = "list",
                                          projections_maha = "list",
                                          projections_suit = "list"))

calibration_ellipsoid <- setClass("calibration_ellipsoid",
                                  slots = c(methods = "character",
                                            data = "list",
                                            variable_sets = "list",
                                            level = "numeric",
                                            results = "data.frame",
                                            selection_criteria = "character",
                                            selected_parameters = "data.frame"))

data_overlap <- setClass("data_overlap",
                         slots = c(data = "data.frame",
                                   method = "character",
                                   level = "numeric",
                                   variables = "S4"))

overlap_ellipsoid <- setClass("overlap_ellipsoid",
                              slots = c(ellipsoids = "list",
                                        spp_data = "list",
                                        background = "list",
                                        overlap = "data.frame"))
