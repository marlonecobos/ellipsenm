# R classes for ellipsenm objects
# May 2019
# Version 0.1.3
# Licence GPL v3

ellipsoid_mod <- setClass("ellipsoid_model",
                          slots = c(method = "character",
                                    occ_predictions = "data.frame",
                                    centroid = "numeric",
                                    covariance_matrix = "matrix",
                                    niche_volume = "numeric",
                                    semi_axes_length = "numeric",
                                    axis_coordinates = "list",
                                    predictions = "RasterLayer",
                                    projections = "list"))

ellipsoid <- setClass("ellipsoid",
                      slots = c(method = "character",
                                centroid = "numeric",
                                covariance_matrix = "matrix",
                                niche_volume = "numeric",
                                semi_axes_length = "numeric",
                                axis_coordinates = "list"))

ellipsoid_cal <- setClass("ellipsoid_calibration",
                          slots = c(centroid_types = "character",
                                    covariance_types = "character",
                                    variable_sets = "list",
                                    results = "data.frame"))

data_overlap <- setClass("data_overlap",
                         slots = c(method = "character",
                                   occurrences = "data.frame",
                                   m_layers = "RasterStack"))

ellipsoid_ove <- setClass("ellipsoid_overlap",
                          slots = c(method = "character",
                                    data = "list",
                                    models = "list",
                                    overlap = "data.frame"))
