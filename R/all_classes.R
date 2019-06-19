# R classes for ellipsenm objects
# May 2019
# Version 0.1.3
# Licence GPL v3

ellipsoid_mod <- setClass("ellipsoid_model",
                          slots = c(data_predictions = "data.frame",
                                    prediction = "RasterLayer",
                                    projections = "list"),
                          contains = "ellipsoid_fit")

ellipsoid <- setClass("ellipsoid_fit",
                      slots = c(method = "character",
                                centroid = "numeric",
                                covariance_matrix = "matrix",
                                level = "numeric",
                                niche_volume = "numeric",
                                semi_axes_length = "numeric",
                                axis_coordinates = "list"))

ellipsoid_cal <- setClass("ellipsoid_calibration",
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

ellipsoid_ove <- setClass("ellipsoid_overlap",
                          slots = c(method = "character",
                                    data = "list",
                                    models = "list",
                                    overlap = "data.frame"))
