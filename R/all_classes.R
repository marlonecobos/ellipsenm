# R classes for ellipsenm objects
# May 2019
# Version 0.0.1
# Licence GPL v3

ellipsoid_model <- setClass("ellipsoid_model",
                            slots = c(method = "character",
                                      occ_predictions = "data.frame",
                                      centroid = "numeric",
                                      cov_matrix = "matrix",
                                      predictions = "RasterLayer",
                                      projections = "list"))

ellipsoid_overlap <- setClass("ellipsoid_overlap",
                              slots = c(method = "character",
                                        data = "list",
                                        models = "list",
                                        overlap = "data.frame"))

ellipsoid_calibration <- setClass("ellipsoid_calibration",
                                  slots = c(centroid_types = "character",
                                            covariance_types = "character",
                                            variable_sets = "list",
                                            results = "data.frame"))
