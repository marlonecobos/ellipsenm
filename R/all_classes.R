# R classes for ellipsenm objects
# May 2019
# Version 0.0.1
# Licence GPL v3

ellipsoid_model <- setClass("ellipsoid_model",
                            slots = c(algorithm = "character",
                                      occ_points = "data.frame",
                                      occ_predictions = "numeric",
                                      bg_points = "data.frame",
                                      bgpredictions = "numeric",
                                      centroid = "numeric",
                                      cov_matrix = "matrix",
                                      predictions = "RasterStack",
                                      projections = "list"))
