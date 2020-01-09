# R classes for ellipsenm objects
# May 2019
# Version 0.3.2
# Licence GPL v3

#' An S4 class to organize data and results of ellipsoid* objects
#' @name ellipsoid
#' @aliases ellipsoid-class ellipsoid_model_sim-class ellipsoid_model_rep-class
#' @slot method ("character"), not in ellipsoid_model_rep.
#' @slot level ("numeric"), not in ellipsoid_model_rep.
#' @slot centroid ("numeric"), not in ellipsoid_model_rep.
#' @slot covariance_matrix ("matrix"), not in ellipsoid_model_rep.
#' @slot niche_volume ("numeric"), not in ellipsoid_model_rep.
#' @slot semi_axes_length ("numeric"), not in ellipsoid_model_rep.
#' @slot axes_coordinates list of matrices
#' @slot ellipsoids list of objects of class ellipsoid, only ellipsoid_model_rep.
#' @slot mahalanobis if ellipsoid_model_sim, ("numeric"); if ellipsoid_model_rep,
#' numeric matrix.
#' @slot suitability if ellipsoid_model_sim, ("numeric"); if ellipsoid_model_rep,
#' numeric matrix.
#' @slot prevalence if ellipsoid_model_sim, ("numeric"); if ellipsoid_model_rep,
#' numeric matrix.
#' @slot prediction_maha object of class Raster*.
#' @slot prediction_suit object of class Raster*.
#' @slot mahalanobis_proj if ellipsoid_model_sim, numeric matrix; if
#' ellipsoid_model_rep, list of numeric matrices.
#' @slot suitability_proj if ellipsoid_model_sim, numeric matrix; if
#' ellipsoid_model_rep, list of numeric matrices.
#' @slot projections_maha if ellipsoid_model_sim, object of class Raster*, if
#' ellipsoid_model_rep, list of Raster* objects.
#' @slot projections_suit if ellipsoid_model_sim, object of class Raster*, if
#' ellipsoid_model_rep, list of Raster* objects.
#' @export
#' @examples
#' showClass("ellipsoid")
#' #' showClass("ellipsoid_model_rep")
#' @rdname ellipsoid
ellipsoid <- setClass("ellipsoid",
                      slots = c(method = "character",
                                level = "numeric",
                                centroid = "numeric",
                                covariance_matrix = "matrix",
                                niche_volume = "numeric",
                                semi_axes_length = "numeric",
                                axes_coordinates = "list"))

#' @rdname ellipsoid
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

#' @rdname ellipsoid
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


#' An S4 class to organize data and results of calibration_ellipsoid objects
#' @name calibration_ellipsoid
#' @aliases calibration_ellipsoid-class
#' @slot method ("character")
#' @slot data list of data.frames with the complete set of occurrences, the
#' training set, and the testing set.
#' @slot variable_sets list of two objects, a RasterStack with variables, and a
#' list of character vectors with sets of variables.
#' @slot level ("numeric")
#' @slot results a data.frame of all results of the model evaluation.
#' @slot selection_criteria ("character")
#' @slot selected_parameters data.frame of selected paramete settings.
#' @export
#' @examples
#' showClass("calibration_ellipsoid")
calibration_ellipsoid <- setClass("calibration_ellipsoid",
                                  slots = c(methods = "character",
                                            data = "list",
                                            variable_sets = "list",
                                            level = "numeric",
                                            results = "data.frame",
                                            selection_criteria = "character",
                                            selected_parameters = "data.frame"))


#' An S4 class to organize data for data_overlap objects
#' @name data_overlap
#' @aliases data_overlap-class
#' @slot data data.frame
#' @slot main_columns ("character")
#' @slot method ("character")
#' @slot level ("numeric")
#' @slot results a data.frame of all results of the model evaluation.
#' @slot variables object of class RasterStack.
#' @export
#' @examples
#' showClass("data_overlap")
data_overlap <- setClass("data_overlap",
                         slots = c(data = "data.frame",
                                   main_columns = "character",
                                   method = "character",
                                   level = "numeric",
                                   variables = "S4"))


#' An S4 class to organize data and results of overlap_ellipsoid objects
#' @name overlap_ellipsoid
#' @aliases overlap_ellipsoid-class
#' @slot ellipsoids list of objects of class ellipsoid.
#' @slot data list of data.frames with the data used for analyses.
#' @slot variable_names ("character")
#' @slot full_background list of matrices representing background for full
#' overlap.
#' @slot full_overlap a data.frame of overlap results considering full background.
#' @slot union_background list of matrices representing background for union
#' overlap.
#' @slot union_overlap data.frame of overlap results considering union background.
#' @slot significance_results list of all overlap results resulted from
#' significance random tests.
#' @export
#' @examples
#' showClass("overlap_ellipsoid")
overlap_ellipsoid <- setClass("overlap_ellipsoid",
                              slots = c(ellipsoids = "list",
                                        data = "list",
                                        variable_names = "character",
                                        full_background = "list",
                                        full_overlap = "data.frame",
                                        union_background = "list",
                                        union_overlap = "data.frame",
                                        significance_results = "list"))
