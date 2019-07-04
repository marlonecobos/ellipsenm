#' Fit ellipsoids based on distinct methods
#'
#' @description ellipsoid_fit helps in finding the centroid and matrix that
#' define an ellipsoid. It uses distinct methods with asumptions that differ
#' from each other.
#'
#' @param data data.frame of occurrence records. Columns must be: species,
#' longitude, and latitude. Optionally, if \code{raster_layers} is not defined,
#' \code{data} must include more columns containing the values of at least two
#' variables to be used for fitting ellipsoid* models.
#' @param species (character) name of the column with the name of the species.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @param raster_layers RasterStack of at least two environmental variables to be
#' extracted using geographic coordinates present in \code{data}. If not provided
#' data must include additional columns containing values of variables to fit
#' ellipsoid* models.
#' @param method (character) method to construct the ellipsoid that characterizes
#' the species ecological niche. Available methods are: "covmat", "mve1", and
#' "mve2". See details of \code{\link{ellipsoid_fit}}. Default = "covmat".
#' @param level (numeric) the confidence level of a pairwise confidence region
#' for the ellipsoid, expresed as percentage. Default = 95.
#' @param replicates (numeric) number of replicates to perform. Default = 1
#' produces a single model using all the data.
#' @param replicate_type (character) type of replicates to perform. Options are:
#' "bootstrap" and "jackknife"; default = "bootstrap". See details. Ignored if
#' \code{replicates} = 1.
#' @param bootstrap_percentage (numeric) percentage of data to be bootstrapped
#' for each replicate. Default = 50. Valid if \code{replicates} > 1 and
#' \code{replicate_type} = "bootstrap".
#' @param projection_variables optional: (character) name of folder containing other
#' folders with raster layers that represent the scenarios for projections; a
#' RasterStack of layers respresenting an only scenario for projection; or a
#' named list of RasterStacks representing multiple scenarios for projection.
#' See details. Default = NULL.
#' @param prediction (character) type of prediction to be made, options are:
#' "suitability", "mahalanobis", and "both". Default = "suitability".
#' @param return_numeric (logical) whether or not to return values of mahalanobis
#' distance and suitability as part of the results (it depends on the type of
#' \code{prediction} selected). Default = FALSE.
#' @param tolerance the tolerance for detecting linear dependencies.
#' Default = 1e-60.
#' @param format (charater) file type for raster outputs to be written in
#' \code{output_directory}. Default = "GTiff". See \code{\link[raster]{writeFormats}}.
#' @param overwrite (logical) whether or not to overwrite exitent results in
#' \code{output_directory}. Default = FALSE.
#' @param color_palette a color palette function to be used in plotting
#' suitability values in an HTML report produced at the end of all analyses.
#' @param output_directory name of the folder were all results will be written.
#' This avoids saturation of the RAM.
#'
#' @return
#'
#' @details
#' \code{replicate_type}
#'
#' \code{projection_variables}
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
#' # creating the model with no replicates
#' ell_model <- ellipsoid_model(data = occurrences, species = "species",
#'                              longitude = "longitude", latitude = "latitude",
#'                              raster_layers = vars, method = "mve1", level = 99,
#'                              replicates = 1, prediction = "suitability",
#'                              return_numeric = TRUE, format = "GTiff",
#'                              overwrite = FALSE, output_directory = "ellipsenm_model")
#'
#' class(ell_model)
#' # check your directory, folder "ellipsenm_model"

ellipsoid_model <- function (data, species, longitude, latitude, raster_layers,
                             method = "covmat", level = 95, replicates = 1,
                             replicate_type = "bootstrap", bootstrap_percentage = 75,
                             projection_variables = NULL, prediction = "suitability",
                             return_numeric = TRUE, tolerance = 1e-60, format = "GTiff",
                             overwrite = FALSE, color_palette = viridis::magma,
                             output_directory = "ellipsenm_model") {
  # -----------
  # detecting potential errors, other potential problems tested in code
  if (missing(data)) {
    stop("Argument occurrences is necessary to perform the analysis.")
  }
  if (missing(species)) {
    stop("Argument species is not defined.")
  }
  if (missing(longitude)) {
    stop("Argument longitude is not defined.")
  }
  if (missing(latitude)) {
    stop("Argument latitude is not defined.")
  }
  if (missing(raster_layers)) {
    variables <- data[, !colnames(data) %in% c(species, longitude, latitude)]
    if (ncol(variables) < 2) {
      stop("If raster_layers is not defined, data must contain information of at least\ntwo variables to fit ellipsoids. See function's help.")
    }
  }
  if (overwrite == FALSE & dir.exists(output_directory)) {
    stop("output_directory already exists, to replace it use overwrite = TRUE.")
  }
  if (overwrite == TRUE & dir.exists(output_directory)) {
    unlink(x = output_directory, recursive = TRUE, force = TRUE)
  }

  # -----------
  # preparing data and variables
  cat("\nPreparing data...\n")
  sp <- as.character(data[1, species])

  if (!missing(raster_layers)) {
    data <- cbind(data, raster::extract(raster_layers, data[, c(longitude, latitude)]))
    raster_base <- raster_layers[[1]]
    nona <- !is.na(raster::values(raster_base))
    variable1 <- raster_layers[[1]]
    variable_names <- names(raster_layers)
    r_values <- na.omit(raster::values(raster_layers))
    nb <- nrow(r_values)
    n_prop <- ifelse(nb > 100000, 0.1, 0.3)
    set.seed(1)
    samp <- sample(nb, ceiling(nb * n_prop))
    r_values <- r_values[samp, ]
    variables <- raster_layers
  } else {
    data <- cbind(data[, c(species, longitude, latitude)], variables)
    variable_names <- colnames(variables)
  }

  n_var <- length(variable_names)
  mpos <- replicates + 1

  # -----------
  # fitting ellipsoids and getting statistics
  if (replicates >= 1) {
    data1 <- data_subsample(data[, -1], replicates, replicate_type, bootstrap_percentage)
  } else {
    stop("Argument replicates needs to be numeric and >= 1, see function's help.")
  }

  cat("\nFitting ellipsoids using occurrence data:\n")
  ellipsoids <- lapply(1:replicates, function(x){
    cat("\tFitting ellipsoid for replicate", x, "of", length(data1), "\n")
    ellipsoid_fit(data1[[x]], longitude, latitude, method, level)
  })
  names(ellipsoids) <- paste0("replicate", 1:replicates)

  if (replicates > 1) {
    ellipsoids <- new("ellipsoid_model_rep",
                      ellipsoids = c(ellipsoids, mmm_ellipsoid(ellipsoids)))
  } else {
    ellipsoids <- ellipsoids[[1]]
  }

  # -----------
  # prediction in calibration area
  cat("\nPreparing raster predictions for calibration area:\n")
  nam_format <- rformat_type(format)
  namer <- paste0(output_directory, "/calibration_", sp, nam_format)
  dir.create(output_directory)

  if (replicates > 1) {
    predictions <- predict(ellipsoids, variables, prediction, return_numeric,
                           tolerance, namer, format, overwrite,
                           force_return = TRUE, return_name = "mean_ellipsoid")
  } else {
    predictions <- predict(ellipsoids, variables, prediction, return_numeric,
                           tolerance, namer, format, overwrite, force_return = TRUE)
  }

  # -----------
  # model projections
  if (!is.null(projection_variables)) {
    cat("\nProducing results for projection  scenario(s):\n")
    projections <- model_projection(predictions, projection_variables, sp,
                                    return_numeric, prediction, tolerance,
                                    format, overwrite, output_directory)
  }

  # -----------
  # returning metadata and preparing needed variables
  cat("\nObtaining and writing metadata of ellipsoid models and predictions:\n")
  cat("\tMetadata for ellipsoid models\n")
  ell_meta <- write_ellmeta(predictions,
                            name = paste0(output_directory, "/ellipsoid_metadata"))

  if (is.null(projection_variables)) {
    if (prediction != "mahalanobis") {
      layer <- predictions@prediction_suit
      if (class(predictions)[1] == "ellipsoid_model_rep") {
        mean_pred <- predictions@ellipsoids[[mpos]]
        prevalences <- predictions@prevalence
      } else {
        mean_pred <- predictions
        prevalences <- data.frame(ellipsoid_model = predictions@prevalence)
      }
      write.csv(prevalences, paste0(output_directory, "/calibration_prevalences.csv"),
                row.names = TRUE)
      if (prediction == "both") {
        layer <- raster::stack(layer, predictions@prediction_maha)
      }
      cat("\tPrevalence in calibration area\n")
    } else {
      layer <- predictions@prediction_maha
      if (class(predictions)[1] == "ellipsoid_model_rep") {
        mean_pred <- predictions@ellipsoids[[mpos]]
      } else {
        mean_pred <- predictions
      }
      prevalences <- vector()
    }
  } else {

  }


  # -----------
  # producing report
  cat("\nAnalyses finished. Producing HTML report...\n")
  save(data, variable_names, variable1, n_var, r_values, ell_meta, mean_pred,
       layer, prevalences, replicates, replicate_type, bootstrap_percentage,
       color_palette, file = paste0(output_directory, "/enm_report_data.RData"))

  report_format(name = paste0(output_directory, "/report_format"))
  if (prediction != "mahalanobis") {
    if (prediction == "both") {
      report("enm_both", name = paste0(output_directory, "/enm_both_report"))
    } else {
      report("enm_suitability", name = paste0(output_directory, "/enm_suitability_report"))
    }
  } else {
    report("enm_mahalanobis", name = paste0(output_directory, "/enm_mahalanobis_report"))
  }

  # -----------
  # returning results
  return(predictions)
}
