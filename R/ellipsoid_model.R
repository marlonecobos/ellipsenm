#' Fit ellipsoids based on distinct methods
#'
#' @description ellipsoid_fit helps in finding the centroid and matrix that
#' define an ellipsoid. It uses distinct methods with asumptions that differ
#' from each other.
#'
#' @param data data.frame of occurrence records. Columns must be: species,
#' longitude, and latitude.
#' @param species (character) name of the column with the name of the species.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @param raster_layers RasterStack of environmental variables to be extracted
#' using geographic coordinates present in \code{data}.
#' @param method (character) method to construct the ellipsoid that characterizes
#' the species ecological niche. Available methods are: "covmat", "mve1", and
#' "mve2". See details of \code{\link{ellipsoid_fit}}. Default = "mve1".
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
#' @param projection_layers optional: (character) name of folder containing other
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
#' \code{output_directory}. Default = = "GTiff". See \code{\link[raster]{writeFormats}}.
#' @param overwrite (logical) whether or not to overwrite exitent results in
#' \code{output_directory}. Default = FALSE.
#' @param output_directory name of the folder were all results will be written.
#' This avoids saturation of the RAM.
#'
#' @return
#'
#' @details
#' \code{replicate_type}
#'
#' \code{projection_layers}
#'
#' @export
#'
#' @examples
#' occurrences <- read.csv(system.file("extdata", "occurrences_comp.csv",
#'                                     package = "ellipsenm"))
#'
#' vars <- raster::stack(list.files(system.file("extdata", package = "ellipsenm"),
#'                                  pattern = "m_bio", full.names = TRUE))


ellipsoid_model <- function (data, species, longitude, latitude, raster_layers,
                             method = "mve1", level = 95, replicates = 1,
                             replicate_type = "bootstrap", bootstrap_percentage = 50,
                             projection_layers = NULL, prediction = "suitability",
                             return_numeric = TRUE, tolerance = 1e-60, format = "GTiff",
                             overwrite = FALSE, output_directory = "ellipsenm_model") {
  # -----------
  # detecting potential errors, other potential problems tested in code
  if (missing(data)) {
    stop("Argument occurrences is necessary to perform the analysis.")
  }
  if (missing(longitude)) {
    stop("Argument longitude is not defined.")
  }
  if (missing(latitude)) {
    stop("Argument latitude is not defined.")
  }
  if (overwrite == FALSE & dir.exists(output_directory)) {
    stop("output_directory already exists, to replace it use overwrite = TRUE.")
  }
  if (overwrite == TRUE & dir.exists(output_directory)) {
    unlink(x = output_directory, recursive = TRUE, force = TRUE)
  }

  # -----------
  # preparing data and variables
  cat("\nPreparing data:\n")
  sp <- as.character(data[1, species])
  data <- cbind(data, raster::extract(raster_layers, data[, c(longitude, latitude)]))

  raster_base <- raster_layers[[1]]
  r_values <- na.omit(raster::values(raster_layers))

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
  # prediction in raster_layers
  cat("\nPreparing raster predictions for calibration area:\n")
  nam_format <- rformat_type(format)
  name <- paste0(output_directory, "/", sp, nam_format)

  predictions <- predict(ellipsoids, raster_layers, prediction, return_numeric,
                         tolerance, name, format, overwrite) # check this part force_raster_in

  if (prediction != "mahalanobis") {
    prevalences <- do.call(cbind, lapply(predictions, function(x) {x@prevalence}))
  } else {
    prevalences <- vector()
  }

  cat("\nObtaining and writing metadata of ellipsoid models:\n")
  ell_meta <- do.call(cbind, lapply(ellipsoids, function(x) {
    c(x@method, x@level, round(x@niche_volume, digits = 2))
  }))
  ell_meta <- as.data.frame(rbind(ell_meta, paste0(c(1:replicates, "mean", "min", "max"),
                                                   "_ell_meta.txt")))
  row.names(ell_meta) <- c("Method", "Level", "Niche volume", "Other characteristics")
  write.csv(ell_meta, paste0(output_directory, "/ellipsoid_metadata_summary.csv"),
            row.names = TRUE)



  #predictions <- predictions[[(replicates + 1)]]
  if (is.null(projection_layers)) {

  } else {
    cclas <- class(projection_layers)[1]
    if (cclas == "character" | cclas == "RasterStack" | cclas == "list") {
      if (cclas == "RasterStack") {
        if (!all(names(raster_layers) == names(projection_layers))) {
          stop("Variable names of projection_layers do not match names of raster_layers.")
        }
      }

      if (cclas == "list") {
        for (i in 1:length(projection_layers)) {
          if (!all(names(raster_layers) == names(projection_layers[[i]]))) {
            stop("Variable names of projection_layers do not match names of raster_layers.")
          }
        }
      }

      if (cclas == "character") {
        dirs <- dir(projection_layers, full.names = TRUE)
        for (i in 1:length(dirs)) {
          namest <- list.files(dirs[i], pattern = format_pl)
          if (!all(names(raster_layers) == namest)) {
            stop("Variable names of projection_layers do not match names of raster_layers.")
          }
        }
      }
    } else {
      stop("Argument projection_layers is not valid, see function's help.")
    }

  }

  # -----------
  # performing analyses
  if (is.null(projection_layers)) {
    ## fitting ellipsoids
    ellipsoid_fit(data, longitude, latitude, method, level, raster_layers)

    ## predicting


    # -----------
    # preparing further results
    ## model statistics (median, range, min, max)
    ## tables
    ## statistis of performance
    ## plots

    # -----------
    # writing results

  } else {
    ## fitting ellipsoids
    ellipsoid_fit(data, longitude, latitude, method, level, raster_layers)

    ## predicting


    # -----------
    # preparing further results
    ## model statistics (median, range, min, max)
    ## tables
    ## statistis of performance
    ## plots

    # -----------
    # writing results and savind objects for report


    ## objects for report
    save(ellipsoids, file = "enm_report_data.RData")
  }


  # -----------
  # producing report
  save(data, variable_names, variable1, n_var, back1, ell_meta, predictions,
       prevalences, file = paste0(output_directory, "/enm_report_data.RData"))

  report_format(name = paste0(output_directory, "/eenm_report_format"))
  report("enm", name = paste0(output_directory, "/eenm_results_report"))

  unlink(paste0(output_directory, "/eenm_report_format.css"))

  # -----------
  # returning results

  return(ellipsoid_m)
}


