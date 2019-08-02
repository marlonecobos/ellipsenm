#' Calibration of ellipsoid-based ecological niche models
#'
#' @description ellipsoid_calibration helps in creating and evaluating multiple
#' candidate ellipsoid envelop models to find parameter settings that produce
#' the best results.
#'
#' @param data (character or list) if character, vector of names of csv files
#' containing all, training, and testing occurrences located in the working
#' directory; if list, object resulted from \code{\link{split_data}}. Columns of
#' tables must include: species, longitude, and latitude.
#' @param species (character) name of the column with the name of the species.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @param variables (character or list) if character, name of a folder containing
#' subfolders of at least one set of variables; if list, object derived from
#' \code{\link{prepare sets}}. Sets of variables must contain at least two layers.
#' @param methods (character) methods to construct the ellipsoid ecological niche
#' models to be tested. Available methods are: "covmat", "mve1", and "mve2".
#' See details of \code{\link{ellipsoid_fit}}.
#' @param level (numeric) the confidence level of a pairwise confidence region
#' for the ellipsoid, expresed as percentage. Default = 95.

ellipsoid_calibration <- function(data, species, longitude, latitude, variables,
                                  methods, level = 95, parallel = FALSE,
                                  overwrite = FALSE,
                                  output_directory = "calibration_results") {
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
  if (missing(variables)) {
    stop("Argument variables is necessary to perform the analysis.")
  }
  if (missing(methods)) {
    stop("Argument methods is necessary to perform the analysis.")
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
  clocc <- class(data)[1]
  if (clocc == "character" | clocc == "list") {
    if (clocc == "character") {
      data <- lapply(data, function(x) {
        read.csv(x)
      })
    }
    sp <- as.character(data[[1]][1, species])

  } else {
    stop("data must be either character or list. See function's help.")
  }




  # -----------
  # producing report
  report_format(name = paste0(output_directory, "/report_format"))
  report(report_type = "calibration", output_directory)

  # -----------
  # returning results
  results <- calibration_ellipsoid(methods = "character",
                                   data = "data.frame",
                                   variable_sets = "list",
                                   level = "numeric",
                                   results = "data.frame",
                                   selection_criteria = "character",
                                   selected_parameters = "character")

  return(results)
}
