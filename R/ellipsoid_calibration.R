ellipsoid_calibration <- function(data, species, longitude, latitude,
                                  variables, methods, overwrite = FALSE,
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
      sp <- as.character(data[[1]][1, species])
    } else {
      sp <- as.character(data[[1]][1, species])
    }

  } else {
    stop("class of data is not compatible. See function's help.")
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
