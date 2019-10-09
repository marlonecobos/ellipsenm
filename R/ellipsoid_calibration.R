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
#' \code{\link{prepare_sets}}. Sets of variables must contain at least two layers.
#' @param format_in (character) if \code{variables} is character, format of the
#' variables found in folders. Default = NULL.
#' @param methods (character) methods to construct the ellipsoid ecological niche
#' models to be tested. Available methods are: "covmat", "mve1", and "mve2".
#' See details of \code{\link{ellipsoid_fit}}.
#' @param level (numeric) the confidence level of a pairwise confidence region
#' for the ellipsoid, expresed as percentage. Default = 95.
#' @param selection_criteria (character) set of criteria to select best models,
#' options are: "S_OR" (statistical significance and low omission) and
#' "S_OR_P" (statistical significance, low omission, and low prevalence).
#' See details. Default = "S_OR_P".
#' @param error (numeric) value from 0 to 100 to represent the percentage of
#' potential error (E) that the data could have due to any source of uncertainty.
#' Default = 5.
#' @param iterations (numeric) number of bootstrap iterations to be performed;
#' default = 500.
#' @param percentage (numeric) percentage of testing data to be used in each
#' bootstrapped process for calculating the partial ROC. Default = 50.
#' @param parallel (logical) whether or not to run analyses in parallel. If
#' defined as TRUE, it will only run in parallel if the number of parameter
#' settings to be tested is equal or larger than the number of cores available.
#' @param overwrite (logical) whether or not to overwrite exitent results in
#' \code{output_directory}. Default = FALSE.
#' @param output_directory (character) name of the folder were results of model
#' calibration and selection will be written.
#'
#' @return
#' A calibration_ellipsoid object with all results and details derived from the
#' calibration process. A folder named \code{output_directory}, containing all
#' results as well as a detailed HTML report will also be created.
#'
#' @export
#'
#' @details
#' Statistical significance is assessed using the partial_ROC test, omission
#' rates refer to the proportion of testing data known to be in suitable areas
#' but predicted as unsuitable, and prevalence is the proportion of geographic
#' and environmental space predicted as suitable.
#'
#' The maximum expected omission rates are 1 - (\code{level} / 100). Thus, if
#' \code{level} is determined as 95, an adequate omission rate should not be
#' higher than 0.5.
#'
#' Good models are expected to have low omission rates and low prevalence. This
#' implies that the model is predicting correctly in smaller areas.
#'
#' @examples
#' # reading data
#' occurrences <- read.csv(system.file("extdata", "occurrences.csv",
#'                                     package = "ellipsenm"))
#' colnames(occurrences)
#'
#' # raster layers of environmental data (this ones are masked to the accessible area)
#' # users must prepare their layers accordingly if using other data
#' vars <- raster::stack(list.files(system.file("extdata", package = "ellipsenm"),
#'                                  pattern = "bio", full.names = TRUE))
#'
#' # preparing training and testing data
#' data_split <- split_data(occurrences, method = "random", longitude = "longitude",
#'                          latitude = "latitude", train_proportion = 0.75)
#'
#' # sets of variables (example)
#' sets <- list(set_1 = c("bio_1", "bio_7", "bio_15"),
#'              set_2 = c("bio_1", "bio_12", "bio_15")) # change as needed
#'
#' variable_sets <- prepare_sets(vars, sets)
#'
#' # methods to create ellipsoids
#' methods <- c("covmat")
#'
#' # model calibration process (Make sure to define your working directory first)
#' calib <- ellipsoid_calibration(data = data_split, species = "species",
#'                                longitude = "longitude", latitude = "latitude",
#'                                variables = variable_sets, methods = methods,
#'                                level = 99, selection_criteria = "S_OR_P",
#'                                error = 5, iterations = 500, percentage = 50,
#'                                output_directory = "calibration_results")

ellipsoid_calibration <- function(data, species, longitude, latitude, variables,
                                  format_in = NULL, methods, level = 95,
                                  selection_criteria = "S_OR_P",
                                  error = 5, iterations = 500, percentage = 50,
                                  parallel = FALSE, overwrite = FALSE,
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
  cl_var <- class(variables)[1]
  if (!cl_var %in% c("character", "list")) {
    stop("variables must of class character or list, see function's help.")
  } else {
    if (cl_var == "character") {
      if (is.null(format_in)) {
        stop("Argument fomat_in cannot be NULL if variables is of class character.")
      }
    }
  }

  # -----------
  # preparing data and variables
  cat("\nPreparing data...\n")
  ## occurrences
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

  ## variables
  if (cl_var == "character") {
    patt <- paste0(rformat_type(format_in), "$")
    vars <- list.files(variables, pattern = patt, full.names = TRUE,
                       recursive = TRUE)
    varss <- list.files(variables, pattern = patt, recursive = TRUE)
    places <- !duplicated(varss)

    dirs <- dir(variables)
    sets <- lapply(1:length(dirs), function(x) {
      vs <- list.files(dirs[x], pattern = patt)
      gsub(patt, "", vs)
    })
    names(sets) <- paste0("set_", 1:length(sets))

    variables <- list(raster_layers = raster::stack(vars[places]),
                      variable_sets = sets)
  }

  # -----------
  # calibration process
  dir.create(output_directory)

  cat("\nRunning model calibration:\n")
  settings <- length(methods) * length(variables[[2]])
  if (parallel == FALSE) {
    calibration <- sapply(1:length(variables[[2]]), function(x) {
      varss <- variables[[1]][[variables[[2]][[x]]]]
      var_vals <- na.omit(raster::values(varss))

      resul <- lapply(1:length(methods), function(y) {
        ## pROC
        xy_train <- data[[2]][, c(longitude, latitude)]
        train_data <- data.frame(xy_train, na.omit(raster::extract(varss, xy_train)))
        occ_vals <- na.omit(raster::extract(varss, data[[3]][, c(longitude, latitude)]))

        ellip <- ellipsoid_fit(train_data, longitude, latitude, methods[y], level)

        pred <- predict(ellip, projection_variables = var_vals, truncate = FALSE)
        pred_occ <- predict(ellip, projection_variables = occ_vals)

        proc <- partial_roc(pred@suitability, pred_occ@suitability, longitude,
                            latitude, error, iterations, percentage)

        ## Omission rate
        om_rate <- sum(pred_occ@suitability == 0) / length(pred_occ@suitability)

        ## Prevalence
        xy_all <- data[[1]][, c(longitude, latitude)]
        all_data <- data.frame(xy_all, na.omit(raster::extract(varss, xy_all)))

        ellip <- ellipsoid_fit(all_data, longitude, latitude, methods[y], level)

        pred <- predict(ellip, projection_variables = var_vals, truncate = FALSE)

        #par <- paste0("Method_", methods[y], "_", names(variables[[2]])[x])
        par <- data.frame(Method = methods[y],
                          Variable_set = names(variables[[2]])[x])
        res <- c(proc[[1]], Omission_rate = om_rate, pred@prevalence[1],
                 pred@prevalence[2])
        nams <- names(res)
        res <- matrix(res, nrow = 1)
        colnames(res) <- nams
        res <- cbind.data.frame(par, res)

        ## Counting
        num <- ifelse(x == 1, 0, (length(methods) * (x - 1))) + y
        cat("\tParameter setting", num, "of", length(settings), "\n")
        return(res)
      })

      return(resul)
    })
  } else {
    cat("\nParallel option is still in development, sequential process will be used.\n")
    calibration <- sapply(1:length(variables[[2]]), function(x) {
      varss <- variables[[1]][[variables[[2]][[x]]]]
      var_vals <- na.omit(raster::values(varss))

      resul <- lapply(1:length(methods), function(y) {
        ## pROC
        xy_train <- data[[2]][, c(longitude, latitude)]
        train_data <- data.frame(xy_train, na.omit(raster::extract(varss, xy_train)))
        occ_vals <- na.omit(raster::extract(varss, data[[3]][, c(longitude, latitude)]))

        ellip <- ellipsoid_fit(train_data, longitude, latitude, methods[y], level)

        pred <- predict(ellip, projection_variables = var_vals, truncate = FALSE)
        pred_occ <- predict(ellip, projection_variables = occ_vals)

        proc <- partial_roc(pred@suitability, pred_occ@suitability, longitude,
                            latitude, error, iterations, percentage)

        ## Omission rate
        om_rate <- sum(pred_occ@suitability == 0) / length(pred_occ@suitability)

        ## Prevalence
        xy_all <- data[[1]][, c(longitude, latitude)]
        all_data <- data.frame(xy_all, na.omit(raster::extract(varss, xy_all)))

        ellip <- ellipsoid_fit(all_data, longitude, latitude, methods[y], level)

        pred <- predict(ellip, projection_variables = var_vals, truncate = FALSE)

        #par <- paste0("Method_", methods[y], "_", names(variables[[2]])[x])
        par <- data.frame(Method = methods[y],
                          Variable_set = names(variables[[2]])[x])
        res <- c(proc[[1]], Omission_rate = om_rate, pred@prevalence[1],
                 pred@prevalence[2])
        nams <- names(res)
        res <- matrix(res, nrow = 1)
        colnames(res) <- nams
        res <- cbind.data.frame(par, res)

        ## Counting
        num <- ifelse(x == 1, 0, (length(methods) * (x - 1))) + y
        cat("\tParameter setting", num, "of", length(settings), "\n")
        return(res)
      })

      return(resul)
    })
  }

  #calibration <- as.data.frame(do.call(rbind, calibration), stringsAsFactors = FALSE)
  calibration <- do.call(rbind, calibration)
  write.csv(calibration, paste0(output_directory, "/calibration_results.csv"),
            row.names = FALSE)

  # -----------
  # parameter selection
  cat("\nSelecting best parameter settings:\n")
  selected <- select_best(calibration, selection_criteria, level, error)

  write.csv(selected, paste0(output_directory, "/selected_parameterizations.csv"),
            row.names = FALSE)

  # -----------
  # producing report
  #cat("\nAnalyses finished. Producing HTML report...\n")
  #save(data, variable_names, variable1, n_var, r_values, ell_meta, mean_pred,
  #     layer, prevalences, replicates, replicate_type, bootstrap_percentage,
  #     color_palette, file = paste0(output_directory, "/enm_report_data.RData"))

  #report_format(name = paste0(output_directory, "/report_format"))
  #report(report_type = "calibration", output_directory)

  # -----------
  # returning results
  results <- calibration_ellipsoid(methods = methods,
                                   data = data,
                                   variable_sets = variables,
                                   level = level,
                                   results = calibration,
                                   selection_criteria = selection_criteria,
                                   selected_parameters = selected)

  return(results)
}
