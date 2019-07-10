#' Helper function to project models in ellipsoid model
#'
#' @param ellipsoid a fitted object of class ellipsoid*.
#' @param projection_variables (RasterStack, list, or character): if RasterStack,
#' a stack of layers respresenting an only scenario for projection; if list, a
#' named list of RasterStacks representing multiple scenarios for projection; if
#' character, name of the folder (in the working directory) containing other
#' folders (scenarios for projection) with raster layers to be used as variables.
#' See details of \code{\link{ellipsoid_model}}.
#' @param prvariables_format (character) if \code{projection_variables} is a list,
#' raster type of variables (raster layers) to be used and located in
#' subdirectories. Default = NULL. See \code{\link[raster]{writeFormats}} for
#' details and options.
#' @param sp_name (character) name of the species for which model(s) will be
#' projected. If not defined "species" is used.
#' @param prediction (character) type of prediction to be made, options are:
#' "suitability", "mahalanobis", and "both". Default = "suitability".
#' @param return_numeric (logical) whether or not to return values of mahalanobis
#' distance and suitability as part of the results (it depends on the type of
#' \code{prediction} selected). If \code{projection_variables} is a RasterStack,
#' default = FALSE; if \code{projection_variables} is a matrix, default = TRUE.
#' For both options the default can be changed. See details.
#' @param tolerance the tolerance for detecting linear dependencies.
#' Default = 1e-60.
#' @param format (charater) raster type of layers to be written.
#' See \code{\link[raster]{writeFormats}} for details and options.
#' @param overwrite (logical) whether or not to overwrite an exitent file with
#' the exact same name. Default = FALSE.
#' @param force_return (logical) whether or not to force returning numeric and
#' raster results for one of the ellipsoids defined by name in \code{return_name}.
#' See details.
#' @param return_name (character) names of the ellipsoid (part of \code{ellipsoid})
#' for which numeric and raster results will be forced to return. Default = NULL.
#' @param output_directory name of the folder were all results will be written.
#' This avoids saturation of the RAM.
#'
#' @return
#' A list with an example of each of the following: raster layers of suitability
#' or mahalanobis distances, depending on the type of \code{prediction} selected;
#' and the prevalence in environemnetal and geographic space for all projection
#' scenarios, if \code{prediction} was other than mahalanobis.
#'
#' Additionally all results of projections will be written in the folder defined
#' in \code{output_directory}.
#'
#' @export

model_projection <- function(ellipsoid, projection_variables, prvariables_format = NULL,
                             sp_name, prediction = "suitability", return_numeric = TRUE,
                             tolerance = 1e-60, format = "GTiff", overwrite = FALSE,
                             force_return = FALSE, return_name = NULL,
                             output_directory = "ellipsenm_projections") {
  # -----------
  # detecting potential errors
  if (!missing(ellipsoid)) {
    cls <- class(ellipsoid)[1]
    if (!cls %in% c("ellipsoid", "ellipsoid_model_sim", "ellipsoid_model_rep")) {
      stop("Argument ellipsoid must be of class ellipsoid*.")
    }
  } else {
    stop("Argument ellipsoid is necessary to perform the analysis.")
  }
  if (missing(projection_variables)) {
    stop("Argument projection_variables needs to be defined.")
  }
  if (class(projection_variables)[1] == "character" & is.null(prvariables_format)) {
    stop("Argument prvariables_format needs to be defined when projection_variables is a character.")
  }
  if(missing(sp_name)) {
    sp_name <- "species"
    warning("\nsp_name not defined, species will be used as generic sp_name.\n")
  }

  # -----------
  # producing projections
  cclas <- class(projection_variables)[1]
  if (cclas == "character" | cclas == "RasterStack" | cclas == "list") {
    # preparing data
    if (cls %in% c("ellipsoid", "ellipsoid_model_sim")) {
      ellv_names <- names(ellipsoid@centroid)
    } else {
      ellv_names <- names(ellipsoid@ellipsoids[[1]]@centroid)
    }
    nam_format <- rformat_type(format)

    if (cclas == "RasterStack") {
      r_values <- na.omit(projection_variables)
      lnames <- "projection"
      if (all(ellv_names == names(projection_variables))) {
        namer <- paste0(output_directory, "/", lnames, "_", sp_name, nam_format)

        predictions <- predict(ellipsoid, projection_variables, prediction, return_numeric,
                               tolerance, namer, format, overwrite)
      } else {
        stop("Variable names of projection_variables do not match variable names in ellipsoid.")
      }
    }

    if (cclas == "list") {
      # preparing further data
      lnames <- names(projection_variables)
      if (is.null(lnames)) {
        lnames <- paste0("scenario", 1:length(projection_variables))
      }
      predictions <- list()

      for (i in 1:length(lnames)) {
        cat("   Projection to", lnames[i], "\n")
        if (all(ellv_names == names(projection_variables[[i]]))) {
          if (i == 1) {r_values <- na.omit(projection_variables[[i]])}
          namer <- paste0(output_directory, "/", lnames, "_", sp_name, nam_format)

          predictions[[i]] <- predict(ellipsoid, projection_variables[[i]], prediction,
                                      return_numeric, tolerance, namer, format, overwrite)
        } else {
          stop("Variable names of projection_variables do not match variable names in ellipsoid.")
        }
      }
    }

    if (cclas == "character") {
      # preparing further data
      dirs <- dir(projection_variables, full.names = TRUE)
      lnames <- dir(projection_variables)
      format_pl <- paste0(rformat_type(prvariables_format), "$")

      predictions <- list()
      for (i in 1:length(dirs)) {
        cat("   Projection to", lnames[i], "\n")
        namest <- gsub(format_pl, "", list.files(dirs[i], pattern = format_pl))

        if (all(ellv_names == namest)) {
          p_layers <- raster::stack(list.files(dirs[i], pattern = format_pl,
                                               full.names = TRUE))
          if (i == 1) {r_values <- na.omit(p_layers)}

          namer <- paste0(output_directory, "/", lnames[i], "_", sp_name, nam_format)

          if (cls == "ellipsoid_model_rep") {
            predictions[[i]] <- predict(ellipsoid, p_layers, prediction, return_numeric,
                                        tolerance, namer, format, overwrite, force_return,
                                        return_name)
          } else {
            predictions[[i]] <- predict(ellipsoid, p_layers, prediction, return_numeric,
                                        tolerance, namer, format, overwrite, force_return)
          }

        } else {
          stop("Variable names of projection_variables do not match variable names in ellipsoid.")
        }
      }
    }
  } else {
    stop("Argument projection_variables is not valid, see function's help.")
  }

  # -----------
  # returning metadata for projection scenarios and returning results
  cat("\nWriting prevalence for projection scenario(s):\n")
  if (cclas == "RasterStack") {
    if (prediction != "mahalanobis") {
      suit <- predictions@prediction_suit
      prevalence <- predictions@prevalence
      if (prediction == "both") {
        maha <- predictions@prediction_maha
      }
      write.csv(prevalence, paste0(output_directory, "/projection_prevalence.csv"),
                row.names = TRUE)
      cat("\tPrevalence in projection scenario\n")
    } else {
      maha <- predictions@prediction_maha
    }
  } else {
    if (prediction != "mahalanobis") {
      suit <- predictions[[1]]@prediction_suit
      prevalence <- predictions[[1]]@prevalence
      if (prediction == "both") {
        maha <- predictions[[1]]@prediction_maha
      }
      prevs <- lapply(1:length(predictions), function(x) {
        write.csv(predictions[[x]]@prevalence,
                  paste0(output_directory, "/", lnames[x],"_prevalence.csv"),
                  row.names = TRUE)
        cat("\tPrevalence in scenario", lnames[x],"\n")
      })
    } else {
      maha <- predictions[[1]]@prediction_maha
    }
  }

  # -----------
  # returning results
  nb <- nrow(r_values)
  n_prop <- ifelse(nb > 100000, 0.1, 0.3)
  set.seed(1)
  samp <- sample(nb, ceiling(nb * n_prop))
  r_values <- r_values[samp, ]


  # -----------
  # returning results
  if (prediction != "mahalanobis") {
    if (prediction == "both") {
      results <- list(s_layer = suit, m_layer = maha, prevalence = prevalence,
                      r_values = r_values, scenarios = lnames)
    } else {
      results <- list(s_layer = suit, prevalence = prevalence,
                      r_values = r_values, scenarios = lnames)
    }
  } else {
    results <- list(m_layer = maha, prevalence = vector(),
                    r_values = r_values, scenarios = lnames)
  }

  return(results)
}
