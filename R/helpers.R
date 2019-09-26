#' Helper function to find raster extention
#' @param format (character) any of the format types allowed for raster objects.
#' See \code{\link[raster]{writeFormats}}
#' @export
#' @return Raster extension according to format type.

rformat_type <- function(format) {
  if (missing(format)) {stop("Argument format needs to be defined.")}
  if (format == "raster") {format1 <- ".grd"}
  if (format == "GTiff") {format1 <- ".tif"}
  if (format == "EHdr") {format1 <- ".bil"}
  if (format == "ascii") {format1 <- ".asc"}
  if (format == "SAGA") {format1 <- ".sdat"}
  if (format == "IDRISI") {format1 <- ".rst"}
  if (format == "CDF") {format1 <- ".nc"}
  if (format == "ENVI") {format1 <- ".envi"}
  if (format == "HFA") {format1 <- ".img"}
  return(format1)
}

#' Helperf function to get and write ellipsoid metadata
#' @param ellipsoid object of class ellipsoid*.
#' @param name (character) name of the file to be written. If the object in
#' \code{ellipsoid} is replicated, names present in slot ellipsoids are added as
#' prefixes to each file. Default = "ellipsoid_metadata".
#' @export
#' @return
#' A summary of ellipsoid metadata as a data.frame.
#'
#' Writes a csv file with the summary of ellipsoid metadata named
#' "metadata_summary.csv" and txt files with the complete ellipsoid
#' metadata (per each element if replicated) in the working directory.

write_ellmeta <- function(ellipsoid, name = "ellipsoid_metadata") {
  if (!missing(ellipsoid)) {
    cls <- class(ellipsoid)[1]
    if (!cls %in% c("ellipsoid", "ellipsoid_model_sim", "ellipsoid_model_rep")) {
      stop("Argument ellipsoid must be of class ellipsoid*.")
    }
  } else {
    stop("Argument ellipsoid is necessary to perform the analysis.")
  }
  name <- gsub("\\\\", "/", name)
  name <- unlist(strsplit(name, "/"))
  ndir <- paste0(paste(name[-length(name)], collapse = "/"), "/")
  namesum <- paste0(ndir, "metadata_summary.csv")
  name <- paste0(name[length(name)], ".txt")

  if (cls %in% c("ellipsoid", "ellipsoid_model_sim")) {
    namesim <- name
    name <- paste0(ndir, name)
    cat("Ellipsoid_metadata\n", file = name, append = TRUE, sep = "")
    cat("\nMethod:\t", ellipsoid@method, file = name, append = TRUE, sep = "")
    cat("\n\nLevel:\t", ellipsoid@level, file = name, append = TRUE, sep = "")
    cat("\n\nCentroid:\n", paste0(names(ellipsoid@centroid), "\t",
                                ellipsoid@centroid, "\n"),
        file = name, append = TRUE, sep = "")
    cat("\nCovariance_matrix:\n",
        paste0(c("", colnames(ellipsoid@covariance_matrix)), collapse = "\t"),
        "\n", file = name, append = TRUE, sep = "")
    suppressWarnings(write.table(ellipsoid@covariance_matrix, sep = "\t", file = name,
                                 append = TRUE, quote = FALSE, col.names = FALSE))
    cat("\nVolume:\t", ellipsoid@niche_volume, file = name, append = TRUE)
    cat("\n\nSemi-axes_length:\n", paste0(names(ellipsoid@semi_axes_length), "\t",
                                  ellipsoid@semi_axes_length, "\n"),
        file = name, append = TRUE, sep = "")
    cat("\nAxes_coordinates:\n", file = name, append = TRUE, sep = "")
    a_cord <- ellipsoid@axes_coordinates
    cords <- lapply(1:length(a_cord), function(x) {
      cat(letters[x], "\n", file = name, append = TRUE, sep = "")
      cat(paste0(c("", colnames(a_cord[[x]])), collapse = "\t"),
          "\n", file = name, append = TRUE, sep = "")
      suppressWarnings(write.table(a_cord[[x]], sep = "\t", file = name, append = TRUE,
                                   quote = FALSE, col.names = FALSE))
    })

    ell_meta <- data.frame(ellipsoid_model = c(ellipsoid@method, ellipsoid@level,
                                               round(ellipsoid@niche_volume, 2),
                                               namesim))

  } else {
    ellipsoid <- ellipsoid@ellipsoids
    nam <- names(ellipsoid)
    if (is.null(nam)) {
      enames <- as.character(1:length(ellipsoid))
      enames1 <- paste0("ellipsoid", enames)
    } else {
      if (length(grep("replicate", nam)) > 0 & length(grep("mean", nam)) > 0) {
        enames <- c(1:(length(nam) - 3), "mean", "min", "max")
        enames1 <- nam
      } else {
        enames <- nam
        enames1 <- nam
      }
    }
    namesim <- paste0(enames, "_", name)
    name <- paste0(ndir, enames, "_", name)
    ell_meta <- lapply(1:length(name), function(x) {
      cat("Ellipsoid_metadata\n", file = name[x], append = TRUE, sep = "")
      cat("\nMethod:\t", ellipsoid[[x]]@method, file = name[x], append = TRUE, sep = "")
      cat("\n\nLevel:\t", ellipsoid[[x]]@level, file = name[x], append = TRUE, sep = "")
      cat("\n\nCentroid:\n", paste0(names(ellipsoid[[x]]@centroid), "\t",
                                    ellipsoid[[x]]@centroid, "\n"),
          file = name[x], append = TRUE, sep = "")
      cat("\nCovariance_matrix:\n",
          paste0(c("", colnames(ellipsoid[[x]]@covariance_matrix)), collapse = "\t"),
          "\n", file = name[x], append = TRUE, sep = "")
      suppressWarnings(write.table(ellipsoid[[x]]@covariance_matrix, sep = "\t", file = name[x],
                                   append = TRUE, quote = FALSE, col.names = FALSE))
      cat("\nVolume:\t", ellipsoid[[x]]@niche_volume, file = name[x], append = TRUE)
      cat("\n\nSemi-axes_length:\n", paste0(names(ellipsoid[[x]]@semi_axes_length), "\t",
                                            ellipsoid[[x]]@semi_axes_length, "\n"),
          file = name[x], append = TRUE, sep = "")
      cat("\nAxes_coordinates:\n", file = name[x], append = TRUE, sep = "")
      a_cord <- ellipsoid[[x]]@axes_coordinates
      cords <- lapply(1:length(a_cord), function(y) {
        cat(letters[y], "\n", file = name[x], append = TRUE, sep = "")
        cat(paste0(c("", colnames(a_cord[[y]])), collapse = "\t"),
            "\n", file = name[x], append = TRUE, sep = "")
        suppressWarnings(write.table(a_cord[[y]], sep = "\t", file = name[x], append = TRUE,
                                     quote = FALSE, col.names = FALSE))
      })
      ellm <- c(ellipsoid[[x]]@method, ellipsoid[[x]]@level,
                round(ellipsoid[[x]]@niche_volume, digits = 2))
      return(ellm)
    })
    ell_meta <- as.data.frame(rbind(do.call(cbind, ell_meta), namesim))
    colnames(ell_meta) <- enames1
  }
  row.names(ell_meta) <- c("Method", "Level", "Volume", "Other_metadata")
  write.csv(ell_meta, namesum, row.names = TRUE)
  return(ell_meta)
}



#' Helper funtion to select best parameter settings
#' @param calibration_table data.frame of results from model calibration in
#' ellipsenm.
#' @param selection_criteria (character) set of criteria to select best models,
#' options are: "S_OR" (statistical significance and low omission) and
#' "S_OR_P" (statistical significance, low omission, and low prevalence).
#' See details. Default = "Sig_OR_Prev".
#' @param level (numeric) the confidence level of a pairwise confidence region
#' for the ellipsoid, expresed as percentage. Default = 95.
#' @param error (numeric) value from 0 to 100 to represent the percentage of
#' potential error (E) that the data could have due to any source of uncertainty.
#' Default = 5.
#' @export
#' @return
#' data.frame with the selected parameter settings according to the argument
#' \code{selection_criteria}.

select_best <- function(calibration_table, selection_criteria = "S_OR_P",
                        level = 95, error = 5) {
  if (selection_criteria %in% c("S_OR", "S_OR_P")) {
    sig <- calibration_table[calibration_table[, 4] <= error / 100, ]
    if (nrow(sig) == 0) {
      sig <- calibration_table[calibration_table[, 4] ==
                                 min(calibration_table[, 4]), ]
      warning("None of the parameter settings resulted in significant models.\n  The ones with the lowest partial ROC values were selected.\n")
    }
    res <- sig[sig[, 6] <= ((100 - level) / 100), ]
    if (nrow(res) == 0) {
      res <- sig[sig[, 6] == min(sig[, 6]), ]
      warning("None of the models had omission rates lower or equal than expected.\n  The ones with the lowest omission rates were selected.\n")
    }
    if (selection_criteria == "S_OR_P") {
      res <- res[res[, 8] == min(res[, 8]), ]
    }
  } else {
    stop("Argument selection_criteria is not valid, see function's help.")
  }
  cat("\tA total of", nrow(res), "paramter settings were selected.\n")
  return(res)
}

#' Helper funtion to create data_overlap objects
#' @param data data.frame of species' occurrence records. Columns must include
#' species, longitude, and latitude.  Optionally, if \code{variables} is a matrix
#' or data.frame, \code{data} must include more columns containing the values of
#' at least two variables to be used for fitting ellipsoid* models.
#' @param method (character) method to construct the ellipsoid that characterizes
#' the species ecological niche. Available methods are: "covmat", "mve1", and
#' "mve2". See details. Default = "covmat".
#' @param level (numeric) the confidence level of a pairwise confidence region
#' for the ellipsoid, expresed as percentage. Default = 95.
#' @param variables (optional) RasterStack, matrix, or data.frame of at least two
#' variables to represent a set of conditions relevant for overlap analyses.
#' @export
#' @return
#' A data_overlap object for posterior use in overlap analyses.

overlap_object <- function(data, species, longitude, latitude, method = "covmat",
                           level = 95, variables = NULL) {
  if (missing(data)) {
    stop("Argument data is needed for creating data_overlap object.")
  }
  if (!is.null(variables)) {
    if (!class(variables)[1] %in% c("RasterStack", "matrix", "data.frame")) {
      stop("Argument variables not valid, see function's help.")
    }
  } else {
    if (ncol(data) <= 3) {
      stop("If variables are not defined, data must include other columns representing environmental data.")
    }
  }
  object <- data_overlap(data = data,
                         main_columns = c(species, longitude, latitude),
                         method = method,
                         level = level)
  slot(object, "variables", check = FALSE) <- variables
  return(object)
}

#' Helper funtion to get attributes from ellipsoid lists
#' @param ellipsoids list of ellipsoid objects.
#' @param attribute (character) name of the attribute to be obtained from elements
#' in \code{ellipsoids}. Options are: method, centroid, covariance_matrix, level,
#' niche_volume, semi_axes_length, and axes_coordinates. Default = "method".
#' @export
#' @return
#' The attribute selected.

get_attribute <- function(ellipsoids, attribute = "method"){
  if (missing(ellipsoids)) {
    stop("Argument ellipsoids is needed, see function's help.")
  }
  if (!attribute %in% c("method", "centroid", "covariance_matrix", "level",
                        "niche_volume", "semi_axes_length", "axes_coordinates")) {
    stop("Invalid attribute, see function's help.")
  }
  attr1 <- lapply(ellipsoids, function(x) {
    attr1 <- slot(x, attribute)
    if(attribute == "axes_coordinates") {attr1 <- do.call(rbind, attr1)}
    return(attr1)
  })
  return(attr1)
}


#' Helper funtion to get data for Montecarlo simulations
#' @param ellipsoids list of ellipsoid objects.
#' @param n_points (numeric) number of random points to be generated.
#' @export
#' @return
#' A total of \code{n_points} created randomly considering the limits of the
#' ellipsoids considered.

hypercube_boundaries <- function(ellipsoids, n_points = 1000000) {
  if (missing(ellipsoids)) {
    stop("Argument ellipsoids is needed, see function's help.")
  }
  axis_list <- get_attribute(ellipsoids, attribute = "axes_coordinates")
  ellipsoid_axis <- do.call(rbind, axis_list)
  min_env <- apply(ellipsoid_axis, 2, min)
  max_env <- apply(ellipsoid_axis, 2, max)
  data_rand <- sapply(1:length(min_env), function (i) {
    runif(n_points, min = min_env[i], max = max_env[i])
  })
  return(data_rand)
}


#' Helper function to compute overlap metrics
#' @param comparison_matrix matrix resulted from using the \code{\link[utils]{combn}}
#' function to obtain the pairwise overlap of all niches to be compared.
#' @param background matrix of environmental variables to be used as background.
#' @param mahalanobis matrix of mahalanobis distances calculated for all niches
#' to be compared.
#' @param suitability matrix of suitability values calculated for all niches to
#' be compared.
#' @export
#' @return
#' A list of two data.frames containing metrics of overlap and other data useful
#' for plotting results.

overlap_metrics <- function(comparison_matrix, background, mahalanobis,
                            suitability) {
  if (missing(comparison_matrix)) {
    stop("Argument comparison_matrix is necessary to perform the analysis")
  }
  if (missing(background)) {
    stop("Argument background is necessary to perform the analysis")
  }
  if (missing(mahalanobis)) {
    stop("Argument mahalanobis is necessary to perform the analysis")
  }
  if (missing(suitability)) {
    stop("Argument suitability is necessary to perform the analysis")
  }
  get_metrics <- lapply(1:ncol(comparison_matrix), function(x) {
    compare <- suitability[, comparison_matrix[, x]]
    compare <- compare > 0
    zeros1 <- which(rowSums(compare) %in% 0)
    compare <- cbind(1:nrow(compare), compare)
    compare <- compare[-zeros1, ]
    union_npoints <- nrow(compare)
    union_prop <- colSums(compare[, -1]) / union_npoints
    union_sum <- rowSums(compare[, -1])
    intersect_id <- which(union_sum == 2)
    np_intesection_global <- length(intersect_id)
    intersection <- np_intesection_global / union_npoints
    prop_size_niche_1_vs_2 <- union_prop[1] / union_prop[2]
    prop_size_niche_2_vs_1 <- union_prop[2] / union_prop[1]

    overlap <- data.frame(total_points = union_npoints,
                          overlaped_points = np_intesection_global,
                          overlap = intersection,
                          prop_size_niche_1_vs_2 = prop_size_niche_1_vs_2,
                          prop_size_niche_2_vs_1 = prop_size_niche_2_vs_1)

    background_in <- data.frame(background[compare[intersect_id, 1], ],
                                mahalanobis[compare[intersect_id, 1],
                                            comparison_matrix[, x]],
                                suitability[compare[intersect_id, 1],
                                            comparison_matrix[, x]])

    return(list(overlap = overlap, background = background_in))
  })
  return(get_metrics)
}
