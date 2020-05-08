#' Overlap of ellipsoid-based ecological niche models
#'
#' @description ellipsoid_overlap performs analyses to measure the degree of
#' overlap between two or more ellipsoid-based ecological niches as in pairwise
#' comparisons. Measures can be done considering the entire ellipsoid volume or
#' sets of environmental conditions (background).
#'
#' @param ... data_overlap objects containing data for individual niches to be
#' compared in overlap analyses. At least two data_overlap objects are needed
#' to perform analyses. These objects can be created with the function
#' \code{\link{overlap_object}}.
#' @param overlap_type (character) type of overlap to be measured. Options are:
#' "all", "full", and "back_union". Default = "all". See details.
#' @param n_points (character) number of random points to be generated for
#' performing Monte-Carlo simulations for full overlap-type measurements.
#' Default = 1000000.
#' @param significance_test (logical) whether or not to perform a test to determine
#' statistical significance of overlap results. See details; default = FALSE.
#' @param replicates (numeric) number of replicates to be performed during the
#' significance test; default = 1000.
#' @param confidence_limit (numeric) confidence limit for the significance test.
#' Default = 0.05
#'
#' @return
#' An object of class \code{\link{overlap_ellipsoid}} containing all results
#' from overlap analyses as well as other information needed for plotting.
#'
#' @usage
#' ellipsoid_overlap(..., overlap_type = "all", n_points = 1000000,
#'                   significance_test = FALSE, replicates = 1000,
#'                   confidence_limit = 0.05)
#'
#' @details
#' Types of overlap are as follows:
#' - "all", performs all types of overlap analyses allowed.
#' - "full", measures overlap of the complete volume of the ellipsoidal niches.
#' - "back_union", meausures overlap of ellipsoidal niches considering only the
#' union of the environmental conditions relevant for the two species (backgrounds).
#'
#' The statistical significance test consist in randomly sampling the background
#' with n = to the number of records of each species and creting ellipsoids with
#' such data. Overlap is measured for each pair of random-ellipsoids according to
#' the \code{overlap_type} selected. The process is repeated \code{replicate}
#' times and the observed overlap value is compared to the values found for
#' pairs of random-ellipsoids. The null hypothesis is that the niches are
#' overlaped and if the observed values are as exterme or more extreme than the
#' lower limit of the values found for the random-ellipsoids, the null hypothesis
#' is rejected. This is, if the observed overlap value is lower than the 95%
#' (or the value defined in \code{confidence_limit}) of the random-derived values
#' of overlap, the niches are considered not-overlapped. If the observed value
#' cannot be distinguished from random, the null hypothesis cannot be rejected.
#' A p-value and the pre-defined \code{confidence_limit} will be added to the
#' overlap matrix when the test is performed. A list with all the overlap results
#' from the analyses with random-ellipsoids will be added to the
#' \code{\link{overlap_ellipsoid}} object returned.
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
#' # preparing data
#' vext <- raster::extent(vars)
#' ext1 <- raster::extent(vext[1], (mean(vext[1:2]) + 0.2), vext[3:4])
#' ext2 <- raster::extent((mean(vext[1:2]) + 0.2), vext[2], vext[3:4])
#'
#' # croping variables and splitting occurrences
#' vars1 <- raster::stack(raster::crop(vars, ext1))
#' vars2 <- raster::stack(raster::crop(vars, ext2))
#'
#' occurrences1 <- occurrences[occurrences$longitude < (mean(vext[1:2]) + 0.2), ]
#' occurrences2 <- occurrences[!occurrences$longitude %in% occurrences1$longitude, ]
#'
#' # preparing overlap objects to perform analyses
#' niche1 <- overlap_object(occurrences1, species =  "species", longitude = "longitude",
#'                          latitude = "latitude", method = "covmat", level = 95,
#'                          variables = vars1)
#'
#' niche2 <- overlap_object(occurrences2, species =  "species", longitude = "longitude",
#'                          latitude = "latitude", method = "covmat", level = 95,
#'                          variables = vars2)
#'
#' # niche overlap analysis
#' overlap <- ellipsoid_overlap(niche1, niche2)

ellipsoid_overlap <- function(..., overlap_type = "all", n_points = 1000000,
                              significance_test = FALSE, replicates = 1000,
                              confidence_limit = 0.05) {
  # -----------
  # detecting potential errors
  if (missing(...)) {
    stop("Argument '...' is necessary to perform the analysis")
  } else {
    plits <- list(...)
    if (length(plits) < 2) {
      stop("At least two ellipsoid* objects are needed to perfrom analyses.")
    }
    cls <- sapply(plits, function (x) {class(x)[1]})
    if (any(cls != "data_overlap")) {
      stop("All objects to be compared must be of class ellipsoid*.")
    }
  }
  if (!overlap_type[1] %in% c("all", "full", "back_union")) {
    stop("Argument 'overlap_type' is not valid, see function's help.")
  }

  # -----------
  # preparing data
  ## background availability
  cat("\nPreparing data...\n")
  back <- sapply(plits, function(x) {
    if (is.null(slot(x, "variables"))) {FALSE} else {TRUE}
  })
  back <- all(back == TRUE)

  if (back == FALSE & overlap_type[1] %in% c("all", "back_union")) {
    message("overlap analysis using background is not possible, only full overlap will be performed.")
    overlap_type <- "full"
  }

  ## arguments
  data <- lapply(plits, function(x) {slot(x, "data")})
  species <- sapply(plits, function(x) {slot(x, "main_columns")[1]})
  longitude <- sapply(plits, function(x) {slot(x, "main_columns")[2]})
  latitude <- sapply(plits, function(x) {slot(x, "main_columns")[3]})
  data <- lapply(1:length(data), function(x) {
    data[[x]][, colnames(data[[x]]) != species[x]]
  })
  method <- sapply(plits, function(x) {slot(x, "method")})
  level <- sapply(plits, function(x) {slot(x, "level")})
  variables <- lapply(plits, function(x) {
    if (back == TRUE) {slot(x, "variables")} else {NULL}
  })
  variable_names <- sapply(1:length(plits), function(x) {
    if (back == TRUE) {
      vars <- slot(plits[[x]], "variables")
      cl <- class(vars)[1]
      if (cl == "RasterStack") {v <- names(vars)} else {v <- colnames(vars)}
    } else {
      v <- colnames(data[[x]])
      v <- v[!v %in% c(longitude[x], latitude[x])]
    }
    return(v)
  })[, 1]
  cls <- sapply(variables, function (x) {class(x)[1]})
  if (back == TRUE & all(cls != cls[1])) {
    stop("Variables from all 'data_overlap' objects must be of the same class.")
  }
  data <- lapply(1:length(data), function(x) {
    if (class(variables[[x]])[1] == "RasterStack") {
      cbind(data[[x]], raster::extract(variables[[x]], data[[x]]))
    } else {
      data[[x]]
    }
  })

  ## ellipsoids
  cat("\nFitting ellipsoids\n")
  ellipsoids <- lapply(1:length(plits), function(x) {
    if (class(variables[[x]])[1] != "RasterStack") {
      ellipsoid_fit(data[[x]], longitude[x], latitude[x], method[x], level[x])
    } else {
      ellipsoid_fit(data[[x]], longitude[x], latitude[x], method[x], level[x],
                    variables[[x]])
    }
  })
  names(ellipsoids) <- paste0("Niche_", 1:length(ellipsoids))

  ## data for niche comparisons
  n_niches <- length(ellipsoids)
  comparison_matrix <- combn(n_niches, 2)

  if (overlap_type[1] %in% c("all", "full")) {
    data_rand <- hypercube_boundaries(ellipsoids, n_points = n_points)
    colnames(data_rand) <- variable_names
  }

  # -----------
  # Full overlap, Montocarlo simulation
  if (overlap_type[1] %in% c("all", "full")) {
    cat("\nPerforming full overlap analyses, please wait...\n")
    mah_suit <- sapply(1:length(ellipsoids), function(x) {
      mh_st <- predict(ellipsoids[[x]], data_rand, "both")
      mh_st <- data.frame(slot(mh_st, "suitability"), slot(mh_st, "mahalanobis"))
      return(mh_st)
    })

    suits <- mah_suit[1,]
    names(suits) <- paste0("Niche_", 1:length(ellipsoids), "_S")
    suits <-  do.call(cbind, suits)

    mahas <- mah_suit[2,]
    names(mahas) <- paste0("Niche_", 1:length(ellipsoids),"_M")
    mahas <-  do.call(cbind, mahas)

    metrics <- overlap_metrics(comparison_matrix, data_rand, mahas, suits)

    over_metrics <- lapply(1:length(metrics), function(x) {metrics[[x]][[1]]})
    over_metrics <- do.call(rbind.data.frame, over_metrics)
    rownames(over_metrics) <- sapply(1:ncol(comparison_matrix), function(x) {
      paste0("Niche_", paste0(comparison_matrix[,x], collapse = "_vs_"))
    })

    over_cordinates <- lapply(1:length(metrics), function(x) {metrics[[x]][[2]]})
    names(over_cordinates) <- rownames(over_metrics)

    results <- overlap_ellipsoid(ellipsoids = ellipsoids,
                                 data = data,
                                 full_background = over_cordinates,
                                 full_overlap = over_metrics,
                                 variable_names = variable_names)
  } else {
    results <- overlap_ellipsoid(ellipsoids = ellipsoids,
                                 data = data,
                                 variable_names = variable_names)
  }

  # -----------
  # Background union overlap
  if(back == TRUE & overlap_type[1] %in% c("all", "back_union")) {
    cat("\nPerforming union background overlap analyses, please wait...\n")
    if (cls[1] == "RasterStack") {
      backg <- lapply(variables, function(x) {na.omit(raster::values(x))})
      background <- unique(do.call(rbind, backg))
    } else {
      backg <- variables
      background <- unique(do.call(rbind, variables))
    }
    colnames(background) <- variable_names

    mah_suit <- sapply(1:length(ellipsoids), function(x) {
      mh_st <- predict(ellipsoids[[x]], background, "both")
      mh_st <- data.frame(slot(mh_st, "suitability"), slot(mh_st, "mahalanobis"))
      return(mh_st)
    })

    suits <- mah_suit[1,]
    names(suits) <- paste0("Niche_", 1:length(ellipsoids), "_S")
    suits <-  do.call(cbind, suits)

    mahas <- mah_suit[2,]
    names(mahas) <- paste0("Niche_", 1:length(ellipsoids),"_M")
    mahas <-  do.call(cbind, mahas)

    metrics <- overlap_metrics(comparison_matrix, background, mahas, suits)

    over_metrics <- lapply(1:length(metrics), function(x) {metrics[[x]][[1]]})
    over_metrics <- do.call(rbind.data.frame, over_metrics)
    rownames(over_metrics) <- sapply(1:ncol(comparison_matrix), function(x) {
      paste0("Niche_", paste0(comparison_matrix[,x], collapse = "_vs_"))
    })

    over_cordinates <- lapply(1:length(metrics), function(x) {metrics[[x]][[2]]})
    names(over_cordinates) <- rownames(over_metrics)

    slot(results, "union_background") <- over_cordinates
    slot(results, "union_overlap") <- over_metrics
  }

  if (significance_test == TRUE) {
    cat("\nPerforming statistical significance analyses, please wait:\n")
    if (overlap_type == "all") {background <- list(data_rand, backg)}
    if (overlap_type == "full") {background <- list(data_rand)}
    if (overlap_type == "back_union") {background <- list(backg)}
    sample_size <- sapply(data, function(y) {nrow(y)})

    random_results <- overlap_random(background, sample_size, method, level,
                                     overlap_type, replicates)

    if (overlap_type[1] %in% c("all", "full")) {
      f_pval <- sapply(1:length(random_results[[1]]), function(y) {
        sum(random_results[[1]][[y]][, 3] <= results@full_overlap[y, 3]) / replicates
      })
      u_pval <- sapply(1:length(random_results[[2]]), function(y) {
        sum(random_results[[2]][[y]][, 3] <= results@union_overlap[y, 3]) / replicates
      })
      f_metrics_s <- data.frame(total_points = results@full_overlap[, 1],
                                overlapped_points = results@full_overlap[, 2],
                                overlap = results@full_overlap[, 3],
                                p_value = f_pval,
                                pre_defined_CL = confidence_limit,
                                prop_size_niche_1_vs_2 = results@full_overlap[, 4],
                                prop_size_niche_2_vs_1 = results@full_overlap[, 5])
      rownames(f_metrics_s) <- rownames(results@full_overlap)

      b_metrics_s <- data.frame(total_points = results@union_overlap[, 1],
                                overlapped_points = results@union_overlap[, 2],
                                overlap = results@union_overlap[, 3],
                                p_value = u_pval,
                                pre_defined_CL = confidence_limit,
                                prop_size_niche_1_vs_2 = results@union_overlap[, 4],
                                prop_size_niche_2_vs_1 = results@union_overlap[, 5])
      rownames(b_metrics_s) <- rownames(results@union_overlap)

      slot(results, "full_overlap") <- f_metrics_s
      slot(results, "union_overlap") <- b_metrics_s
    }
    if (overlap_type[1] == "full") {
      f_pval <- sapply(1:length(random_results[[1]]), function(y) {
        sum(random_results[[1]][[y]][, 3] <= results@full_overlap[y, 3]) / replicates
      })
      f_metrics_s <- data.frame(total_points = results@full_overlap[, 1],
                                overlapped_points = results@full_overlap[, 2],
                                overlap = results@full_overlap[, 3],
                                p_value = f_pval,
                                pre_defined_CL = confidence_limit,
                                prop_size_niche_1_vs_2 = results@full_overlap[, 4],
                                prop_size_niche_2_vs_1 = results@full_overlap[, 5])
      rownames(f_metrics_s) <- rownames(results@full_overlap)
      slot(results, "full_overlap") <- f_metrics_s
    }
    if (back == TRUE & overlap_type[1] == "back_union") {
      u_pval <- sapply(1:length(random_results[[1]]), function(y) {
        sum(random_results[[1]][[y]][, 3] <= results@union_overlap[y, 3]) / replicates
      })
      b_metrics_s <- data.frame(total_points = results@union_overlap[, 1],
                                overlapped_points = results@union_overlap[, 2],
                                overlap = results@union_overlap[, 3],
                                p_value = u_pval,
                                pre_defined_CL = confidence_limit,
                                prop_size_niche_1_vs_2 = results@union_overlap[, 4],
                                prop_size_niche_2_vs_1 = results@union_overlap[, 5])
      rownames(b_metrics_s) <- rownames(results@union_overlap)
      slot(results, "union_overlap") <- b_metrics_s
    }
    slot(results, "significance_results") <- random_results
  }

  cat("\nProcess finished\n")

  return(results)
}
