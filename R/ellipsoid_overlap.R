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
#' performing Monte-Carlo simulations (full overlap). Default = 1000000.
#'
#' @return
#' A overlap_ellipsoid object containing all results from overlap analyses as
#' well as other information needed for plotting.
#'
#' @details
#' Types of overlap are as follows:
#'
#' all = performs all types of overlap analyses allowed.
#'
#' full = measures overlap of the complete volume of the ellipsoidal niches.
#'
#' back_union = meausures overlap of ellipsoidal niches considering only the
#' union of the environmental conditions relevant for the two species (backgrounds).
#'
#' @export
#'
#' @examples
#' # data
#'
#' # full overlap analysis
#'
#' # overlap measured in the union of the two species background

ellipsoid_overlap <- function(..., overlap_type = "all",
                              n_points = 1000000) {
  # -----------
  # detecting potential errors
  if (missing(...)) {
    stop("Argument ... is necessary to perform the analysis")
  } else {
    plits <- list(...)
    if (length(plits) < 2) {
      stop("ellipsoid* objects to be compare must be two or more.")
    }
    cls <- sapply(plits, function (x) {class(x)[1]})
    if (any(cls != "data_overlap")) {
      stop("All objects to be compared must be of class ellipsoid*.")
    }
  }
  if (!overlap_type[1] %in% c("all", "full", "back_union")) {
    stop("Argument overlap_type is not valid, see function's help.")
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
    message("overlap_type using background is not possible, only full overlap will be performed.")
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
    stop("Variables from all data_overlap objects must be of the same class.")
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
    ellipsoid_fit(data[[x]], longitude[x], latitude[x], method[x], level[x],
                  variables[[x]])
  })

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
      bg_vars <- lapply(variables, function(x) {na.omit(raster::values(x))})
      bg_vars <- unique(do.call(rbind, bg_vars))
    } else {
      bg_vars <- unique(do.call(rbind, variables))
    }
    colnames(bg_vars) <- variable_names

    mah_suit <- sapply(1:length(ellipsoids), function(x) {
      mh_st <- predict(ellipsoids[[x]], bg_vars, "both")
      mh_st <- data.frame(slot(mh_st, "suitability"), slot(mh_st, "mahalanobis"))
      return(mh_st)
    })

    suits <- mah_suit[1,]
    names(suits) <- paste0("Niche_", 1:length(ellipsoids), "_S")
    suits <-  do.call(cbind, suits)

    mahas <- mah_suit[2,]
    names(mahas) <- paste0("Niche_", 1:length(ellipsoids),"_M")
    mahas <-  do.call(cbind, mahas)

    metrics <- overlap_metrics(comparison_matrix, bg_vars, mahas, suits)

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
  cat("\nProcess finished\n")

  return(results)
}
