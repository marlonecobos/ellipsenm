#' Overlap of ellipsoid-based ecological niche models
#'
#' @description
#'
#' @param ... data_overlap objects containing data for individual niches to be
#' compared in overlap analyses. At least two data_overlap objects are needed
#' to perform analyses. These objects can be created with the function
#' \code{\link{overlap_object}}.
#' @param overlap_type (character) type of overlap to be measured. Options are:
#' "all", "full", and "back_union". Default = "all". See details.
#' @param n_points
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
  back <- sapply(plits, function(x) {
    if (is.null(slot(x, "variables"))) {FALSE} else {TRUE}
  })
  back <- all(back == TRUE)

  if (back == FALSE & overlap_type[1] %in% c("all", "back_union")) {
    warning("overlap_type using background is not possible, only full overlap will be performed.")
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

  ## ellipsoids
  ellipsoids <- lapply(1:length(plits), function(x) {
    ell <- ellipsoid_fit(data[[x]], longitude[x], latitude[x], method[x],
                         level[x], variables[[x]])
    if (back == TRUE) {
      ell <- predict(ell, variables, "both", return_numeric = TRUE)
    }
    return(ell)
  })

  ## data for niche comparisons
  n_niches <- length(ellipsoids)
  comparison_matrix <- combn(n_niches, 2)

  if (overlap_type[1] %in% c("all", "full")) {
    data_rand <- hypercube_boundaries(ellipsoids, n_points = n_points)
  }

  # -----------
  # Full overlap, Montocarlo simulation
  if (!overlap_type[1] %in% c("all", "full")) {
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
                                 spp_data = data,
                                 full_background = over_cordinates,
                                 full_overlap = over_metrics)
  }


  # -----------
  # Background union overlap
  if(back == TRUE & !overlap_type[1] %in% c("all", "back_union")){
    cls <- sapply(variables, function (x) {class(x)[1]})
    a

    bg_vars <- na.omit(getValues(bg_vars))
    in_EllipsoidMh <- inEllipsoidMh(emd = ellipsoid_metadata,
                                    env_data = bg_vars,level)

    in_Ellipsoid <- in_EllipsoidMh[1,]
    names(in_Ellipsoid) <- paste0("Niche_",1:length(envdata_list))
    in_Ellipsoid <-  as.data.frame(in_Ellipsoid)
    in_Ellipsoid <- as.matrix(in_Ellipsoid)
    maha_distM <- in_EllipsoidMh[2,]
    names(maha_distM) <- paste0("Niche_",1:length(envdata_list),"_Mh")
    maha_distM <-  as.data.frame(maha_distM)
    maha_distM <- as.matrix(maha_distM)

    niche_metricBG <- overlap_metrics(comparison_matrix, bg_vars,
                                      maha_distM, in_Ellipsoid)

    nicheover_metricBG <- lapply(1:length(niche_metricBG), function(x){
      niche_metricBG[[x]][[1]]
    })

    nicheover_cordinatesBG <- lapply(1:length(niche_metricBG), function(x){
      niche_metricBG[[x]][[2]]
    })


    nicheover_metricBG <- do.call(rbind.data.frame,nicheover_metricBG)

    rownames(nicheover_metricBG) <- sapply(1:dim(comparison_matrix)[2],function(x)
      paste0("Niche_",paste0(comparison_matrix[,x],collapse="_vs_")))

    names(nicheover_cordinatesBG) <- rownames(nicheover_metricBG)

    overlap_results <- list(nicheover_metrics, nicheover_metricBG)
    names(overlap_results) <- c("Full_Overlap","Background_Context_Overlap")
    nicheover_cordinatesAll <- list(nicheover_cordinates,nicheover_cordinatesBG)
    names(nicheover_cordinatesAll) <- c("Full_Overlap","Background_Context_Overlap")
    results <-  list(overlap_results=overlap_results,
                     overlap_coordinates= nicheover_cordinatesAll,
                     ellipsoid_metadata=ellipsoid_metadata)

  }


  return(results)

}
