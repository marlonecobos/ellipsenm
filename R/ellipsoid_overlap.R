ellipsoid_overlap <- function(..., overlap_type = "all",
                              n_points = 1000000) {
  # -----------
  # detecting potential errors
  if (missing(...)) {
    stop("Argument ... is necessary to perform the analysis")
  } else {
    ells <- list(...)
    if (length(ells) < 2) {
      stop("ellipsoid* objects to be compare must be two or more.")
    }
    cls <- sapply(ells, function (x) {class(x)[1]})
    if (any(!cls %in% c("ellipsoid", "ellipsoid_model_sim", "ellipsoid_model_rep"))) {
      stop("All objects to be compared must be of class ellipsoid*.")
    }
  }
  if (!overlap_type[1] %in% c("all", "full", "back_union")) {
    stop("Argument overlap_type is not valid, see function's help.")
  }

  # -----------
  # preparing data
  ## ellipsoids
  ellipsoids <- lapply(list(...), function(x) {
    if (class(x)[1] %in% c("ellipsoid", "ellipsoid_model_sim")) {
      return(x)
    } else {
      ell <- slot(x, "ellipsoids")[["mean_ellipsoid"]]
      ell <- new("ellipsoid_model_sim",
                 method =  slot(ell, "method"),
                 level = slot(ell, "level"),
                 centroid = slot(ell, "centroid"),
                 covariance_matrix = slot(ell, "covariance_matrix"),
                 niche_volume = slot(ell, "niche_volume"),
                 semi_axes_length = slot(ell, "semi_axes_length"),
                 axes_coordinates = slot(ell, "axes_coordinates"))
      slot(ell, "suitability", check = FALSE) <- slot(x, "suitability")[, "mean_ellipsoid"]
      slot(ell, "mahalanobis", check = FALSE) <- slot(x, "mahalanobis")[, "mean_ellipsoid"]
      return(ell)
    }
  })

  ## data for niche comparisons
  n_niches <- length(ellipsoids)
  comparison_matrix <- combn(n_niches, 2)

  if (!overlap_type[1] %in% c("all", "full")) {
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

    metrics <- overlap_metrics(comparison_matrix, env_data = data_rand,
                               maha_distM, suits)

    over_metrics <- lapply(1:length(metrics), function(x) {metrics[[x]][[1]]})
    over_metrics <- do.call(rbind.data.frame, over_metrics)
    rownames(over_metrics) <- sapply(1:ncol(comparison_matrix), function(x) {
      paste0("Niche_", paste0(comparison_matrix[,x], collapse = "_vs_"))
    })

    over_cordinates <- lapply(1:length(metrics), function(x) {metrics[[x]][[2]]})
    names(over_cordinates) <- rownames(over_metrics)

    results <-  list(overlap_results=nicheover_metrics,
                     overlap_coordinates= nicheover_cordinates,
                     ellipsoid_metadata=ellipsoid_metadata)

    results <- overlap_ellipsoid(ellipsoids = ellipsoids,
                                 spp_data = "list",
                                 background = "list",
                                 overlap = "data.frame")
  }


  # -----------
  # Background union overlap
  if(!is.null(bg_vars) && class(bg_vars) == "RasterStack"){
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
