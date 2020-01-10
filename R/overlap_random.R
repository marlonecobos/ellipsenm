#' Overlap analyses for random sampled data
#'
#' @description overlap_random performs analyses of overlap on ellipsoidal objects
#' derived from random samples of their background.
#'
#' @param background list containing the backgrounds to be used. If \code{overlap_type}
#' = "all", background for "full" overlap must be first, and background for "union"
#' overlap, second.
#' @param sample_size (numeric) vector of sample sizes for the species under
#' analyses that should be sampled in each replicate.
#' @param method (character) vector of methods to construct the ellipsoids that
#' characterize the species ecological niches.
#' @param level (numeric) vector of the confidence levels of a pairwise
#' confidence region for the ellipsoids, expresed as percentage.
#' @param overlap_type (character) type of overlap to be measured. Options are:
#' "all", "full", and "back_union". Default = "all". See details.
#' @param replicates (numeric) number of replicates to be performed during the
#' significance test; default = 1000.
#'
#' @return
#' A list of data.frames with the overlap metrics calculated for each replicate,
#' according to the \code{overlap_type} selected.
#'
#' @usage
#' overlap_random(background, sample_size, method, level, overlap_type = "all",
#'                replicates = 1000)
#'
#' @export
#'
#' @details
#' Details about the \code{method} can be revised in \code{\link{ellipsoid_fit}}.
#'
#' Details about the overlap analyses can be checked in \code{\link{ellipsoid_overlap}}.

overlap_random <- function(background, sample_size, method, level,
                           overlap_type = "all", replicates = 1000) {
  # -----------
  # detecting potential errors
  if (missing(background)) {
    stop("Argument 'background' is necessary to perform the analysis")
  }
  if (missing(sample_size)) {
    stop("Argument 'sample_size' is necessary to perform the analysis")
  }
  if (missing(method)) {
    stop("Argument 'method' is necessary to perform the analysis")
  }
  if (missing(level)) {
    stop("Argument 'level' is necessary to perform the analysis")
  }

  if (overlap_type[1] %in% c("all", "full")) {nf <- 1; nb <- 2} else {nb <- 1}

  b_over <- lapply(1:replicates, function(x) {
    # -----------
    # Full overlap, Montocarlo simulation
    if (overlap_type[1] %in% c("all", "full")) {
      set.seed(x)
      sdata <- lapply(1:length(sample_size), function (y) {
        cbind(longitude = rep(0, sample_size[y]), latitude = rep(0, sample_size[y]),
              background[[nf]][sample(nrow(background[[nf]]), sample_size[y]), ])
      })

      ## ellipsoids
      ellipsoids <- lapply(1:length(sdata), function(y) {
        ellipsoid_fit(sdata[[y]], "longitude", "latitude", method[y], level[y])
      })
      names(ellipsoids) <- paste0("Niche_", 1:length(ellipsoids))

      ## data for niche comparisons
      n_niches <- length(ellipsoids)
      comparison_matrix <- combn(n_niches, 2)

      data_rand <- background[[nf]]

      ## predictions
      mah_suit <- sapply(1:length(ellipsoids), function(y) {
        mh_st <- predict(ellipsoids[[y]], data_rand, "both")
        mh_st <- data.frame(slot(mh_st, "suitability"), slot(mh_st, "mahalanobis"))
        return(mh_st)
      })

      suits <- mah_suit[1,]
      names(suits) <- paste0("Niche_", 1:length(ellipsoids), "_S")
      suits <-  do.call(cbind, suits)

      mahas <- mah_suit[2,]
      names(mahas) <- paste0("Niche_", 1:length(ellipsoids),"_M")
      mahas <-  do.call(cbind, mahas)

      metrics <- overlap_metrics(comparison_matrix, data_rand, mahas, suits,
                                 return_background = FALSE)

      over_metrics <- lapply(1:length(metrics), function(y) {metrics[[y]][[1]]})
      over_metrics <- do.call(rbind.data.frame, over_metrics)
      rownames(over_metrics) <- sapply(1:ncol(comparison_matrix), function(y) {
        paste0("Niche_", paste0(comparison_matrix[, y], collapse = "_vs_"))
      })
    }

    # -----------
    # Background union overlap
    if(overlap_type[1] %in% c("all", "back_union")) {
      set.seed(x)
      sdata <- lapply(1:length(sample_size), function (y) {
        cbind(longitude = rep(0, sample_size[y]), latitude = rep(0, sample_size[y]),
              background[[nb]][[y]][sample(nrow(background[[nb]][[y]]), sample_size[y]), ])
      })

      ## ellipsoids
      ellipsoids <- lapply(1:length(sdata), function(y) {
        ellipsoid_fit(sdata[[y]], "longitude", "latitude", method[y], level[y])
      })
      names(ellipsoids) <- paste0("Niche_", 1:length(ellipsoids))

      ## data for niche comparisons
      n_niches <- length(ellipsoids)
      comparison_matrix <- combn(n_niches, 2)

      back <- unique(do.call(rbind, background[[nb]]))

      ## predictions
      mah_suit <- sapply(1:length(ellipsoids), function(y) {
        mh_st <- predict(ellipsoids[[y]], back, "both")
        mh_st <- data.frame(slot(mh_st, "suitability"), slot(mh_st, "mahalanobis"))
        return(mh_st)
      })

      suits <- mah_suit[1,]
      names(suits) <- paste0("Niche_", 1:length(ellipsoids), "_S")
      suits <-  do.call(cbind, suits)

      mahas <- mah_suit[2,]
      names(mahas) <- paste0("Niche_", 1:length(ellipsoids),"_M")
      mahas <-  do.call(cbind, mahas)

      metrics <- overlap_metrics(comparison_matrix, back, mahas, suits,
                                 return_background = FALSE)

      over_metricsb <- lapply(1:length(metrics), function(y) {metrics[[y]][[1]]})
      over_metricsb <- do.call(rbind.data.frame, over_metricsb)
      rownames(over_metricsb) <- sapply(1:ncol(comparison_matrix), function(y) {
        paste0("Niche_", paste0(comparison_matrix[, y], collapse = "_vs_"))
      })
    }
    cat("\treplicate", x, "of", replicates, "\n")

    if (overlap_type[1] == "all") {
      return(list(full = over_metrics, back = over_metricsb))
    }
    if (overlap_type[1] == "full") {return(list(full = over_metrics))}
    if (overlap_type[1] == "back_union") {return(list(back = over_metricsb))}
  })

  overs <- rownames(b_over[[1]][[1]])

  if (overlap_type[1] == "all") {
    fo <- lapply(overs, function(x) {
      f <- do.call(rbind, lapply(b_over, function (y) {y[[1]][x, ]}))
      rownames(f) <- paste0("Replicate_", 1:nrow(f)); return(f)
    })
    uo <- lapply(overs, function(x) {
      u <- do.call(rbind, lapply(b_over, function (y) {y[[2]][x, ]}))
      rownames(u) <- paste0("Replicate_", 1:nrow(u)); return(u)
    })
    names(fo) <- overs
    names(uo) <- overs
    return(list(full_random = fo, union_random = uo))
  }
  if (overlap_type[1] == "full") {
    fo <- lapply(overs, function(x) {
      f <- do.call(rbind, lapply(b_over, function (y) {y[[1]][x, ]}))
      rownames(f) <- paste0("Replicate_", 1:nrow(f)); return(f)
    })
    names(fo) <- overs
    return(list(full_random = fo))
  }
  if (overlap_type[1] == "back_union") {
    uo <- lapply(overs, function(x) {
      u <- do.call(rbind, lapply(b_over, function (y) {y[[1]][x, ]}))
      rownames(u) <- paste0("Replicate_", 1:nrow(u)); return(u)
    })
    names(uo) <- overs
    return(list(union_random = uo))
  }
}
