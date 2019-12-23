#' Plot results of niche overlap analyses
#'
#' @description plot_overlap helps to create two or thre dimentional representations
#' of the overlap results obtained with the function \code{\link{ellipsoid_overlap}}.
#'
#' @param object overlap_ellipsoid object resulted from using the function
#' \code{\link{ellipsoid_overlap}}.
#' @param niches (numeric) pair of integer numbers denoting the niches to be
#' plotted. Default = c(1, 2).
#' @param niche_col colors to be used to plot ellipsoids of niches to be
#' compared. Default = c("blue", "red").
#' @param data (logical) whether or not to plot points of species data. Default
#' = TRUE.
#' @param data_col colors to be used to plot data points of niches to be
#' compared. Default = c("blue", "red").
#' @param background (logical) whether or not to plot background points. Default
#' = TRUE.
#' @param background_type (character) type of background to be plotted. Options
#' are: "full", "back_union". See \code{\link{ellipsoid_overlap}}.
#' @param proportion (numeric) proportion of background to be plotted. Default = 0.3.
#' @param background_col color ramp to be used for coloring background points.
#' Default = viridis::viridis.
#' @param change_labels (logical) whether or not to change axes label.
#' Default = FALSE.
#' @param xlab (character) lable of x axis. Default = "".
#' @param ylab (character) lable of y axis. Default = "".
#' @param zlab (character) lable of z axis. Default = "".
#' @param legend (logical) whether or not to add a simple legend. Default = TRUE.
#'
#' @return A plot of the niches to be compared in environmental space.
#'
#' @usage
#' plot_overlap(object, niches = c(1, 2), niche_col = c("blue", "red"),
#'              data = TRUE, data_col = c("blue", "red"),
#'              background = FALSE, background_type, proportion = 0.3,
#'              background_col = viridis::viridis, change_labels = FALSE,
#'              xlab = "", ylab = "", zlab = "", legend = TRUE)
#'
#' @export
#'
#' @examples
#' # Preparing example
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
#' occurrences1 <- occurrences[occurrences$longitude < (mean(vext[1:2]) + 0.2),]
#' occurrences2 <- occurrences[!occurrences$longitude %in% occurrences1$longitude,]
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
#'
#' # Now the plots
#' # plotting only ellipsoids
#' plot_overlap(overlap)
#'
#' # plotting ellispodis and background for full overlap
#' plot_overlap(overlap, background = TRUE, proportion = 0.6, background_type = "full")
#'
#' # plotting ellispodis and background for overlap based on accessible environments
#' plot_overlap(overlap, background = TRUE,  proportion = 1, background_type = "back_union")

plot_overlap <- function(object, niches = c(1, 2), niche_col = c("blue", "red"),
                         data = TRUE, data_col = c("blue", "red"),
                         background = FALSE, background_type, proportion = 0.3,
                         background_col = viridis::viridis, change_labels = FALSE,
                         xlab = "", ylab = "", zlab = "", legend = TRUE) {

  # -----------
  # detecting potential errors
  if (missing(object)) {
    stop("Argument object is necessary to perform the analysis.")
  }
  if (length(niche_col) < length(object@ellipsoids)) {
    message("Number of niches to plot exceeds number of colors, using automatic selection.")
    data_col <- rainbow(length(object@ellipsoids))
    niche_col <- rainbow(length(object@ellipsoids))
  }
  if (background == TRUE & missing(background_type)) {
    stop("Argument background_type needs to be defined if background = TRUE.")
  }

  # -----------
  # preparing data
  var_names <- object@variable_names
  iter <- niches[1]:niches[2]
  backs <- paste0("Niche_", niches[1], "_vs_", niches[2])

  # -----------
  # plotting
  if(length(var_names) > 2) {
    if (data == TRUE) {
      points <- lapply(iter, function(x) {
        sp_data <- object@data[[x]][, var_names]
        if (x == niches[1]) {
          if (change_labels == TRUE) {
            rgl::plot3d(sp_data[, 1:3], col = data_col[x], size = 6, xlab = xlab,
                        ylab = ylab, zlab = zlab)
          } else {
            rgl::plot3d(sp_data[, 1:3], col = data_col[x], size = 6)
          }
        } else {
          rgl::plot3d(sp_data[, 1:3], col = data_col[x], size = 6, add = TRUE)
        }
      })
    }

    if (background == TRUE) {
      if (background_type == "full") {
        mh_sort <- object@full_background[[backs]]
      }
      if (background_type == "back_union") {
        mh_sort <- object@union_background[[backs]]
      }

      mh_sort <- mh_sort[order(mh_sort[, "Niche_1_S"]), ]
      mh_sort <- mh_sort[sample(ceiling(nrow(mh_sort) * proportion)), ]
      mh_sort <- mh_sort[order(mh_sort[, "Niche_1_S"]), 1:3]

      if (data == TRUE) {
        rgl::plot3d(mh_sort, col = background_col(nrow(mh_sort)), add = TRUE)
      } else {
        if (change_labels == TRUE) {
          rgl::plot3d(mh_sort, col = background_col(nrow(mh_sort)), xlab = xlab,
                      ylab = ylab, zlab = zlab)
        } else {
          rgl::plot3d(mh_sort, col = background_col(nrow(mh_sort)))
        }
      }
    }

    ellipsoides <- lapply(iter, function(x) {
      centroid <- object@ellipsoids[[x]]@centroid[1:3]
      cov_mat <- object@ellipsoids[[x]]@covariance_matrix[1:3, 1:3]
      level <- object@ellipsoids[[x]]@level / 100
      ell <- rgl::ellipse3d(cov_mat, centre = centroid, level = level)
      if (change_labels == TRUE) {
        rgl::wire3d(ell, col = niche_col[x], alpha = 0.5, xlab = xlab, ylab = ylab,
                    zlab = zlab)
      } else {
        rgl::wire3d(ell, col = niche_col[x], alpha = 0.5)
      }

    })

    if (legend == TRUE) {
      rgl::legend3d("topright", legend = paste("Niche", niches[1:2]),
                    lty = 1, col = niche_col, inset = 0.02, bty = "n")
    }
  } else {
    el1 <- lapply(iter, function(x) {
      centroid <- object@ellipsoids[[x]]@centroid[1:2]
      cov_mat <- object@ellipsoids[[x]]@covariance_matrix[1:2, 1:2]
      level <- object@ellipsoids[[x]]@level / 100
      ellipse::ellipse(x = cov_mat, centre = centroid, level = level)
    })

    xlim <- range(unlist(lapply(el1, function(x) {x[, 1]})))
    ylim <- range(unlist(lapply(el1, function(x) {x[, 2]})))

    par(mar = c(4, 4, 1, 1))
    if (background == TRUE) {
      if (background_type == "full") {
        mh_sort <- object@full_background[[backs]]
      }
      if (background_type == "back_union") {
        mh_sort <- object@union_background[[backs]]
      }

      mh_sort <- mh_sort[order(mh_sort[, "Niche_1_S"]), ]
      mh_sort <- mh_sort[sample(ceiling(nrow(mh_sort) * proportion)), ]
      mh_sort <- mh_sort[order(mh_sort[, "Niche_1_S"]), 1:2]

      if (change_labels == TRUE) {
        plot(mh_sort, col = background_col(nrow(mh_sort)), xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab)
      } else {
        plot(mh_sort, col = background_col(nrow(mh_sort)), xlim = xlim, ylim = ylim,
             xlab = var_names[1], ylab = var_names[2])
      }
    }

    if (data == TRUE) {
      points <- lapply(iter, function(x) {
        sp_data <- object@data[[x]][, var_names]
        if (x == niches[1]) {
          if (background == TRUE) {
            points(sp_data[, 1:2], pch = 19, col = data_col[x])
          } else {
            if (change_labels == TRUE) {
              plot(sp_data[, 1:2], pch = 19, col = data_col[x], xlim = xlim, ylim = ylim,
                   xlab = xlab, ylab = ylab)
            } else {
              plot(sp_data[, 1:2], pch = 19, col = data_col[x], xlim = xlim, ylim = ylim,
                   xlab = var_names[1], ylab = var_names[2])
            }
          }
        } else {
          points(sp_data[, 1:2], pch = 19, col = data_col[x])
        }
      })
    }

    ellipsoides <- lapply(iter, function(x) {
      lines(el1[[x]], col = niche_col[x], lwd = 1.5)
    })

    if (legend == TRUE) {
      legend("topright", legend = paste("Niche", niches[1:2]), lty = 1,
             col = niche_col, bty = "n", horiz = TRUE)
    }
  }

}
