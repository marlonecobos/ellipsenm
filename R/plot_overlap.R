#' Plot results of niche overlap analyses
#'
#' @description plot_overlap
#'
#' @param object overlap_ellipsoid object resulted from using the function
#' \code{\link{ellipsoid_overlap}}.
#' @param niches (numeric) pair of integer numbers denoting the niches to be
#' plotted. Default = c(1, 2).
#' @param data (logical) whether or not to plot points of species data. Default
#' = TRUE.
#' @param col colors to be used to plot data and ellipsoids of niches to be
#' compared. Default = c("blue", "red").
#' @param background (logical) whether or not to plot background points. Default
#' = TRUE.
#' @param background_type (character) type of background to be plotted. Options
#' are: "full", "back_union". See \code{\link{ellipsoid_overlap}}.
#' @param proportion (numeric) proportion of background to be plotted. Default = 0.3.
#' @param background_col color ramp to be used for coloring background points.
#' Default = viridis::viridis.
#' @param xlab (character) lable of x axis. Default = "".
#' @param ylab (character) lable of y axis. Default = "".
#' @param zlab (character) lable of z axis. Default = "".
#' @param legend (logical) whether or not to add a simple legend. Default = TRUE.
#'
#' @return A plot of the niches to be compared in environmental space.
#'
#' @export
#'
#' @examples
#' # data
#'
#' #simple plot

plot_overlap <- function(object, niches = c(1, 2), data = TRUE,
                         col = c("blue", "red"), background = FALSE,
                         background_type, proportion = 0.3,
                         background_col = viridis::viridis, xlab = "", ylab = "",
                         zlab = "", legend = TRUE) {

  # -----------
  # detecting potential errors
  if (missing(object)) {
    stop("Argument object is necessary to perform the analysis.")
  }
  if (length(col) < length(object@ellipsoids)) {
    message("Number of niches to plot exceeds number of colors, using automatic selection.")
    col <- rainbow(length(object))
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
          rgl::plot3d(sp_data[, 1:3], col = col[x], size = 6, xlab = xlab,
                      ylab = ylab, zlab = zlab)
        } else {
          rgl::plot3d(sp_data[, 1:3], col = col[x], size = 6, add = TRUE)
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
        rgl::plot3d(mh_sort, col = background_col(nrow(mh_sort)), xlab = xlab,
                    ylab = ylab, zlab = zlab)
      }
    }

    ellipsoides <- lapply(iter, function(x) {
      centroid <- object@ellipsoids[[x]]@centroid[1:3]
      cov_mat <- object@ellipsoids[[x]]@covariance_matrix[1:3, 1:3]
      level <- object@ellipsoids[[x]]@level / 100
      ell <- rgl::ellipse3d(cov_mat, centre = centroid, level = level)
      rgl::wire3d(ell, col = col[x], alpha = 0.5, xlab = xlab, ylab = ylab,
                  zlab = zlab)
    })

    if (legend == TRUE) {
      rgl::legend3d("topright", legend = paste("Niche", niches[1:2]),
                    lty = 1, col = col, inset = 0.02, bty = "n")
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

      plot(mh_sort, col = background_col(nrow(mh_sort)), xlim = xlim, ylim = ylim,
           xlab = var_names[1], ylab = var_names[2])
    }

    if (data == TRUE) {
      points <- lapply(iter, function(x) {
        sp_data <- object@data[[x]][, var_names]
        if (x == niches[1]) {
          if (background == TRUE) {
            points(sp_data[, 1:2], pch = 19, col = col[x])
          } else {
            plot(sp_data[, 1:2], pch = 19, col = col[x], xlim = xlim, ylim = ylim,
                 xlab = var_names[1], ylab = var_names[2])
          }
        } else {
          points(sp_data[, 1:2], pch = 19, col = col[x])
        }
      })
    }

    ellipsoides <- lapply(iter, function(x) {
      lines(el1[[x]], col = col[x], lwd = 1.5)
    })

    if (legend == TRUE) {
      legend("topright", legend = paste("Niche", niches[1:2]), lty = 1,
             col = col, bty = "n", horiz = TRUE)
    }
  }

}
