plot_ellipsoids3d <- function(object, level = 0.95, ocurrences = NULL,
                              col = c("black", "red", "darkgreen"),
                              background_points = NULL, proportion = 0.3,
                              background_col = NULL) {

  # -----------
  # detecting potential errors
  if (missing(object)) {
    stop("Argument object is necessary to perform the analysis.")
  }

  # -----------
  # plotting
  if (!is.null(sp_coords)) {
    points <- lapply(1:length(sp_coords), function(x) {
      if (x == 1) {
        plot3d(sp_coords[[x]], col = col[x], size = 6)
      } else {
        plot3d(sp_coords[[x]], col = col[x], size = 6, add = TRUE)
      }
    })
  }

  ellipsoides <- lapply(1:length(object), function(x) {
    meta_data <- object[[x]]
    ell <- ellipse3d(meta_data$covariance[, 1:3],
                     centre = meta_data$centroid[1:3], level = level)

    wire3d(ell, col = col[x], alpha = 0.1)
  })

  if(!is.null(background_points)) {
    mh_sort <- background_points[order(background_points[, 4],decreasing = F), ]
    mh_sort <- mh_sort[sample(ceiling(nrow(mh_sort) * proportion)), ]
    mh_sort <- mh_sort[order(mh_sort[, 4], decreasing = F), ]

    if(is.null(background_col)){
      rainB1 <- rev(c("#002dff","#00ff5d", "#fff000","#ffd200", "#ff0000"))
      colfunc <- colorRampPalette(rainB1)
      print(rainB1)
      col1 <- rev(colfunc(nrow(mh_sort)))
    } else {
      col1 <- background_col
    }

    if (!is.null(sp_coords)) {
      plot3d(mh_sort, col = col1, add = TRUE)
    } else {
      plot3d(mh_sort, col = col1)
    }
  }

  ellipsoides <- lapply(1:length(object), function(x) {
    meta_data <- object[[x]]
    ell <- ellipse3d(meta_data$covariance[, 1:3],
                     centre = meta_data$centroid[1:3],level = level)

    wire3d(ell, col = col[x], alpha = 0.5)
  })
}
