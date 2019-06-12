#' Convex hull polygons based on occurrences
#'
#' @description convex_area helps in creating convex polygons based on
#' occurrences, buffering polygons, masking raster layers, and writing results
#' if needed.
#'
#' @param data data.frame of occurrence records containing at least longitude
#' and latitude columns.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @param split (logical) if TRUE, a distance (for hierarchical clustering) or
#' a number (for K-means clustering) is used to separate distinct chunks of
#' data Recommended when the species of interest has a disjunct distribution.
#' Default = FALSE.
#' @param n_kmeans (numeric) if \code{split} = TRUE, number of clusters in which
#' the species occurrences will be grouped using the "k-means" method.
#' Default = NULL.
#' @param buffer_distance (numeric) distance in km to be used to create a buffer
#' for the convex hull. Default = NULL.
#' @param raster_layers optional RasterStack to be used in restricting the
#' resulting SpatialPolygon and in preparing variables for further processing.
#' Default = NULL.
#' @param clip (logical) whether or not to clip polygons considering boundaries
#' of layers in \code{raster_layers}. Using this option increases the time of
#' processing considerably. Default = FALSE.
#' @param mask (logical) whether or not to mask the \code{raster_layers} to the
#' area created with thsi function. Default = FALSE.
#' @param save (logical) whether or not to write the results in the working
#' directory. Default = FALSE.
#' @param name (character) name of a folder to be written with the results if
#' \code{save} = TRUE. Default = "calib_area_convex".
#'
#' @return
#' If raster layers are masked, a lits containing a SpatialPolygonDataFrame and
#' a RasterStack of the masked layers.
#'
#' If raster layers are not masked, a SpatialPolygonDataFrame.
#'
#' If \code{save} = TRUE, results are written in a folder named as in \code{name}.
#'
#' @export
#'
#' @examples
#' # data
#' occurrences <- read.csv(system.file("extdata", "occurrences_comp.csv",
#'                                     package = "ellipsenm"))
#'
#' # producing simple convex polygons
#' cx_area <- convex_area(data = occs, longitude = "LONGITUDE",
#'                        latitude = "LATITUDE")
#'
#' # producing convex polygons with buffers
#' cx_area1 <- convex_area(data = occs, longitude = "LONGITUDE",
#'                         latitude = "LATITUDE", buffer_distance = 50)
#'
#' # producing convex polygons splitted considering clusters
#' cx_area2 <- convex_area(data = occs, longitude = "LONGITUDE",
#'                         latitude = "LATITUDE", split = TRUE, n_kmeans = 3,
#'                         buffer_distance = 50)
#'
#' # producing convex polygons, masking layers
#' vars <- raster::stack(list.files(system.file("extdata", package = "ellipsenm"),
#'                                  pattern = "m_bio", full.names = TRUE))
#'
#' cx_area3 <- convex_area(data = occs, longitude = "LONGITUDE",
#'                         latitude = "LATITUDE", buffer_distance = 100,
#'                         raster_layers = vars, mask = TRUE)
#'
#' # producing convex polygons, masking layers, and saving results
#' cx_area4 <- convex_area(data = occurrences, longitude = "longitude",
#'                        latitude = "latitude", buffer_distance = 100,
#'                        raster_layers = vars, mask = TRUE, save = TRUE,
#'                        name = "buff_area")

convex_area <- function(data, longitude, latitude, split = FALSE,
                        n_kmeans = NULL, buffer_distance = NULL,
                        raster_layers = NULL, clip = FALSE, mask = FALSE,
                        save = FALSE, name = "calib_area_convex") {
  # -----------
  # detecting potential errors
  if (missing(data)) {
    stop("Argument occurrences is necessary to perform the analysis")
  }
  if (missing(longitude)) {
    stop("Argument longitude is not defined.")
  }
  if (missing(latitude)) {
    stop("Argument latitude is not defined.")
  }

  # -----------
  # preparing data
  wgs <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  WGS84 <- sp::CRS(wgs)
  occ_sp <- sp::SpatialPointsDataFrame(coords = data[, c(longitude, latitude)],
                                       data = data, proj4string = WGS84)

  if (!is.null(raster_layers)) {
    if (is.na(raster::crs(raster_layers))) {raster::crs(raster_layers) <- wgs}
  }

  # -----------
  # Defining clusters if needed
  if (split == TRUE) {
    occ_sp <- cluster_split(occ_sp, n_kmeans = n_kmeans)
  } else {
    occ_sp@data <- data.frame(occ_sp@data, clusters = 1)
  }

  ## splitting points
  if (sum(unique(occ_sp@data$clusters)) > 1) {
    occ_sps <- split(occ_sp, occ_sp@data$clusters)
  } else {
    occ_sps <- list(occ_sp)
  }

  # -----------
  # creating convex hull(s) from points
  hulls_buffer <- lapply(1:length(occ_sps), function(x) {
    if (length(occ_sps) > 1) {
      coord <- as.data.frame(sp::coordinates(occ_sps[[x]]@coords))
    } else {
      coord <- as.data.frame(sp::coordinates(occ_sps[[x]]))
    }

    if (nrow(coord) > 2) {
      covexhull <- chull(coord) # convex hull from points
      coord_pol <- coord[c(covexhull, covexhull[1]),] # defining coordinates
      poly_list <- list(sp::Polygons(list(sp::Polygon(coord_pol)), ID = 1))
      hulls <- sp::SpatialPolygons(poly_list, proj4string = WGS84)
    } else {
      hulls <- sp::SpatialPointsDataFrame(coords = coord, data = coord,
                                          proj4string = WGS84)
    }

    return(hulls)
  })

  # -----------
  # creating buffers based on a user-defined distance if needed
  if (!is.null(buffer_distance)) {
    buffer_distance <- buffer_distance / 111.32

    hulls_buffer <- lapply(hulls_buffer, function(x) {
      hb <- suppressWarnings(rgeos::gBuffer(x, width = buffer_distance))
    })

    # union buffered polygons
    hulls_buffer <- do.call(sp::rbind.SpatialPolygons,
                            c(hulls_buffer, list(makeUniqueIDs = TRUE)))

    hulls_buffer <- raster::disaggregate(hulls_buffer)
  } else {
    hulls_buffer <- hulls_buffer[[1]]
  }

  # -----------
  # Clipping with area of interest
  if (!is.null(raster_layers) & clip == TRUE) {
    polygons <- raster_poly(raster_layers[[1]])
    hulls_buffer <- suppressWarnings(rgeos::gIntersection(hulls_buffer, polygons,
                                                          byid = FALSE,
                                                          drop_lower_td = TRUE))
  }

  hulls_buffer <- sp::SpatialPolygonsDataFrame(hulls_buffer,
                                               data = data.frame(
                                                 RD = rep(1, length(hulls_buffer))),
                                               match.ID = FALSE)

  # -----------
  # masking
  if (!is.null(raster_layers) & mask == TRUE) {
    mask_var <- mask(crop(raster_layers, hulls_buffer), hulls_buffer)

    hulls_buffer <- list(calibration_area = hulls_buffer, masked_variables = mask_var)
  }

  # -----------
  # saving results
  if (save == TRUE) {
    if (!is.null(raster_layers) & mask == TRUE) {
      save_areas(name, hulls_buffer[[1]], area_type = "convex_area",
                 raster_layers = hulls_buffer[[2]])
    } else {
      save_areas(name, hulls_buffer, area_type = "convex_area")
    }
  }

  return(hulls_buffer)
}


#' Split data based on clusters
#'
#' @param data data.frame of occurrence records containing at least longitude
#' and latitude columns.
#' @param cluster_method (character) name of the method to be used for clustering
#' the occurrences. Options are "hierarchical" and "k-means"; default =
#' "hierarchical". Note that this parameter is ignored when \code{split} = FALSE.
#' See details \code{\link[ellipsenm]{cluster_split}}.
#' @param split_distance (numeric) distance in km that will be considered as the
#' limit of connectivity among polygons created with clusters of occurrences.
#' This parameter is used when \code{cluster_method} = "hierarchical" and
#' \code{split} = TRUE. Default = NULL.
#' @param n_kmeans (numeric) if \code{split} = TRUE, number of clusters in which
#' the species occurrences will be grouped when using the "k-means"
#' \code{cluster_method}. Default = NULL.
#'
#' @details The \code{cluster_method} must be chosen based on the spatial
#' configuration of the species occurrences. Both methods make distinct assumptions
#' and one of them may perform better than the other depending on the spatial
#' pattern of the data.
#'
#' The k-means method, for example, perfomrs better when the following assumptions
#' are fulfilled: Clusters are spatially grouped—or “spherical” and Clusters are
#' of a similar size. Owing to the nature of the hierarchical clustering algorithm
#' it may take more time than the k-means method. Both methods make assumptions
#' and they may work well on some data sets, and fail on others.
#'
#' Another important factor to consider is that the k-means method allways starts
#' with a random choice of cluster centers, thus it may end in different results
#' on different runs. That may be problematic when trying to replicate your
#' methods. With hierarchical clustering, most likely the same clusters can be
#' obtained if the process is repeated.
#'
#' For more information on these clustering methods see Aggarwal and Reddy (2014)
#' \url{https://goo.gl/RQ2ebd}.
#'
#' @export

cluster_split <- function(data, cluster_method = "k-means",
                          split_distance = 250, n_kmeans = NULL) {

  if (cluster_method == "hierarchical" | cluster_method == "k-means") {
    # split groups of points based on the split distance

    if (cluster_method == "hierarchical") {
      ## defining a hierarchical cluster method for the occurrences
      cluster_method <- hclust(dist(data.frame(rownames = 1:nrow(data@data),
                                               x = sp::coordinates(data)[, 1],
                                               y = sp::coordinates(data)[, 2])),
                               method = "complete")

      ## defining wich points are clustered based on the user-defined distance
      cluster_vector <- cutree(cluster_method, h = (split_distance / 111.32))
    } else {
      set.seed(1) # to get always the same answer when using the same data

      ## identifying clusters from occurrences
      cluster_vector <- kmeans(as.matrix(sp::coordinates(data)), n_kmeans)$cluster
    }

  } else {
    stop("Options of cluster_method are: \n \"hierarchical\" or \"k-means\"")
  }

  ## Join results to occurrences
  data@data <- data.frame(data@data, clusters = cluster_vector)

  return(data)
}



