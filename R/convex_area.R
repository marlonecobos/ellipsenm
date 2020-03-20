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
#' @param cluster_method (character) name of the method to be used for clustering
#' the occurrences. Options are "hierarchical" and "k-means"; default =
#' "hierarchical". Note that this parameter is ignored when \code{split} = FALSE.
#' See details \code{\link{cluster_split}}.
#' @param split_distance (numeric) distance in km that will be considered as the
#' limit of connectivity among polygons created with clusters of occurrences.
#' This parameter is used when \code{cluster_method} = "hierarchical" and
#' \code{split} = TRUE. Default = NULL.
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
#' @usage
#' convex_area(data, longitude, latitude, split = FALSE,
#'             cluster_method = "hierarchical", split_distance = NULL,
#'             n_kmeans = NULL, buffer_distance = NULL,
#'             raster_layers = NULL, clip = FALSE, mask = FALSE,
#'             save = FALSE, name = "calib_area_convex")
#'
#' @export
#'
#' @examples
#' # reading data
#' occurrences <- read.csv(system.file("extdata", "occurrences.csv",
#'                                     package = "ellipsenm"))
#'
#' # producing simple convex polygons
#' cx_area <- convex_area(data = occurrences, longitude = "longitude",
#'                        latitude = "latitude")
#'
#' sp::plot(cx_area)
#' points(occurrences[, 2:3])
#'
#' # producing convex polygons with buffers
#' cx_area1 <- convex_area(data = occurrences, longitude = "longitude",
#'                         latitude = "latitude", buffer_distance = 50)
#'
#' sp::plot(cx_area1)
#' points(occurrences[, 2:3])
#'
#' # producing convex polygons splitted considering clusters
#' cx_area2 <- convex_area(data = occurrences, longitude = "longitude",
#'                         latitude = "latitude", split = TRUE,
#'                         cluster_method = "k-means", n_kmeans = 2,
#'                         buffer_distance = 5)
#'
#' sp::plot(cx_area2)
#' points(occurrences[, 2:3])
#'
#' # producing convex polygons, masking layers
#' vars <- raster::stack(list.files(system.file("extdata", package = "ellipsenm"),
#'                                  pattern = "bio", full.names = TRUE))
#'
#' cx_area3 <- convex_area(data = occurrences, longitude = "longitude",
#'                         latitude = "latitude", buffer_distance = 50,
#'                         raster_layers = vars, mask = TRUE)
#'
#' raster::plot(cx_area3$masked_variables[[1]])
#' sp::plot(cx_area3$calibration_area, add = TRUE)
#' points(occurrences[, 2:3])
#'
#' # producing convex polygons, masking layers, and saving results
#' cx_area4 <- convex_area(data = occurrences, longitude = "longitude",
#'                         latitude = "latitude", buffer_distance = 50,
#'                         raster_layers = vars, mask = TRUE, save = TRUE,
#'                         name = "convex_area")
#'
#' # check directory
#' dir()

convex_area <- function(data, longitude, latitude, split = FALSE,
                        cluster_method = "hierarchical", split_distance = NULL,
                        n_kmeans = NULL, buffer_distance = NULL,
                        raster_layers = NULL, clip = FALSE, mask = FALSE,
                        save = FALSE, name = "calib_area_convex") {
  # -----------
  # detecting potential errors
  if (missing(data)) {
    stop("Argument 'data' is necessary to perform the analysis")
  }
  if (missing(longitude)) {
    stop("Argument 'longitude' is not defined.")
  }
  if (missing(latitude)) {
    stop("Argument 'latitude' is not defined.")
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
    if (!is.null(n_kmeans) | !is.null(split_distance)) {
      if (cluster_method == "hierarchical" & is.null(split_distance)) {
        stop("If 'cluster_method' = 'hierarchical', 'split_distance' must be defined")
      }
      if (cluster_method == "k-means" & is.null(n_kmeans)) {
        stop("If 'cluster_method' = 'k-means', 'n_kmeans' must be defined")
      }
      occ_sp <- cluster_split(occ_sp, cluster_method, split_distance, n_kmeans)
    } else {
      stop("If 'split' = TRUE, 'split_distance' or 'n_kmeans' must be defined.")
    }
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
      coord <- as.data.frame(occ_sps@coords[[x]])
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

    #hulls_buffer <- raster::disaggregate(hulls_buffer)
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
    mask_var <- raster::mask(raster::crop(raster_layers, hulls_buffer), hulls_buffer)

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
