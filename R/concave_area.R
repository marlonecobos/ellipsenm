#' Concave hull polygons based on occurrences
#'
#' @description concave_area helps in creating concave polygons based on
#' occurrences, buffering polygons, masking raster layers, and writing results
#' if needed.
#'
#' @param data data.frame of occurrence records containing at least longitude
#' and latitude columns.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @param length_threshold (numeric) distance in km to define when a segment
#' length stops being considered for further detalization. Higher values result
#' in simpler shapes. Default = 5.
#' @param split (logical) if TRUE, a distance (for hierarchical clustering) or
#' a number (for K-means clustering) is used to separate distinct chunks of
#' data Recommended when the species of interest has a disjunct distribution.
#' Default = FALSE.
#' @param n_kmeans (numeric) if \code{split} = TRUE, number of clusters in which
#' the species occurrences will be grouped using the "k-means" method.
#' Default = NULL.
#' @param buffer_distance (numeric) distance in km to be used to create a buffer
#' for the concave hull. Default = NULL.
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
#' \code{save} = TRUE. Default = "calib_area_concave".
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
#' # reading data
#' occurrences <- read.csv(system.file("extdata", "occurrences.csv",
#'                                     package = "ellipsenm"))
#'
#' # producing simple concave polygons
#' cv_area <- concave_area(data = occurrences, longitude = "longitude",
#'                        latitude = "latitude")
#'
#' sp::plot(cv_area)
#' points(occurrences[, 2:3])
#'
#' # producing concave polygons with buffers
#' cv_area1 <- concave_area(data = occurrences, longitude = "longitude",
#'                         latitude = "latitude", buffer_distance = 50)
#'
#' sp::plot(cv_area1)
#' points(occurrences[, 2:3])
#'
#' # producing concave polygons splitted considering clusters
#' cv_area2 <- concave_area(data = occurrences, longitude = "longitude",
#'                         latitude = "latitude", split = TRUE, n_kmeans = 2,
#'                         buffer_distance = 5)
#'
#' sp::plot(cv_area2)
#' points(occurrences[, 2:3])
#'
#' # producing concave polygons, masking layers
#' vars <- raster::stack(list.files(system.file("extdata", package = "ellipsenm"),
#'                                  pattern = "bio", full.names = TRUE))
#'
#' cv_area3 <- concave_area(data = occurrences, longitude = "longitude",
#'                         latitude = "latitude", buffer_distance = 50,
#'                         raster_layers = vars, mask = TRUE)
#'
#' raster::plot(cv_area3$masked_variables[[1]])
#' sp::plot(cv_area3$calibration_area, add = TRUE)
#' points(occurrences[, 2:3])
#'
#' # producing concave polygons, masking layers, and saving results
#' cv_area4 <- concave_area(data = occurrences, longitude = "longitude",
#'                         latitude = "latitude", buffer_distance = 50,
#'                         raster_layers = vars, mask = TRUE, save = TRUE,
#'                         name = "concave_area")
#'
#' # check directory
#' dir()


concave_area <- function(data, longitude, latitude, length_threshold = 5,
                         split = FALSE, n_kmeans = NULL, buffer_distance = NULL,
                         raster_layers = NULL, clip = FALSE, mask = FALSE,
                         save = FALSE, name = "calib_area_concave") {
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
  # creating concave hull(s) from points
  hulls_buffer <- lapply(1:length(occ_sps), function(x) {
    if (length(occ_sps) > 1) {
      coord <- as.data.frame(occ_sps@coords[[x]])
    } else {
      coord <- as.data.frame(sp::coordinates(occ_sps[[x]]))
    }

    if (dim(coord)[1] > 2) {
      length_threshold <- length_threshold / 111.32
      sppoints <- sf::st_multipoint(as.matrix(coord))
      sf_points <- sf::st_sf(sf::st_sfc(sppoints, crs = WGS84@projargs))
      concavehull <- concaveman::concaveman(sf_points,
                                            length_threshold = length_threshold)
      hulls <- sf::as_Spatial(sf::st_zm(concavehull$polygons))
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
    mask_var <- raster::mask(raster::crop(raster_layers, hulls_buffer), hulls_buffer)

    hulls_buffer <- list(calibration_area = hulls_buffer, masked_variables = mask_var)
  }

  # -----------
  # saving results
  if (save == TRUE) {
    if (!is.null(raster_layers) & mask == TRUE) {
      save_areas(name, hulls_buffer[[1]], area_type = "concave_area",
                 raster_layers = hulls_buffer[[2]])
    } else {
      save_areas(name, hulls_buffer, area_type = "concave_area")
    }
  }

  return(hulls_buffer)
}

