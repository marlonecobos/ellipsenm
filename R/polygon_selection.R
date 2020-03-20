#' Selection of spatial polygons based on occurrences
#'
#' @description polygon_selection helps in selecting polygons that intersect
#' with occurrences. Buffering polygons, masking raster layers, and writing
#' results can also be performed if needed.
#'
#' @param data data.frame of occurrence records containing at least longitude
#' and latitude columns.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @param polygons SpatialPolygon* object.
#' @param buffer_distance (numeric) distance in km to be used to create a buffer
#' for the concave hull. Default = NULL.
#' @param raster_layers optional RasterStack to be used in restricting the
#' resulting SpatialPolygon and in preparing variables for further processing.
#' Default = NULL.
#' @param clip (logical) whether or not to clip polygons considering boundaries
#' of layers in \code{raster_layers}. Using this option increases the time of
#' processing considerably. Default = FALSE.
#' @param dissolve (logical) whether or not to disolve polygons. Default = FALSE.
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
#' @usage
#' polygon_selection(data, longitude, latitude, polygons,
#'                   buffer_distance = NULL, raster_layers = NULL,
#'                   clip = FALSE, dissolve = FALSE, mask = FALSE,
#'                   save = FALSE, name = "calib_area_poligons")
#'
#' @export
#'
#' @examples
#' # reading data
#' occurrences <- read.csv(system.file("extdata", "occurrences.csv",
#'                                     package = "ellipsenm"))
#'
#' spatial_poly <- #HERE YOUR SPATIAL POLYGON OBJECT
#'
#' # selecting polygons
#' sp_area <- polygon_selection(data = occurrences, longitude = "longitude",
#'                              latitude = "latitude", polygons = spatial_poly)
#'
#' sp::plot(sp_area)
#' points(occurrences[, 2:3])
#'
#' # selecting polygons with buffers
#' sp_area1 <- polygon_selection(data = occurrences, longitude = "longitude",
#'                               latitude = "latitude", polygons = spatial_poly,
#'                               buffer_distance = 50)
#'
#' sp::plot(sp_area1)
#' points(occurrences[, 2:3])
#'
#' # selecting polygons, masking layers
#' vars <- raster::stack(list.files(system.file("extdata", package = "ellipsenm"),
#'                                  pattern = "bio", full.names = TRUE))
#'
#' sp_area3 <- polygon_selection(data = occurrences, longitude = "longitude",
#'                               latitude = "latitude", polygons = spatial_poly,
#'                               buffer_distance = 50, raster_layers = vars,
#'                               mask = TRUE)
#'
#' raster::plot(sp_area3$masked_variables[[1]])
#' sp::plot(sp_area3$calibration_area, add = TRUE)
#' points(occurrences[, 2:3])
#'
#' # selecting polygons, masking layers, and saving results
#' sp_area4 <- polygon_selection(data = occurrences, longitude = "longitude",
#'                               latitude = "latitude", polygons = spatial_poly,
#'                               buffer_distance = 50, raster_layers = vars,
#'                               mask = TRUE, save = TRUE, name = "p_selection")
#'
#' # check directory
#' dir()

polygon_selection <- function(data, longitude, latitude, polygons,
                              buffer_distance = NULL, raster_layers = NULL,
                              clip = FALSE, dissolve = FALSE, mask = FALSE,
                              save = FALSE, name = "calib_area_poligons") {
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
  if (missing(polygons)) {
    stop("Argument 'polygons' is not defined.")
  }

  # -----------
  # preparing data
  wgs <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  WGS84 <- sp::CRS(wgs)
  occ_sp <- sp::SpatialPointsDataFrame(coords = data[, c(longitude, latitude)],
                                       data = data, proj4string = WGS84)

  if (is.na(sp::proj4string(polygons))) {
    sp::proj4string(polygons) <- WGS84
  } else {
    if (sp::proj4string(polygons) != wgs) {
      polygons <- sp::spTransform(polygons, WGS84)
    }
  }

  if (!is.null(raster_layers)) {
    if (is.na(raster::crs(raster_layers))) {raster::crs(raster_layers) <- wgs}
  }

  # -----------
  # selection of polygons
  hulls_buffer <- polygons[occ_sp, ]

  # -----------
  # creating buffers based on a user-defined distance if needed
  if (!is.null(buffer_distance)) {
    buffer_distance <- buffer_distance / 111.32

    hulls_buffer <- suppressWarnings(rgeos::gBuffer(hulls_buffer,
                                                    width = buffer_distance))
  } else {
    hulls_buffer <- hulls_buffer
  }

  # -----------
  # Clipping with area of interest
  if (!is.null(raster_layers) & clip == TRUE) {
    polygons <- raster_poly(raster_layers[[1]])
    hulls_buffer <- suppressWarnings(rgeos::gIntersection(hulls_buffer, polygons,
                                                          byid = FALSE,
                                                          drop_lower_td = TRUE))
  }

  if (class(hulls_buffer)[1] != "SpatialPolygonsDataFrame") {
    hulls_buffer <- sp::SpatialPolygonsDataFrame(hulls_buffer,
                                                 data = data.frame(
                                                   RD = rep(1, length(hulls_buffer))),
                                                 match.ID = FALSE)
  }

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
      save_areas(name, hulls_buffer[[1]], area_type = "polygon_selection",
                 raster_layers = hulls_buffer[[2]])
    } else {
      save_areas(name, hulls_buffer, area_type = "polygon_selection")
    }
  }

  return(hulls_buffer)
}

