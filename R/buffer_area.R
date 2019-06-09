#' Polygons based on buffered occurrences
#'
#' @description buffer_area helps in creating buffer polygons based on occurrences,
#' masking raster layers, and writing results if needed.
#'
#' @param data data.frame of occurrence records containing at least longitude
#' and latitude columns.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @param buffer_distance (numeric) distance in km to be used to create a buffer
#' for the convex hull. Default = 100.
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
#' # producing polygons
#' b_area <- buffer_area(data = occurrences, longitude = "longitude",
#'                       latitude = "latitude", buffer_distance = 100)
#'
#' # producing polygons and masking layers
#' vars <- raster::raster(system.file("extdata", "m_bio1.tif",
#'                                    package = "ellipsenm"))
#'
#' b_area1 <- buffer_area(data = occurrences, longitude = "longitude",
#'                        latitude = "latitude", buffer_distance = 100,
#'                        raster_layers = vars, mask = TRUE)
#'
#' # producing polygons, masking layers, and saving results
#' b_area2 <- buffer_area(data = occurrences, longitude = "longitude",
#'                        latitude = "latitude", buffer_distance = 100,
#'                        raster_layers = vars, mask = TRUE, save = TRUE,
#'                        name = "buff_area")

buffer_area <- function(data, longitude, latitude, buffer_distance = 100,
                        raster_layers = NULL, clip = FALSE, mask = FALSE,
                        save = FALSE, name = "calib_area_buffer") {
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
  # creating a buffer based on a user-defined distance
  buffer_distance <- buffer_distance / 111.32

  buff_area <- suppressWarnings(rgeos::gBuffer(occ_sp, width = buffer_distance))

  buff_area <- raster::disaggregate(buff_area)

  # -----------
  # getting only area of interest
  if (!is.null(raster_layers) & clip == TRUE) {
    polygons <- raster_poly(raster_layers[[1]])
    buff_area <- rgeos::gIntersection(polygons, buff_area, byid = FALSE,
                                      drop_lower_td = TRUE)
  }

  buff_area <- sp::SpatialPolygonsDataFrame(buff_area,
                                            data = data.frame(
                                              RD = rep(1, length(buff_area))),
                                            match.ID = FALSE)

  # -----------
  # masking
  if (!is.null(raster_layers) & mask == TRUE) {
    mask_var <- mask(crop(raster_layers, buff_area), buff_area)

    buff_area <- list(calibration_area = buff_area, masked_variables = mask_var)
  }

  # -----------
  # saving results
  if (save == TRUE) {
    if (!is.null(raster_layers) & mask == TRUE) {
      save_areas(name, buff_area[[1]], area_type = "buffer_area",
                 raster_layers = buff_area[[2]])
    } else {
      save_areas(name, buff_area, area_type = "buffer_area")
    }
  }

  return(buff_area)
}


#' Helper function to obtain polygon from RasterLayer
#'
#' @param raster_layer RasterLayer of a region of interest.
#'
#' @export

raster_poly <- function(raster_layer) {
  suppressPackageStartupMessages(library(raster))
  wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

  if (is.na(raster_layer@crs)) {
    raster::crs(raster_layer) <- wgs84
  } else {
    raster_layer <- raster::projectRaster(raster_layer, crs = wgs84)
  }

  raster_layer <- raster_layer > min(raster_layer[], na.rm = TRUE)
  polygons <- raster::rasterToPolygons(raster_layer, dissolve = TRUE)

  return(polygons)
}


#' Helper function to save results from calibration area creation
#'
#' @param name (character) name of a folder to be written with the results if
#' \code{save} = TRUE.
#' @param area_polygon SpatialPolygonDataFrame to be written in \code{name}.
#' @param area_type (character) type of polygona (area) to be written.
#' @param raster_layers optional RasterStack of layers to be written in \code{name}.
#'
#' @export

save_areas <- function(name, area_polygon, area_type, raster_layers = NULL) {
  dir.create(name)

  if (!is.null(raster_layers)) {
    rgdal::writeOGR(area_polygon, name, area_type,
                    driver = "ESRI Shapefile")
    rnames <- paste0(name, "/", names(raster_layers), ".tif")
    r <- lapply(1:length(rnames), function(x) {
      raster::writeRaster(raster_layers[[x]], filename = rnames[x],
                          format = "GTiff")
    })
  } else {
    rgdal::writeOGR(area_polygon, name, area_type,
                    driver = "ESRI Shapefile")
  }
}