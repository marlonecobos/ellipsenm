buffer_area <- function(data, buffer_distance = 100, raster_layer = NULL,
                        save_shp = FALSE, name = "calibration_area_buffer") {
  # -----------
  # detecting potential errors
  if (missing(data)) {
    stop("Argument occurrences is necessary to perform the analysis")
  }

  # -----------
  # making a spatial object from coordinates
  WGS84 <- sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  occ_sp <- sp::SpatialPointsDataFrame(coords = occ[, c(longitude, latitude)],
                                       data = occ, proj4string = WGS84)

  # -----------
  # getting a polygon from raster_layer
  if (!is.null(raster_layer)) {
    suppressPackageStartupMessages(library(raster))

    if (is.na(raster_layer@crs)) {
      raster::crs(raster_layer) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    } else {
      raster_layer <- raster::projectRaster(raster_layer,
                                            crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    }

    raster_layer[] <- !is.na(raster_layer[])
    polygon <- raster::rasterToPolygons(raster_layer)
  }

  # -----------
  # creating a buffer based on a user-defined distance
  occ_sp <- occ_sp[polygon, ] # for keeping only records in area of interest

  buffer_distance <- nose ######################

  buff_area <- rgeos::gBuffer(occ_sp, width = buffer_distance)

  buff_area <- raster::disaggregate(buff_area)

  # -----------
  # getting only area of interest
  polygon <- rgeos::gUnaryUnion(polygon) # check if needed ####################

  clip_area <- rgeos::gIntersection(polygon, buff_area, byid = FALSE, # check byid ################
                                    drop_lower_td = TRUE)
}
