#' Spatial thinning of occurrence data
#'
#' @description Clean duplicated longitude and latitude data using threshold distance
#'              which is a distance (in grades) between points to be considered
#'              duplicates.
#' @param data data.frame containing longitude and latitude columns.
#' @param longitude (character) vector of the column name of longitude.
#' @param latitude (character) vector of the column name of latitude.
#' @param thin_class (numerc) distance in kilometers to thin the data.
#' @param raster_layer RasterLayer
#' @param thin_distance (numeric)
#'
#' @return Returns a data.frame with coordinate data from species
#'
#' @export
#'
#' @examples
#' # GBIF search
#' occurrences <- read.csv(system.file("extdata", "occurrences.csv",
#'                         package = "kuenm"))
#' thin_occurrences <- thin_data(occurrences, longitude = "longitude",
#'                               latitude = "latitude", thin_class = 5)
#' # Check the dimensions of  data
#' dim(occurrences)
#' dim(thin_occurrences)

thin_data <- function(data, longitude, latitude, thin_class = NULL,
                      raster_layer = NULL, thin_distance = 0) {
  # -----------
  # detecting potential errors
  if (missing(data)) {
    stop("Argument data is not defined.")
  }
  if (missing(longitude)) {
    stop("Argument longitude is not defined.")
  }
  if (missing(latitude)) {
    stop("Argument latitude is not defined.")
  }

  # -----------
  # erasing NAs
  data <- data[!is.na(data[, longitude]) & !is.na(data[, latitude]), ]

  # -----------
  # erasing duplicates if raster_layer != NULL & thin_distance == 0
  if (thin_distance == 0) {
    if (!is.null(raster_layer)) {
      ids <- raster::cellFromXY(raster_layer, dat[, c(longitude, latitude)])
      data <- data[!duplicated(ids), ]
    }
  }

  # -----------
  # thinning (also earasing duplicates if raster_layer == NULL)
  WGS84 <- sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  classes <- unique(data[, thin_class]) # to get unique thinning classes
  dat_sp <- list()

  for (i in 1:length(thin_distance)) {
    dat <- data[data[, thin_class] == classes[i], ] # to get data for each class

    # to project points
    dat_s <- sp::SpatialPointsDataFrame(dat[, c(longitude, latitude)], dat,
                                        proj4string = WGS84)
    centroid <- rgeos::gCentroid(dat_s, byid = FALSE) # to get centroid of area
    AEQD <- sp::CRS(paste0("+proj=aeqd +lat_0=", centroid@coords[2], " +lon_0=",
                          centroid@coords[1], " +x_0=0 +y_0=0 +ellps=WGS84",
                          " +datum=WGS84 +units=m +no_defs"))
    dat_s <- sp::spTransform(dat_s, AEQD)

    # to remove duplicates in spatial points
    dat_s <- sp::remove.duplicates(dat_s, zero = thin_distance[i] * 1000)
    dat_sp[[i]] <- sp::spTransform(dat_s, WGS84)@data # to reproject and get data
  }

  dat_sp <- do.call(rbind, dat_sp) # to bind all thinned data

  return(dat_sp)
}
