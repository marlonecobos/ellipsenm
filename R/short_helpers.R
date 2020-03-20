#' Helper function to project points
#' @param data data.frame of occurrence records containing at least longitude and
#'             latitude columns.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @export
#' @return
#' A SpatialPointsDataFrame object projected to azimuthal equidistant projection.

wgs_aeqd <- function(data, longitude, latitude) {
  WGS84 <- sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

  dat_s <- sp::SpatialPointsDataFrame(data[, c(longitude, latitude)], data,
                                      proj4string = WGS84)
  centroid <- rgeos::gCentroid(dat_s, byid = FALSE) # to get centroid of area
  AEQD <- sp::CRS(paste0("+proj=aeqd +lat_0=", centroid@coords[2], " +lon_0=",
                         centroid@coords[1], " +x_0=0 +y_0=0 +ellps=WGS84",
                         " +datum=WGS84 +units=m +no_defs"))
  dat_s <- sp::spTransform(dat_s, AEQD)

  return(dat_s)
}


#' Helper function to find raster extention
#' @param format (character) any of the format types allowed for raster objects.
#' See \code{\link[raster]{writeFormats}}
#' @export
#' @return Raster extension according to format type.

rformat_type <- function(format) {
  if (missing(format)) {stop("Argument 'format' needs to be defined.")}
  if (format == "raster") {format1 <- ".grd"}
  if (format == "GTiff") {format1 <- ".tif"}
  if (format == "EHdr") {format1 <- ".bil"}
  if (format == "ascii") {format1 <- ".asc"}
  if (format == "SAGA") {format1 <- ".sdat"}
  if (format == "IDRISI") {format1 <- ".rst"}
  if (format == "CDF") {format1 <- ".nc"}
  if (format == "ENVI") {format1 <- ".envi"}
  if (format == "HFA") {format1 <- ".img"}
  return(format1)
}


#' Split occurrences randomly in training and testing data
#' @description occ_randsplit splits a set of occurrences to obtain training and
#' testing data randomly.
#' @param data matrix or data.frame with the occurrences to be split. Columns
#' may vary but species, longitude, and latitue are recommended.
#' @param train_proportion (numeric) proportion (from 0 to 1) of data to be used
#' as training occurrences. The remaining data will be used for testing.
#' Default = 0.5.
#' @export
#' @return
#' List with all occurrences (all), training occurrences (train), and testing
#' (test) occurrences.

occ_randsplit <- function(data, train_proportion = 0.5) {
  ndata <- nrow(data)
  ids <- sample(ndata, size = round(train_proportion * ndata))
  data1 <- list(all = data, train = data[ids, ], test = data[-ids, ])

  return(data1)
}


#' Helper function to obtain polygon from RasterLayer
#' @param raster_layer RasterLayer of a region of interest.
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


#' Helper function to calculate niche volume
#' @param n_dimensions (numeric) number of dimensions to be considered.
#' @param semi_axes_length (numeric) length of ellipsoid axes.
#' @export

ellipsoid_volume <- function (n_dimensions, semi_axes_length) {
  term1 <- 2 * pi^(n_dimensions / 2)
  term2 <- n_dimensions * gamma(n_dimensions / 2)
  term3 <- prod(semi_axes_length)
  term4 <- (term1 / term2) * term3

  return(term4)
}


#' Helper function to calculate quantiles
#' @param n_data (numeric) total number of data to be considered.
#' @param level (numeric) percentage of data to be considered when creating the
#' ellipsoid that characterizes the species ecological niche. Default = 95.
#' @export

ndata_quantile <- function(n_data, level) {
  n <- floor(n_data * level)
  if (n > n_data) {n <- n_data}
  return(n)
}


#' Helper funtion to get attributes from ellipsoid lists
#' @param ellipsoids list of ellipsoid objects.
#' @param attribute (character) name of the attribute to be obtained from elements
#' in \code{ellipsoids}. Options are: method, centroid, covariance_matrix, level,
#' niche_volume, semi_axes_length, and axes_coordinates. Default = "method".
#' @export
#' @return
#' The attribute selected.

get_attribute <- function(ellipsoids, attribute = "method"){
  if (missing(ellipsoids)) {
    stop("Argument 'ellipsoids' is needed, see function's help.")
  }
  if (!attribute %in% c("method", "centroid", "covariance_matrix", "level",
                        "niche_volume", "semi_axes_length", "axes_coordinates")) {
    stop("Invalid 'attribute', see function's help.")
  }
  attr1 <- lapply(ellipsoids, function(x) {
    attr1 <- slot(x, attribute)
    if(attribute == "axes_coordinates") {attr1 <- do.call(rbind, attr1)}
    return(attr1)
  })
  return(attr1)
}


#' Helper funtion to get data for Montecarlo simulations
#' @param ellipsoids list of ellipsoid objects.
#' @param n_points (numeric) number of random points to be generated.
#' @export
#' @return
#' A total of \code{n_points} created randomly considering the limits of the
#' ellipsoids considered.

hypercube_boundaries <- function(ellipsoids, n_points = 1000000) {
  if (missing(ellipsoids)) {
    stop("Argument 'ellipsoids' is needed, see function's help.")
  }
  axis_list <- get_attribute(ellipsoids, attribute = "axes_coordinates")
  ellipsoid_axis <- do.call(rbind, axis_list)
  min_env <- apply(ellipsoid_axis, 2, min)
  max_env <- apply(ellipsoid_axis, 2, max)
  data_rand <- sapply(1:length(min_env), function (i) {
    runif(n_points, min = min_env[i], max = max_env[i])
  })
  return(data_rand)
}
