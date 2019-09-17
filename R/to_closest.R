#' Move occurrences to closest pixel with environmental data
#'
#' @description to_closest helps in changing the longitude and latitude values
#' of occurrences with no environmental data, so they move to the closest pixel
#' of a raster layer that contains relevant information. This process prevents
#' NAs in future analyses.
#'
#' @param data data.frame or matrix of occurrence records. Columns must include
#' longitude and latitude. Other columns are optional and wont be changed.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @param raster_layer RasterLayer to be used as a reference.
#' @param limit_distance (numeric) maximun distance in km at which an occurrence
#' could be to be moved. Records farther than this distance wont be moved.
#'
#' @return
#' A data.frame with the corrected coordinates and four additional columns.
#' The first of the new columns indicates the condition of the coordinates:
#' Correct, if it was not moved because it was on a pixel with data; Moved, if
#' it was moved to the nearest pixel; and Not_moved, if it was not moved because
#' the occurrence was farther than the \code{limit_distance} to the closest pixel.
#' The second new column indicates the distance to the closest pixel with data.
#' The other two additional columns will contain the initial longitudes and
#' latitudes.
#'
#' @export
#'
#' @examples
#' data <- read.csv(system.file("extdata", "occurrences.csv",
#'                              package = "ellipsenm"))
#'
#' var <- raster::raster(list.files(system.file("extdata", package = "ellipsenm"),
#'                                  pattern = "bio", full.names = TRUE)[1])
#'
#' raster::plot(var)
#'
#' out <- data.frame(as.character(data$species[1]), rbind(c(-103, 23), c(-99, 22.5),
#'                                                        c(-103.5, 25.6), c(-99, 27.3)))
#' colnames(out) <- colnames(data)
#'
#' data <- rbind.data.frame(data, out)
#'
#' points(data[, 2:3])
#'
#' data1 <- to_closest(data, longitude = "longitude", latitude = "latitude",
#'                     raster_layer = var, limit_distance = 70)
#'
#' points(data1[, 2:3], col = "red")

to_closest <- function(data, longitude, latitude, raster_layer, limit_distance) {
  # -----------
  # detecting potential errors
  if (missing(data)) {
    stop("Argument data is necessary to perform the analysis")
  }
  if (missing(longitude)) {
    stop("Argument longitude is not defined.")
  }
  if (missing(latitude)) {
    stop("Argument latitude is not defined.")
  }
  if (missing(raster_layer)) {
    stop("Argument raster_layer is not defined.")
  }
  if (missing(limit_distance)) {
    stop("Argument limit_distance is not defined.")
  }

  # -----------
  # preparing data
  xy <- data[, c(longitude, latitude)]
  vals <- raster::extract(raster_layer, xy)

  tomove <- which(is.na(vals))
  xyout <- xy[tomove, ]

  xyras <- raster::rasterToPoints(raster_layer)[, 1:2]
  dists <- raster::pointDistance(xyout, xyras, lonlat = TRUE)

  condition <- rep("Correct", nrow(data))
  distss <- rep(0, nrow(data))

  limdist <- limit_distance * 1000

  # -----------
  # running process
  cat("\nMoving occurrences to closest pixels:\n")
  for (i in 1:nrow(xyout)) {
    mindis <- min(dists[i, ])[1]
    if (mindis <= limdist) {
      xyin <- xyras[dists[i, ] == mindis, ]
      if (class(xyin)[1] == "matrix") {
        xyin <- xyin[1, ]
      }
      data[tomove[i], longitude] <- xyin[1]
      data[tomove[i], latitude] <- xyin[2]
      condition[tomove[i]] <- "Moved"
      distss[tomove[i]] <- mindis / 1000
      cat("\tOccurrence", i, "of", nrow(xyout), "moved\n")
    } else {
      condition[tomove[i]] <- "Not_moved"
      distss[tomove[i]] <- mindis / 1000
      cat(paste0("\tOccurrence ", i," of ", nrow(xyout), " was not moved because it is more than\n\t",
                 limit_distance, " km apart from the closest pixel with environmental values\n"))
    }
  }
  data <- data.frame(data, condition = condition, distance_km = distss,
                     initial_lon = xy[, 1], initial_lat = xy[, 2],
                     stringsAsFactors = FALSE)

  return(data)
}
