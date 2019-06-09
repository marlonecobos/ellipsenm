#' Prepare species and environmental data for overlap analysis
#'
#' @description spdata_overlap prepares species and environmental data to perform
#' niche overlap analyses.
#'
#' @param method (character) method to contruct the ellipsoid that characterizes
#' the species ecological niche. Available methods are: "covmat" and "mve".
#' See details. Default = "covmat".
#' @param data data.frame of occurrence records. Columns must be species,
#' longitude, and latitude.
#' @param raster_layers a RasterStack of environmental data to be used in niche
#' overlap ananlyses.
#'
#' @details
#' Methods details are as follows:
#'
#' "covmat"
#'
#' "mve"
#'
#' @export
#'
#' @examples
#' # data
#' met <- "covmat"
#' occurrences <- read.csv(system.file("extdata", "occurrences_thin.csv",
#'                                     package = "ellipsenm"))
#' layers <- raster::raster(system.file("extdata", "m_bio1.tif",
#'                                      package = "ellipsenm"))
#'
#' sp_over <- spdata_overlap(method = met, data = occurrences,
#'                           raster_layers = layers)

spdata_overlap <- function(method, data, raster_layers) {

  if (!missing(method)) {
    if (is.character(method) == FALSE) {
      stop("method needs to be character.")
    }
  } else {
    stop("Argument method needs to be defined.")
  }

  if (!missing(data)) {
    if (is.data.frame(data) == FALSE) {
      stop(paste("data needs to be a data.frame with at least three columns:",
                 "\nspecies, longitude, and latitude"))
    }
  } else {
    stop("Argument data needs to be defined.")
  }

  if (!missing(raster_layers)) {
    if (class(raster_layers)[1] != "RasterStack") {
      stop("raster_layers needs to be a RasterStack.")
    }
  } else {
    stop("Argument raster_layers needs to be defined.")
  }

  results <- data_overlap(method = method,
                          occurrences = data,
                          m_layers = raster_layers)

  return(results)
}
