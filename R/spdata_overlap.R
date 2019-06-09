#' Prepare species and environmental data for overlap analysis
#'
#' @description

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
