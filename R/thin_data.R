#' Spatial thinning of occurrence data
#'
#' @description thin_data rarefies spatially occurrence data using a distance
#' in kilometers. Multiple distances are alowed if distinct classes
#' are defined.
#'
#' @param data data.frame of occurrence records containing at least longitude and
#' latitude columns.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @param thin_class (character) name of optional column with numeric values to
#' represent classes for thinning data using distinct distances.
#' @param raster_layer optional RasterLayer to define duplicates based on its
#' cell size.
#' @param thin_distance (numeric) distance in kilometers to thin the data.
#' Default = 0. If distinct classes are defined in \code{thin_class}, a vector
#' of length equal to the number of classes is required.
#' @param save (logical) whether or not to save the results in the working
#' directory. Default = FALSE.
#' @param name (character) if \code{save} = TRUE, name of the csv file to be
#' written; format (.csv) is atomatically added. Default = "occurrences_thin".
#'
#' @return
#' A data.frame with thinned data and messages about how many records were erased
#' and kept.
#'
#' @usage
#' thin_data(data, longitude, latitude, thin_class = NULL,
#'           raster_layer = NULL, thin_distance = 0, save = FALSE,
#'           name = "occurrences_thin")
#'
#' @export
#'
#' @examples
#' # reading data
#' occurrences <- read.csv(system.file("extdata", "occurrences_gbif.csv",
#'                                     package = "ellipsenm"))
#' colnames(occurrences)
#'
#' # simple thinning based on one distance
#' thin_occurrences <- thin_data(occurrences, longitude = "longitude",
#'                               latitude = "latitude", thin_distance = 5)
#'
#' # thinning using a raster (only duplicates are erased)
#' r_layer <- raster::raster(system.file("extdata", "bio_1.tif",
#'                                        package = "ellipsenm"))
#'
#' thin_occurrences1 <- thin_data(occurrences, longitude = "longitude",
#'                                latitude = "latitude", raster_layer = r_layer)
#'
#' # thinning with different distances according to distinct classes
#' occurrences1 <- read.csv(system.file("extdata", "occurrences_spa_cor.csv",
#'                                     package = "ellipsenm"))
#' colnames(occurrences1)
#' occurrences1$thin_factor # to check classes (1 = 5 km, 2 = 10 km, 3 = 15 km)
#'
#' thin_occurrences2 <- thin_data(occurrences1, longitude = "longitude",
#'                                latitude = "latitude", thin_class = "thin_factor",
#'                                thin_distance = c(5, 10, 15))
#'
#' # saving thinned data
#' thin_occurrences2 <- thin_data(occurrences1, longitude = "longitude",
#'                                latitude = "latitude", thin_class = "thin_factor",
#'                                thin_distance = c(5, 10, 15), save = TRUE,
#'                                name = "thin_data")
#'
#' # Check dimensions of original data and results
#' dim(occurrences)
#' dim(thin_occurrences)
#' dim(thin_occurrences1)
#' dim(occurrences1)
#' dim(thin_occurrences2)

thin_data <- function(data, longitude, latitude, thin_class = NULL,
                      raster_layer = NULL, thin_distance = 0, save = FALSE,
                      name = "occurrences_thin") {
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
  n <- nrow(data) # to obtain total number of records

  # -----------
  # erasing duplicates if raster_layer != NULL & thin_distance == 0
  if (length(thin_distance) == 1) {
    if (thin_distance == 0) {
      if (!is.null(raster_layer)) {
        ids <- raster::cellFromXY(raster_layer, data[, c(longitude, latitude)])
        data <- data[!duplicated(ids), ]
      }
    }
  }

  # -----------
  # thinning (also earasing duplicates if raster_layer == NULL)
  WGS84 <- sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  cat("\nOriginal number of non-NA records:\t", n, "\n")

  if (is.null(thin_class)) {
    # to project points
    dat_s <- wgs_aeqd(data, longitude, latitude)

    # to remove duplicates in spatial points
    dat_s <- sp::remove.duplicates(dat_s, zero = thin_distance * 1000)
    dat_sp <- sp::spTransform(dat_s, WGS84)@data # to reproject and get data

    cat("\nNumber of erased records:\t", n - nrow(dat_sp), "\n")

  } else {
    classes <- sort(unique(data[, thin_class])) # to get unique thinning classes
    dat_sp <- list()

    if (length(classes) != length(thin_distance)) {
      warning(paste("Length of classes does not match length of thin_distance, the",
                    "first\nvalue of thin_distance will be used for all classes."))
      thin_distance <- thin_distance[1]
    }

    for (i in 1:length(thin_distance)) {
      dat <- data[data[, thin_class] == classes[i], ] # to get data for each class

      # to project points
      dat_s <- wgs_aeqd(dat, longitude, latitude)

      # to remove duplicates in spatial points
      dat_s <- sp::remove.duplicates(dat_s, zero = thin_distance[i] * 1000)
      dat_sp[[i]] <- sp::spTransform(dat_s, WGS84)@data # to reproject and get data

      # to report cleaning results
      cat("\nOriginal number of records for class", paste0(i, ":\t"),
          nrow(dat))
      cat("\nNumber of erased records for class", paste0(i, ":\t"),
          nrow(dat) - nrow(dat_sp[[i]]), "\n")
    }

    dat_sp <- do.call(rbind, dat_sp) # to bind all thinned data
  }

  cat("\nTotal number of thinned records:\t", nrow(dat_sp), "\n")

  # -----------
  # writing results
  if (save == TRUE) {
    cat("\nOccurrences were written in the working directory.\n")
    write.csv(dat_sp, file = paste0(name, ".csv"), row.names = FALSE)
  }

  return(dat_sp)
}

