#' Split occurrence data into training and testing data
#'
#' @description split_data splits occurrences into training and testing data
#' based on distinct methods.
#'
#' @param data data.frame of occurrence records containing at least longitude
#' and latitude columns.
#' @param method (character) method for selecting training and testing data.
#' Options are: "random" and "block"; default = "random".
#' @param longitude (character) if \code{method} = "block", name of the column
#' with longitude data.
#' @param latitude (character) if \code{method} = "block", name of the column
#' with latitude data.
#' @param train_proportion (numeric) proportion (from 0 to 1) of data to be used
#' as training occurrences. The remaining data will be used for testing.
#' Default = 0.5 if \code{method} = "random", or 0.75 if \code{method} = "block".
#' @param raster_layer optional RasterLayer to prepare background data if
#' \code{method} = "block".
#' @param background_n (numeric) optional number of coordinates to be extracted
#' using the \code{raster_layer}. Default = 10000.
#' @param save (logical) whether or not to save the results in the working
#' directory. Default = FALSE.
#' @param name (character) if \code{save} = TRUE, name of the csv files to be
#' written (comon name for all files). A suffix will be added depending on the
#' type of data: complete set, training set, or testing set of occurrences.
#' Format (.csv) is automatically added; default = "occurrences".
#'
#' @return
#' A list containing all, training, and testing occurrences. If \code{save} =
#' TRUE, three csv files will be written in the working directory according to
#' the name defined in \code{name} plus the suffix _all for all records, _train
#' for the training set, and _test for the testing set.
#'
#' If \code{method} = "block", an additional data.frame containing all data and
#' an extra column with IDs for each block will be added to the resulted list.
#' If \code{save} = TRUE, this data.frame will be written with the suffix _block.
#' If a raster layer is given in \code{raster_layer}, background coordinates
#' will be returned as part of this list. Data will be named as bg_all, bg_train,
#' bg_test, and bg_block, for all, training, testing, and all background with
#' assigned blocks, respectively.
#'
#' @usage
#' split_data(data, method = "random", longitude, latitude,
#'            train_proportion, raster_layer = NULL,
#'            background_n = 10000, save = FALSE, name = "occurrences")
#'
#' @export
#'
#' @examples
#' # reading data
#' occurrences <- read.csv(system.file("extdata", "occurrences.csv",
#'                                     package = "ellipsenm"))
#'
#' # random split 50% for trainig and 50% for testing
#' data_split <- split_data(occurrences, train_proportion = 0.5)
#'
#' names(data_split)
#' lapply(data_split, head)
#' lapply(data_split, dim)
#'
#' # random split 70% for trainig and 30% for testing
#' data_split1 <- split_data(occurrences, train_proportion = 0.7)
#'
#' names(data_split1)
#' lapply(data_split1, head)
#' lapply(data_split1, dim)
#'
#' # split 75% for trainig and 25% for testing using blocks
#' data_split2 <- split_data(occurrences, method = "block", longitude = "longitude",
#'                           latitude = "latitude", train_proportion = 0.75)
#'
#' names(data_split2)
#' lapply(data_split2, head)
#' lapply(data_split2, dim)
#'
#' # split data using blocks and preparing background
#' r_layer <- raster::raster(system.file("extdata", "bio_1.tif",
#'                                       package = "ellipsenm"))
#'
#' data_split3 <- split_data(occurrences, method = "block", longitude = "longitude",
#'                           latitude = "latitude", train_proportion = 0.75,
#'                           raster_layer = r_layer)
#'
#' # saving data
#' data_split4 <- split_data(occurrences, train_proportion = 0.7, save = TRUE,
#'                           name = "occs")
#'
#' # cheking directory
#' dir()

split_data <- function(data, method = "random", longitude, latitude,
                       train_proportion, raster_layer = NULL,
                       background_n = 10000, save = FALSE, name = "occurrences") {

  # -----------
  # detecting potential errors and defining missing proportion
  if (missing(data)) {
    stop("Argument 'data' is not defined.")
  }
  if (method == "block") {
    if (missing(longitude)) {
      stop("Argument 'longitude' is not defined.")
    }
    if (missing(latitude)) {
      stop("Argument 'latitude' is not defined.")
    }
  }

  if (missing(train_proportion)) {
    if (method == "random") train_proportion <- 0.5
    if (method == "block") train_proportion <- 0.75
  }

  # -----------
  # processing
  occ <- na.omit(data)

  if (method == "random") {
    files <- occ_randsplit(occ, train_proportion)
  }

  if (method == "block") {
    if (!train_proportion %in% c(0.25, 0.5, 0.75)) {
      stop(paste("train_proportion for block method is restricted to 0.25, 0.5,\n",
              "or 0.75. See function's help."))
    }

    if (!is.null(raster_layer)) {
      files <- occ_blocksplit(occ, longitude, latitude, train_proportion,
                              raster_layer, background_n)
    } else {
      files <- occ_blocksplit(occ, longitude, latitude, train_proportion)
    }
  }

  # -----------
  # writing results
  if (save == TRUE) {
    names <- paste0(name, "_", names(files), ".csv")
    wrt <- sapply(1:length(files), function(x){
      write.csv(files[[x]], file = names[x], row.names = FALSE)
    })

    cat("\nOccurrences were written in the working directory.\n")
  }

  return(files)
}
