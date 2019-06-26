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
    stop("Argument data is not defined.")
  }
  if (method == "block") {
    if (missing(longitude)) {
      stop("Argument longitude is not defined.")
    }
    if (missing(latitude)) {
      stop("Argument latitude is not defined.")
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



#' Split occurrences randomly in training and testing data
#'
#' @description occ_randsplit splits a set of occurrences to obtain training and
#' testing data randomly.
#'
#' @param data matrix or data.frame with the occurrences to be split. Columns
#' may vary but species, longitude, and latitue are recommended.
#' @param train_proportion (numeric) proportion (from 0 to 1) of data to be used
#' as training occurrences. The remaining data will be used for testing.
#' Default = 0.5.
#'
#' @export
#'
#' @return
#' List with all occurrences (all), training occurrences (train), and testing
#' (test) occurrences.

occ_randsplit <- function(data, train_proportion = 0.5) {
  ndata <- nrow(data)
  ids <- sample(ndata, size = round(train_proportion * ndata))
  data1 <- list(all = data, train = data[ids, ], test = data[-ids, ])

  return(data1)
}



#' Split occurrences in training and testing data using blocks
#'
#' @description occ_blocksplit splits a set of occurrences to obtain training
#' and testing data using blocks.
#'
#' @param data matrix or data.frame with the occurrences to be split. Columns
#' may vary but species, longitude, and latitue are recommended.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @param raster_layer optional RasterLayer to prepare background data.
#' @param background_n (numeric) optional number of coordinates to be extracted
#' using the \code{raster_layer}. Default = 10000.
#' @param train_proportion (numeric) proportion of data to be used as training
#' occurrences. Available options are 0.25, 0.5, and 0.75. The remaining data
#' will be used for testing. Default = 0.75.
#'
#' @export
#'
#' @return
#' List with all occurrences (all), training occurrences (train), testing
#' (test) occurrences, and all occurrences with assigned blocks (block).
#'
#' If a raster layer is given in \code{raster_layer}, background coordinates
#' will be returned as part of this list. Data will be named as bg_all, bg_train,
#' bg_test, and bg_block, for all, training, testing, and all background with
#' assigned blocks, respectively.

occ_blocksplit <- function(data, longitude, latitude, train_proportion = 0.75,
                           raster_layer = NULL, background_n = 10000) {
  ndata <- nrow(data)

  # -----------
  # occurrences split
  n1 <- ceiling(nrow(data) / 2)
  n2 <- floor(nrow(data) / 2)
  n3 <- ceiling(n1 / 2)
  n4 <- ceiling(n2 / 2)
  grp_a <- data[order(data[, latitude]), ][1:n1, ]
  grp_b <- data[rev(order(data[, latitude])), ][1:n2, ]
  grp1 <- grp_a[order(grp_a[, longitude]), ][1:(n3), ]
  grp2 <- grp_a[!rownames(grp_a) %in% rownames(grp1), ]
  grp3 <- grp_b[order(grp_b[, longitude]), ][1:(n4), ]
  grp4 <- grp_b[!rownames(grp_b) %in% rownames(grp3), ]

  # -----------
  # background split
  if (!is.null(raster_layer)) {
    back <- as.data.frame(sp::coordinates(raster_layer)[!is.na(raster_layer[]), ])

    if (nrow(back) > background_n) {
      back <- back[sample(nrow(back), background_n), ]
    }
    colnames(back) <- c(longitude, latitude)

    bvert <- mean(c(max(grp1[, longitude]), min(grp2[, longitude])))
    tvert <- mean(c(max(grp3[, longitude]), min(grp4[, longitude])))
    horz <- mean(c(max(grp_a[, latitude]), min(grp_b[, latitude])))
    bggrp1 <- back[back[, latitude] <= horz & back[, longitude] < bvert, ]
    bggrp2 <- back[back[, latitude] < horz & back[, longitude] >= bvert, ]
    bggrp3 <- back[back[, latitude] > horz & back[, longitude] <= tvert, ]
    bggrp4 <- back[back[, latitude] >= horz & back[, longitude] > tvert, ]
  }

  # -----------
  # preparing data
  ## occurrences
  if (nrow(grp1) > 0) grp1$grp <- 1
  if (nrow(grp2) > 0) grp2$grp <- 2
  if (nrow(grp3) > 0) grp3$grp <- 3
  if (nrow(grp4) > 0) grp4$grp <- 4

  data_b <- rbind(grp1, grp2, grp3, grp4)
  colnames(data_b) <- c(colnames(data), "Block")

  ## object to return data
  if (!is.null(raster_layer)) {
    ## background
    if (nrow(bggrp1) > 0) bggrp1$grp <- 1
    if (nrow(bggrp2) > 0) bggrp2$grp <- 2
    if (nrow(bggrp3) > 0) bggrp3$grp <- 3
    if (nrow(bggrp4) > 0) bggrp4$grp <- 4

    back_b <- rbind(bggrp1, bggrp2, bggrp3, bggrp4)
    colnames(back_b) <- c(colnames(back), "Block")

    data1 <- list(all = data,
                  train = data_b[data_b$Block != 4, colnames(data_b) != "Block"],
                  test = data_b[data_b$Block == 4, colnames(data_b) != "Block"],
                  block = data_b,
                  bg_all = back,
                  bg_train = back_b[back_b$Block != 4, colnames(back_b) != "Block"],
                  bg_test = back_b[back_b$Block == 4, colnames(back_b) != "Block"],
                  bg_block = back_b)
  } else {
    data1 <- list(all = data,
                  train = data_b[data_b$Block != 4, colnames(data_b) != "Block"],
                  test = data_b[data_b$Block == 4, colnames(data_b) != "Block"],
                  block = data_b)
  }

  return(data1)
}


