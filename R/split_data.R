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
#'
#' @export
#'
#' @examples
#' occurrences <- read.csv(system.file("extdata", "occurrences_comp.csv",
#'                                     package = "ellipsenm"))
#'
#' # random split 50% for trainig and 50% for testing
#' data_split <- split_data(data = occs, train.proportion = 0.5)

split_data <- function(data, method = "random", longitude, latitude,
                       train_proportion, save = FALSE, name = "occurrences") {

  # -----------
  # detecting potential errors and defining missing proportion
  if (missing(data)) {
    stop("Argument data is not defined.")
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

    files <- occ_blocksplit(occ, longitude, latitude, train_proportion)
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
  ids <- sample(ndata, size = round(train.proportion * ndata))
  data1 <- list(all = data, train = data[ids, ], test = data[-ids, ])

  return(data1)
}



#' Split occurrences in training and testing data using blocks
#'
#' @description occ_blocksplit splits a set of occurrences to obtain training and
#' testing data randomly.
#'
#' @param data matrix or data.frame with the occurrences to be split. Columns
#' may vary but species, longitude, and latitue are recommended.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @param train_proportion (numeric) proportion of data to be used as training
#' occurrences. Available options are 0.25, 0.5, and 0.75. The remaining data
#' will be used for testing. Default = 0.75.
#'
#' @export
#'
#' @return
#' List with all occurrences (all), training occurrences (train), and testing
#' (test) occurrences.



occ_blocksplit <- function(data, longitude, latitude, train_proportion = 0.75) {
  ndata <- nrow(data)

  n1 <- ceiling(nrow(data) / latitude)
  n2 <- floor(nrow(data) / latitude)
  n3 <- ceiling(n1 / latitude)
  n4 <- ceiling(n2 / latitude)
  grpA <- data[order(data[, latitude]), ][1:n1, ]
  grpB <- data[rev(order(data[, latitude])), ][1:n2, ]
  grp1 <- grpA[order(grpA[, longitude]), ][1:(n3), ]
  grp2 <- grpA[!rownames(grpA) %in% rownames(grp1), ]
  grp3 <- grpB[order(grpB[, longitude]), ][1:(n4), ]
  grp4 <- grpB[!rownames(grpB) %in% rownames(grp3), ]

  if (nrow(grp1) > 0) grp1$grp <- 1
  if (nrow(grp2) > 0) grp2$grp <- 2
  if (nrow(grp3) > 0) grp3$grp <- 3
  if (nrow(grp4) > 0) grp4$grp <- 4

  data_b <- rbind(grp1, grp2, grp3, grp4)
  colnames(data_b) <- c(colnames(data), "Block")

  data1 <- list(all = data,
                train = data_b[data_b$Block != 4, colnames(data_b) != "Block"],
                test = data_b[data_b$Block == 4, colnames(data_b) != "Block"],
                block = data_b)

  return(data1)
}


