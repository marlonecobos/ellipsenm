#' Split occurrence data into training and testing data
#'
#' @description split_data splits occurrences into training and testing data
#' based on distinct methods.
#'
#' @param data data.frame of occurrence records containing at least longitude
#' and latitude columns.
#' @param train_proportion (numeric) proportion (from 0 to 1) of data to be used
#' as training occurrences. The remaining data will be used for testing.
#' Default = 0.5.
#' @param method (character) method for selecting training and testing data.
#' Options are: "random" and "block"; default = "random".
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
#' @export
#'
#' @examples
#' occurrences <- read.csv(system.file("extdata", "occurrences_thin.csv",
#'                                     package = "kuenm"))
#'
#' # random split 50% for trainig and 50% for testing
#' data_split <- split_data(data = occs, train.proportion = 0.5)

split_data <- function(data, train_proportion = 0.5, method = "random",
                       save = FALSE, name = "occurrences") {

  # -----------
  # detecting potential errors
  if (missing(data)) {
    stop("Argument data is not defined.")
  }

  # -----------
  # processing
  occ <- na.omit(data)

  if (method == "random") {
    files <- occ_randsplit(occ, train_proportion)
  }

  if (method == "block") {
    files <- occ_blocksplit(occ, train_proportion)
  }

  # -----------
  # writing results
  if (save == TRUE) {
    names <- paste0(name, c("_all", "_train", "_test"), ".csv")
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



