#' Omission rates calculation for single models
#'
#' @description omission_rate calculates omission rates of predictions of
#' ecological niche models based on one or multiple user-specified thresholds.
#'
#' @param prdiction RasterLayer of the continuous prediction to be evaluated.
#' @param ellipsenm_result (logical) whether or not the prediction resulted from
#' the \code{\link[ellipsenm]} package.
#' @param threshold (numeric vector) value(s) from 0 to 100 that will be used as
#' thresholds; default = 5.
#' @param train_occurrences a numerical matrix containing coordinates of the
#' occurrence data used to create the ecological niche model to be evaluated;
#' columns must be: longitude and latitude.
#' @param test_occurrences a numerical matrix containing coordinates of the
#' occurrences used to test the ecological niche model to be evaluated; columns
#' must be: longitude and latitude.
#'
#' @return A named numeric numeric vector with the results.
#'
#' @export
#'
#' @examples
#' # single threshold
#' model <- raster::raster(list.files(system.file("extdata", package = "kuenm"),
#'                                    pattern = "sp_model.tif", full.names = TRUE))
#' thres <- 5
#' octr <- read.csv(list.files(system.file("extdata", package = "kuenm"),
#'                             pattern = "sp_train.csv", full.names = TRUE))
#' octe <- read.csv(list.files(system.file("extdata", package = "kuenm"),
#'                             pattern = "sp_test.csv", full.names = TRUE))
#'
#' om_rate <- omission_rate(model, threshold = thres,
#'                          occ.tra = octr, occ.test = octe)
#'
#' # multiple thresholds
#' thres1 <- c(5, 10, 20)
#'
#' om_rate <- omission_rate(model, threshold = thres1,
#'                          occ.tra = octr, occ.test = octe)

omission_rate <- function(prdiction, ellipsenm_result, threshold = 5,
                          train_occurrences = NULL, test_occurrences = NULL) {

  if(model@data@min == model@data@max) {
    warning("\nModel imput has no variability, omission rate will return NA.\n")

    om_rate <- NA
  }else {

    suit_val_cal <- na.omit(raster::extract(model, occ.tra))
    suit_val_eval <- na.omit(raster::extract(model, occ.test))

    om_rate <- vector()

    for (i in 1:length(threshold)) {
      val <- ceiling(length(occ.tra[, 1]) * threshold[i] / 100) + 1
      omi_val_suit <- sort(suit_val_cal)[val]
      om_rate[i] <- as.numeric(length(suit_val_eval[suit_val_eval < omi_val_suit]) / length(suit_val_eval))
    }

    names(om_rate) <- paste("om_rate_", threshold, "%", sep = "")
    return(om_rate)
  }
}
