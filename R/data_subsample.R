#' Helper function to subsample occurrence data
#'
#' @param data data.frame of occurrence records. Columns must be longitude and
#' latitude.
#' @param replicates (numeric) number of subsamples to be produced. If
#' \code{replicate_type} = "jackknife" the maximum number of subsamples produced
#' equal the number of records in \code{data}. Default = 1 returns the complete
#' set of data as the only element of a list.
#' @param replicate_type (character) name of the sampling method to be used.
#' Options are "bootstrap" or "jackknife". Default = "bootstrap".
#' @param bootstrap_percentage (numeric) percentage of data to be bootstrapped
#' for each replicate. Default = 50. Valid if \code{replicates} > 1 and
#' \code{replicate_type} = "bootstrap".
#'
#' @export

data_subsample <- function(data, replicates = 1, replicate_type = "bootstrap",
                           bootstrap_percentage = 50) {

  if (replicates == 1) {
    samples <- list(data)

  } else {
    ndata <- nrow(data)

    if (replicate_type == "bootstrap" | replicate_type == "jackknife") {

      if (replicate_type == "bootstrap") {
        part <- round(ndata * bootstrap_percentage / 100)
        samples <- lapply(1:replicates, function(x) {
          set.seed(x)
          data[sample(ndata, part), ]
        })

      } else {
        if (replicates > ndata) {
          warning(paste0("Maximum number of replicates under \"jackknife\" replicate_type",
                         "\ncannot be > ", ndata, ", only ", ndata, " replicates will be produced."))
          replicates <- ndata
        }

        samples <- lapply(1:replicates, function(x) {
          data[-x, ]
        })
      }

    } else {
      stop("Argument replicate_type needs to be either \"bootstrap\" or \"jackknife\".")
    }
  }

  return(samples)
}
