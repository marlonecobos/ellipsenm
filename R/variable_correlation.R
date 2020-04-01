#' Evaluates correlation among variables
#'
#' @description variable_correlation helps in evaluating correlation among
#' distinct variables.
#'
#' @param variables RasterStack, RasterBrick, or matrix. If matrix, columns
#' represent distinct variables for analysis, otherwise, group of raster layers.
#' @param sample_size (numeric) sample size to be taken from all variables;
#' default = 10000.
#' @param correlation_limit (numeric) absolute value of correlation limit;
#' default = 0.8.
#' @param save (logical) whether or not to save the results; default = FALSE.
#' @param name (character) name of the csv files to be writen;
#' default = "correlation".
#' @param corrplot (logical) whether or not to plot the results; default = FALSE.
#' @param magnify_to (numeric) optional value to be used to magnify all values
#' with absolute correlations above \code{correlation_limit}. Default = NULL.
#' @param ... other arguments to be passed to \code{\link[corrplot]{corrplot}}.
#' Arguments "type", "tl.col", and "tl.srt" are fixed.
#'
#' @return
#' A correlation matrix. If argument \code{corrplot} = TRUE correlation values
#' are shown in a plot.
#'
#' @usage
#' variable_correlation(variables, sample_size = 10000, correlation_limit = 0.8,
#'                      save = FALSE, name = "correlation", corrplot = FALSE,
#'                      magnify_to = NULL, ...)
#'
#' @details
#' If \code{magnify_to} is defined and \code{save} = TRUE, an additional csv
#' file named as "\code{name}_magnified.csv" will be written.
#'
#' @export
#'
#' @examples
#' # raster layers of environmental data
#' vars <- raster::stack(list.files(system.file("extdata", package = "ellipsenm"),
#'                                  pattern = "bio", full.names = TRUE))
#'
#' # simple correlation matrix
#' cors <- variable_correlation(variables, sample_size = 5000)
#'
#' # correlation matrix and plot (values correlated above |0.8| are magnified)
#' cors <- variable_correlation(variables, sample_size = 5000, corrplot = TRUE,
#'                              magnified = 2)
#'
#' # to save results check arguments "save" and "name"

variable_correlation <- function(variables, sample_size = 10000,
                                 correlation_limit = 0.8, save = FALSE,
                                 name = "correlation", corrplot = FALSE,
                                 magnify_to = NULL, ...) {
  # -----------
  # detecting potential errors
  if (missing(variables)) {
    stop("Argument 'variables' is necessary to perform the analysis")
  }
  var_class <- class(variables)[1]
  if (!var_class %in% c("matrix", "RasterStack", "RasterBrick")) {
    stop("'variables' must be either 'matrix', 'RasterStack', or 'RasterBrick'")
  }

  # -----------
  # preparing data
  if (var_class != "matrix") {
    ## getting data from the variables
    variables <- na.omit(raster::values(variables))

    ## sample of 10000 values if more pixels exist (optional)
    if (nrow(variables) > sample_size) {
      variables <- variables[sample(1:nrow(variables), sample_size), ]
    }
  }

  # -----------
  # analyses
  ## correlation matrix calculation
  correlation_matrix <- cor(variables)

  ## detecting correlated varaibles more easily
  correlation_matrix1 <- correlation_matrix # making other table with results

  max_cor <- correlation_limit # maximum value of correlation allowed

  if (!is.null(magnify_to)) {
    mv <- magnify_to
    for (i in 1:dim(correlation_matrix1)[2]) { #correlated values will turn into 2 for easier detection
      for (j in 1:dim(correlation_matrix1)[1]) {
        correlation_matrix1[j, i] <- ifelse(correlation_matrix1[j, i] < -max_cor,
                                            -mv, correlation_matrix1[j, i])
        correlation_matrix1[j, i] <- ifelse(correlation_matrix1[j, i] > max_cor,
                                            mv, correlation_matrix1[j, i])
      }
    }
  }

  # -----------
  # write and plot
  # saving correlation matrix
  if (save == TRUE) {
    write.csv(correlation_matrix, paste0(name, ".csv"), row.names = TRUE)

    if (!is.null(magnify_to)) {
      name <- paste0(name, "_magnified.csv")
      write.csv(correlation_matrix1, name, row.names = TRUE)
    }
  }

  # plotting
  if (corrplot == TRUE) {
    if (!is.null(magnify_to)) {
      cor_mat <- correlation_matrix1/mv
    } else {
      cor_mat <- correlation_matrix
    }
    corrplot::corrplot(cor_mat, type = "upper", tl.col = "black",
                       tl.srt = 45, ...)
  }

  # -----------
  # return results
  return(correlation_matrix)
}
