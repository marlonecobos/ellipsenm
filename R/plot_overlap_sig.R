#' Plot results of significance tests for niche overlap analyses
#'
#' @description plot_overlap_sig helps to create a plot to interpret results
#' from statistical significance tests in overlap analyses.
#'
#' @param object overlap_ellipsoid object resulted from using the function
#' \code{\link{ellipsoid_overlap}}.
#' @param niches (numeric) pair of integer numbers denoting the niches for
#' which results will be plotted. Default = c(1, 2).
#' @param main (character) plot main title, default = "".
#' @param xlab (character) label of x axis, default = "Overlap".
#' @param ylab (character) label of y axis, default = "Frequency".
#' @param breaks type of breaks as in \code{\link[graphics]{hist}}, default = 30.
#' @param bar_col color of histogram bars, default = "gray65".
#' @param cl_col color of confidence limit line, default = "darkgreen".
#' @param observed_col color of line representing observed value of overlap,
#' default = "blue".
#' @param cl_lty type of line for confidence limit, default = 2.
#' @param observed_lty type of line for observed value of overlap, default = 1.
#' @param cl_lwd thickness of line for confidence limit, default = 2.
#' @param observed_lwd thickness of line for observed value of overlap,
#' default = 2.
#' @param xlim limits of x axis, the default, NULL, derives limits from data.
#'
#' @return A plot with results from tests of statistical significance.
#'
#' @usage
#' plot_overlap_sig(object, niches = c(1, 2), main = "",
#'                  xlab = "Overlap", ylab = "Frequency",
#'                  breaks = 30, bar_col = "gray65",
#'                  cl_col = "darkgreen", observed_col = "blue",
#'                  cl_lty = 2, observed_lty = 1, cl_lwd = 2,
#'                  observed_lwd = 2, xlim = NULL)
#'
#' @export

plot_overlap_sig <- function(object, niches = c(1, 2), main = "",
                             xlab = "Overlap", ylab = "Frequency",
                             breaks = 30, bar_col = "gray65",
                             cl_col = "darkgreen", observed_col = "blue",
                             cl_lty = 2, observed_lty = 1, cl_lwd = 2,
                             observed_lwd = 2, xlim = NULL) {

  # tests
  if (missing(object)) {
    stop("Argument 'object' is necessary to perform the analysis.")
  }

  # preparing data
  ## niches
  rown <- paste0("Niche_", paste0(niches, collapse = "_vs_"))

  ## data to plot
  clp <- object@union_overlap[rown, "pre_defined_CL"]
  null_val <- object@significance_results$union_random[[rown]]$overlap
  cl <- quantile(null_val, clp)
  obs_val <- object@union_overlap[rown, "overlap"]

  # xlim if not defined
  if (is.null(xlim)) {
    xlim <- range(c(null_val, obs_val))
  }

  # the plot
  hist(null_val, breaks = breaks, main = main, xlab = xlab, ylab = ylab,
       xlim = xlim)
  abline(v = cl, col = cl_col, lwd = cl_lwd, lty = cl_lty)
  abline(v = obs_val, col = observed_col,
         lwd = observed_lwd, lty = observed_lty)
  box(bty = "l")
}

