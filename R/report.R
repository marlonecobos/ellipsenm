#' Creation of HTML reports from ellipsenm main processes
#'
#' @description report creates an HTML file that summarizes all outputs from
#' the main processes that the ellipsenm package performs.
#'
#' @param repor_type (character) type of report to be produced. Options are:
#' "enm", "calibration", and "overlap".
#' @param name (character) name of the HTML file to be produced.
#'
#' @return
#' An HTML file summarizing results from the main processes that the ellipsenm
#' package performs.
#'
#' @details This function is used along with the \code{\link{ellipsoid_model}},
#' ellipsoid_calibration, and ellipsoid_overlap functions.
#'
#' @export

report <- function(report_type, name = "results_report") {
  # -----------
  # detecting potential errors, other potential problems tested in code
  if (missing(report_type)) {
    stop("Argument report_type is missing.")
  }

  # -----------
  # preparing report
  if (report_type == "enm" | report_type == "calibration" | report_type == "overlap") {
    if (report_type == "enm") {
      enm_report(name = name)
    }

    if (report_type == "calibration") {
      calibration_report(name = name)
    }

    if (report_type == "overlap") {
      overlap_report(name = name)
    }
  } else {
    stop("Argument report_type is not valid, please see function's help.")
  }

  # -----------
  # rendering
  rmarkdown::render(paste0(name, ".Rmd"), "html_document", quiet = TRUE)
  unlink(paste0(name, ".Rmd"))

  # -----------
  # reporting
  cat(paste0("A file named", paste0(name, ".html"), "has been produced.\n",
             "Check your working directory:\t", getwd()))
}
