#' Creation of HTML reports from ellipsenm main processes
#'
#' @description report creates an HTML file that summarizes all outputs from
#' the main processes that the ellipsenm package performs.
#'
#' @param repor_type (character) type of report to be produced. Options are:
#' "calibration", "enm_suitability", "enm_mahalanobis", "enm_both", and "overlap".
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
  if (report_type %in% c("calibration", "enm_suitability", "enm_mahalanobis",
                       "enm_both", "overlap")) {

    if (report_type == "calibration") {
      suppressMessages(
        file.copy(from = system.file("extdata", "enm_calibration_report.Rmd",
                                     package = "ellipsenm"), to = name)
      )
    }

    if (report_type == "enm_suitability") {
      suppressMessages(
        file.copy(from = system.file("extdata", "enm_suitability_report.Rmd",
                                     package = "ellipsenm"), to = name)
      )

    }

    if (report_type == "enm_mahalanobis") {
      suppressMessages(
        file.copy(from = system.file("extdata", "enm_mahalanobis_report.Rmd",
                                     package = "ellipsenm"), to = name)
      )
    }

    if (report_type == "enm_both") {
      suppressMessages(
        file.copy(from = system.file("extdata", "enm_both_report.Rmd",
                                     package = "ellipsenm"), to = name)
      )
    }

    if (report_type == "overlap") {
      suppressMessages(
        file.copy(from = system.file("extdata", "niche_overlap_report.Rmd",
                                     package = "ellipsenm"), to = name)
      )
    }
  } else {
    stop("Argument report_type is not valid, please see function's help.")
  }

  # -----------
  # rendering
  suppressWarnings(rmarkdown::render(paste0(name, ".Rmd"), "all", quiet = TRUE))

  # -----------
  # reporting
  cat(paste0("\nA file named", paste0(name, ".html"), "has been produced.\n",
             "Check your output directory in:\t", getwd()))
}
