#' Creation of HTML reports from ellipsenm main processes
#'
#' @description report creates an HTML file that summarizes all outputs from
#' the main processes that the ellipsenm package performs.
#'
#' @param repor_type (character) type of report to be produced. Options are:
#' "calibration", "enm", and "overlap".
#' @param prediction (character) if report_type = "enm". Type of prediction that
#' was performed, oprions are: "suitability", "mahalanobis", and "both".
#' @param projected (logical) if report_type = "enm", whether or not the model
#' was projected to distinct scenarios.
#' @param output_directory (character) name of the directory where the HTML file
#' will be produced. Default = working directory.
#'
#' @return
#' An HTML file summarizing results from the main processes that the ellipsenm
#' package performs. Some 3D plots from the rgl package may appear.
#'
#' @details This function is used along with the \code{\link{ellipsoid_model}},
#' ellipsoid_calibration, and ellipsoid_overlap functions. To obtain nice reports
#' a css file with format instructions should exist. See \code{\link{report_format}}.
#'
#' @export

report <- function(report_type, prediction, projected, output_directory = getwd()) {
  # -----------
  # detecting potential errors
  if (missing(report_type)) {
    stop("Argument report_type is missing.")
  }
  if (!report_type %in% c("calibration", "enm", "overlap")) {
    stop("Argument report_type is not valid, please see function's help.")
  }
  if (report_type == "enm") {
    if (missing(prediction)) {
      stop("Argument prediction is missing.")
    }
    if (!prediction %in% c("mahalanobis", "both", "suitability")) {
      stop("Argument prediction is not valid, please see function's help.")
    }
    if (missing(projected)) {
      stop("Argument projected is missing.")
    }
  }

  # -----------
  # preparing report
  # report types for enm
  if (report_type == "enm") {
    if (projected == TRUE) {
      if (prediction != "mahalanobis") {
        if (prediction == "both") {
          name <- paste0(output_directory, "/enm_both_pr_report.Rmd")
          suppressMessages(
            file.copy(from = system.file("extdata", "enm_both_pr_report.Rmd",
                                         package = "ellipsenm"), to = name)
          )
        } else {
          name <- paste0(output_directory, "/enm_suitability_pr_report.Rmd")
          suppressMessages(
            file.copy(from = system.file("extdata", "enm_suitability_pr_report.Rmd",
                                         package = "ellipsenm"), to = name)
          )
        }
      } else {
        name <- paste0(output_directory, "/enm_mahalanobis_pr_report.Rmd")
        suppressMessages(
          file.copy(from = system.file("extdata", "enm_mahalanobis_pr_report.Rmd",
                                       package = "ellipsenm"), to = name)
        )
      }
    } else {
      if (prediction != "mahalanobis") {
        if (prediction == "both") {
          name <- paste0(output_directory, "/enm_both_report.Rmd")
          suppressMessages(
            file.copy(from = system.file("extdata", "enm_both_report.Rmd",
                                         package = "ellipsenm"), to = name)
          )
        } else {
          name <- paste0(output_directory, "/enm_suitability_report.Rmd")
          suppressMessages(
            file.copy(from = system.file("extdata", "enm_suitability_report.Rmd",
                                         package = "ellipsenm"), to = name)
          )
        }
      } else {
        name <- paste0(output_directory, "/enm_mahalanobis_report.Rmd")
        suppressMessages(
          file.copy(from = system.file("extdata", "enm_mahalanobis_report.Rmd",
                                       package = "ellipsenm"), to = name)
        )
      }
    }
  }

  # report for calibration
  if (report_type == "calibration") {
    name <- paste0(output_directory, "/enm_calibration_report.Rmd")
    suppressMessages(
      file.copy(from = system.file("extdata", "enm_calibration_report.Rmd",
                                   package = "ellipsenm"), to = name)
    )
  }

  # report type for overlap
  if (report_type == "overlap") {
    name <- paste0(output_directory, "/niche_overlap_report.Rmd")
    suppressMessages(
      file.copy(from = system.file("extdata", "niche_overlap_report.Rmd",
                                   package = "ellipsenm"), to = name)
    )
  }

  # -----------
  # rendering
  suppressWarnings(rmarkdown::render(name, "all", quiet = TRUE))

  # -----------
  # reporting
  cat(paste0("\nA file named ", paste0(gsub(".Rmd$", "", name), ".html"),
             " has been produced.\n", "Check your output directory in:\t", getwd()))
}
