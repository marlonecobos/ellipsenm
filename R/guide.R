#' Creation of HTML reports from ellipsenm main processes
#'
#' @description report creates an HTML file that summarizes all outputs from
#' the main processes that the ellipsenm package performs.
#'
#' @param guide_type (character) type of report to be produced. Options are:
#' "calibration", "enm", "whole_process", "overlap".
#' @param output_directory (character) name of the directory where the HTML file
#' will be produced. Default = working directory.
#'
#' @return
#' An Rmarkdown file that serves as a guide to perform any of the three processes
#' defined in \code{guide_type}. The Rmarkdown will be automatically open after
#' running the function and all processes can be performed there.
#'
#' @export

guide <- function(guide_type, output_directory = getwd()) {
  # -----------
  # detecting potential errors
  if (missing(guide_type)) {
    stop("Argument guide_type is missing.")
  }
  if (!guide_type %in% c("calibration", "enm", "whole_process", "overlap")) {
    stop("Argument guide_type is not valid, please see function's help.")
  }

  # -----------
  # preparing report
  # report types for enm
  if (guide_type == "enm") {
    name <- paste0(output_directory, "/enm_guide.Rmd")
    suppressMessages(
      file.copy(from = system.file("extdata", "enm_guide.Rmd",
                                   package = "ellipsenm"), to = name)
    )
  }

  # report for calibration
  if (guide_type == "calibration") {
    name <- paste0(output_directory, "/calibration_guide.Rmd")
    suppressMessages(
      file.copy(from = system.file("extdata", "calibration_guide.Rmd",
                                   package = "ellipsenm"), to = name)
    )
  }

  # report for the whole process
  if (guide_type == "whole_process") {
    name <- paste0(output_directory, "/whole_process_guide.Rmd")
    suppressMessages(
      file.copy(from = system.file("extdata", "whole_process_guide.Rmd",
                                   package = "ellipsenm"), to = name)
    )
  }

  # report type for overlap
  if (guide_type == "overlap") {
    name <- paste0(output_directory, "/niche_overlap_guide.Rmd")
    suppressMessages(
      file.copy(from = system.file("extdata", "niche_overlap_guide.Rmd",
                                   package = "ellipsenm"), to = name)
    )
  }

  # -----------
  # reporting
  file.edit(name)

  cat(paste0("\nA file named ", name, " has been produced.\n",
             "Check your output directory in:\t", getwd()))
}
