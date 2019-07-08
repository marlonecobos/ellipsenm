#' Helper function to format function reports
#'
#' @param col (character) hexadecimal code for the main color to be used in the
#' html report. Default = "#0000FF"; blue.
#' @param name (character) name of the CSS file to be written. File format does
#' not need to be added. Default = "eenm_report_format".
#'
#' @return CSS file containing HTML formatting instructions.
#'
#' @export
#'
#' @examples
#' # example with green color
#' report_format(col = "#15AB0A")

report_format <- function(col = "#0000FF", name = "eenm_report_format") {
  # html style
  sink(paste0(name, ".css"))
  cat("#main .nav-pills > li.active > a,\n#main .nav-pills > li.active > a:hover,
#main .nav-pills > li.active > a:focus {\n\tbackground-color: ")
  cat(col)
  cat(";\n}\n\n#main .nav-pills > li > a:hover {\n\tbackground-color: ")
  cat(col)
  cat(";\n}\n\nh1,h2,h3,h4,h5,h6,legend{\n\tcolor: ")
  cat(col)
  cat(";\n}\n\n#nav-top span.glyphicon {\n\tcolor: ")
  cat(col)
  cat(";\n}\n\n#table-of-contents header{\n\tcolor: ")
  cat(col)
  cat(";\n}\n\n#table-of-contents h2{\n\tbackground-color: ")
  cat(col)
  cat(";\n}\n\n#main a {\n\tbackground-image: linear-gradient(180deg,#d64a70,#d64a70);
\tcolor:#c7254e;\n}\n\na:hover{\n\tcolor:#3d1308\n}\n\na:visited{\n\tcolor:#3d1308\n}")
  sink()
}
