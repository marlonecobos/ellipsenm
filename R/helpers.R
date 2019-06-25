#' Helper function to find raster extention
#' @param format (character) any of the format types allowed for raster objects.
#' See \code{\link[raster]{writeFormats}}
#' @export
#' @return Raster extension according to format type.

rformat_type <- function(format) {
  if (missing(format)) {
    stop("Argument format needs to be defined.")
  }

  if (format == "raster") {format1 <- ".grd"}
  if (format == "GTiff") {format1 <- ".tif"}
  if (format == "EHdr") {format1 <- ".bil"}
  if (format == "ascii") {format1 <- ".asc"}
  if (format == "SAGA") {format1 <- ".sdat"}
  if (format == "IDRISI") {format1 <- ".rst"}
  if (format == "CDF") {format1 <- ".nc"}
  if (format == "ENVI") {format1 <- ".envi"}
  if (format == "HFA") {format1 <- ".img"}

  return(format1)
}


#' Helperf function to write ellipsoid metadata
#' @param ellipsoid object of class ellipsoid*.
#' @param name (character) name of the file to be written. If the object in
#' \code{ellipsoid} is replicated, number of replicate and mean, min, and max,
#' are added as prefixes to each file. Default = "ell_meta.txt"
#' @export
#' @return Writes a txt file with the ellipsoid metadata in the working directory.

write_ellmeta <- function(ellipsoid, name = "ell_meta.txt") {
  if (!missing(ellipsoid)) {
    cls <- class(ellipsoid)[1]
    if (!cls %in%
        c("ellipsoid", "ellipsoid_model_sim", "ellipsoid_model_rep")) {
      stop("Argument ellipsoid must be of class ellipsoid*.")
    }
  }else {
    stop("Argument ellipsoid is necessary to perform the analysis.")
  }

  name <- gsub("\\\\", "/", name)
  name <- unlist(strsplit(name, "/"))
  ndir <- paste0(paste(name[-length(name)], collapse = "/"), "/")
  name <- name[length(name)]

  if (cls %in% c("ellipsoid", "ellipsoid_model_sim")) {
    name <- paste0(ndir, name)
    cat("Ellipsoid metadata\n", file = name, append = TRUE, sep = "")
    cat("\nMethod:\t", ellipsoid@method, file = name, append = TRUE, sep = "")
    cat("\n\nLevel:\t", ellipsoid@level, file = name, append = TRUE, sep = "")
    cat("\n\nCentroid:\n", paste0(names(ellipsoid@centroid), "\t",
                                ellipsoid@centroid, "\n"),
        file = name, append = TRUE, sep = "")
    cat("\nCovariance matrix:\n",
        paste0(c("  ", colnames(ellipsoid@covariance_matrix)), collapse = "\t"),
        "\n", file = name, append = TRUE, sep = "")
    suppressWarnings(write.table(ellipsoid@covariance_matrix, sep = "\t", file = name,
                                 append = TRUE, quote = FALSE, col.names = FALSE))
    cat("\nVolume:\t", ellipsoid@niche_volume, file = name, append = TRUE)
    cat("\n\nSemi-axes length:\n", paste0(names(ellipsoid@semi_axes_length), "\t",
                                  ellipsoid@semi_axes_length, "\n"),
        file = name, append = TRUE, sep = "")
    cat("\nAxes coordinates:\n", file = name, append = TRUE, sep = "")
    a_cord <- ellipsoid@axes_coordinates
    cords <- lapply(1:length(a_cord), function(x) {
      cat(letters[x], "\n", file = name, append = TRUE, sep = "")
      cat(paste0(c("  ", colnames(a_cord[[x]])), collapse = "\t"),
          "\n", file = name, append = TRUE, sep = "")
      suppressWarnings(write.table(a_cord[[x]], sep = "\t", file = name, append = TRUE,
                                   quote = FALSE, col.names = FALSE))
    })

  } else {
    ellipsoid <- ellipsoid@ellipsoids
    enames <- c(1:(length(ellipsoid) - 3), "mean", "min", "max")
    name <- paste0(ndir, enames, "_", name)
    rep_meta <- lapply(1:length(name), function(x) {
      cat("Ellipsoid metadata\n", file = name[x], append = TRUE, sep = "")
      cat("\nMethod:\t", ellipsoid[[x]]@method, file = name[x], append = TRUE, sep = "")
      cat("\n\nLevel:\t", ellipsoid[[x]]@level, file = name[x], append = TRUE, sep = "")
      cat("\n\nCentroid:\n", paste0(names(ellipsoid[[x]]@centroid), "\t",
                                    ellipsoid[[x]]@centroid, "\n"),
          file = name[x], append = TRUE, sep = "")
      cat("\nCovariance matrix:\n",
          paste0(c("  ", colnames(ellipsoid[[x]]@covariance_matrix)), collapse = "\t"),
          "\n", file = name[x], append = TRUE, sep = "")
      suppressWarnings(write.table(ellipsoid[[x]]@covariance_matrix, sep = "\t", file = name[x],
                                   append = TRUE, quote = FALSE, col.names = FALSE))
      cat("\nVolume:\t", ellipsoid[[x]]@niche_volume, file = name[x], append = TRUE)
      cat("\n\nSemi-axes length:\n", paste0(names(ellipsoid[[x]]@semi_axes_length), "\t",
                                            ellipsoid[[x]]@semi_axes_length, "\n"),
          file = name[x], append = TRUE, sep = "")
      cat("\nAxes coordinates:\n", file = name[x], append = TRUE, sep = "")
      a_cord <- ellipsoid[[x]]@axes_coordinates
      cords <- lapply(1:length(a_cord), function(y) {
        cat(letters[y], "\n", file = name[x], append = TRUE, sep = "")
        cat(paste0(c("  ", colnames(a_cord[[y]])), collapse = "\t"),
            "\n", file = name[x], append = TRUE, sep = "")
        suppressWarnings(write.table(a_cord[[y]], sep = "\t", file = name[x], append = TRUE,
                                     quote = FALSE, col.names = FALSE))
      })
    })
  }
}
