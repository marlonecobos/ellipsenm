#' Helper function to find raster extention
#' @param format (character) any of the format types allowed for raster objects.
#' See \code{\link[raster]{writeFormats}}
#' @export
#' @return Raster extension according to format type.

rformat_type <- function(format) {
  if (missing(format)) {stop("Argument format needs to be defined.")}
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

#' Helperf function to get and write ellipsoid metadata
#' @param ellipsoid object of class ellipsoid*.
#' @param name (character) name of the file to be written. If the object in
#' \code{ellipsoid} is replicated, number of replicates and mean, min, and max,
#' are added as prefixes to each file. The suffix "ell_meta.txt" will be added
#' to \code{name}. Default = "".
#' @export
#' @return
#' A summary of ellipsoid metadata as a data.frame.
#'
#' Writes a csv file with the summary of ellipsoid metadata named
#' "*_metadata_summary.csv" and txt files with the complete ellipsoid
#' metadata in the working directory.

write_ellmeta <- function(ellipsoid, name = "") {
  if (!missing(ellipsoid)) {
    cls <- class(ellipsoid)[1]
    if (!cls %in% c("ellipsoid", "ellipsoid_model_sim", "ellipsoid_model_rep")) {
      stop("Argument ellipsoid must be of class ellipsoid*.")
    }
  } else {
    stop("Argument ellipsoid is necessary to perform the analysis.")
  }

  name <- gsub("\\\\", "/", name)
  name <- unlist(strsplit(name, "/"))
  ndir <- paste0(paste(name[-length(name)], collapse = "/"), "/")
  namesum <- paste0(ndir, name[length(name)], "_metadata_summary.csv")
  name <- paste0(name[length(name)], "_ell_meta.txt")

  if (cls %in% c("ellipsoid", "ellipsoid_model_sim")) {
    namesim <- name
    name <- paste0(ndir, name)
    cat("Ellipsoid_metadata\n", file = name, append = TRUE, sep = "")
    cat("\nMethod:\t", ellipsoid@method, file = name, append = TRUE, sep = "")
    cat("\n\nLevel:\t", ellipsoid@level, file = name, append = TRUE, sep = "")
    cat("\n\nCentroid:\n", paste0(names(ellipsoid@centroid), "\t",
                                ellipsoid@centroid, "\n"),
        file = name, append = TRUE, sep = "")
    cat("\nCovariance_matrix:\n",
        paste0(c("", colnames(ellipsoid@covariance_matrix)), collapse = "\t"),
        "\n", file = name, append = TRUE, sep = "")
    suppressWarnings(write.table(ellipsoid@covariance_matrix, sep = "\t", file = name,
                                 append = TRUE, quote = FALSE, col.names = FALSE))
    cat("\nVolume:\t", ellipsoid@niche_volume, file = name, append = TRUE)
    cat("\n\nSemi-axes_length:\n", paste0(names(ellipsoid@semi_axes_length), "\t",
                                  ellipsoid@semi_axes_length, "\n"),
        file = name, append = TRUE, sep = "")
    cat("\nAxes_coordinates:\n", file = name, append = TRUE, sep = "")
    a_cord <- ellipsoid@axes_coordinates
    cords <- lapply(1:length(a_cord), function(x) {
      cat(letters[x], "\n", file = name, append = TRUE, sep = "")
      cat(paste0(c("", colnames(a_cord[[x]])), collapse = "\t"),
          "\n", file = name, append = TRUE, sep = "")
      suppressWarnings(write.table(a_cord[[x]], sep = "\t", file = name, append = TRUE,
                                   quote = FALSE, col.names = FALSE))
    })

    ell_meta <- data.frame(ellipsoid_model = c(ellipsoid@method, ellipsoid@level,
                                               round(ellipsoid@niche_volume, 2),
                                               namesim))

  } else {
    ellipsoid <- ellipsoid@ellipsoids
    nam <- names(ellipsoid)
    if (is.null(nam)) {
      enames <- as.character(1:length(ellipsoid))
      enames1 <- paste0("ellipsoid", enames)
    } else {
      if (length(grep("replicate", nam)) > 0 & length(grep("mean", nam)) > 0) {
        enames <- c(1:(length(nam) - 3), "mean", "min", "max")
        enames1 <- nam
      } else {
        enames <- nam
        enames1 <- nam
      }
    }
    namesim <- paste0(enames, "_", name)
    name <- paste0(ndir, enames, "_", name)
    ell_meta <- lapply(1:length(name), function(x) {
      cat("Ellipsoid_metadata\n", file = name[x], append = TRUE, sep = "")
      cat("\nMethod:\t", ellipsoid[[x]]@method, file = name[x], append = TRUE, sep = "")
      cat("\n\nLevel:\t", ellipsoid[[x]]@level, file = name[x], append = TRUE, sep = "")
      cat("\n\nCentroid:\n", paste0(names(ellipsoid[[x]]@centroid), "\t",
                                    ellipsoid[[x]]@centroid, "\n"),
          file = name[x], append = TRUE, sep = "")
      cat("\nCovariance_matrix:\n",
          paste0(c("", colnames(ellipsoid[[x]]@covariance_matrix)), collapse = "\t"),
          "\n", file = name[x], append = TRUE, sep = "")
      suppressWarnings(write.table(ellipsoid[[x]]@covariance_matrix, sep = "\t", file = name[x],
                                   append = TRUE, quote = FALSE, col.names = FALSE))
      cat("\nVolume:\t", ellipsoid[[x]]@niche_volume, file = name[x], append = TRUE)
      cat("\n\nSemi-axes_length:\n", paste0(names(ellipsoid[[x]]@semi_axes_length), "\t",
                                            ellipsoid[[x]]@semi_axes_length, "\n"),
          file = name[x], append = TRUE, sep = "")
      cat("\nAxes_coordinates:\n", file = name[x], append = TRUE, sep = "")
      a_cord <- ellipsoid[[x]]@axes_coordinates
      cords <- lapply(1:length(a_cord), function(y) {
        cat(letters[y], "\n", file = name[x], append = TRUE, sep = "")
        cat(paste0(c("", colnames(a_cord[[y]])), collapse = "\t"),
            "\n", file = name[x], append = TRUE, sep = "")
        suppressWarnings(write.table(a_cord[[y]], sep = "\t", file = name[x], append = TRUE,
                                     quote = FALSE, col.names = FALSE))
      })
      ellm <- c(ellipsoid[[x]]@method, ellipsoid[[x]]@level,
                round(ellipsoid[[x]]@niche_volume, digits = 2))
      return(ellm)
    })
    ell_meta <- as.data.frame(rbind(do.call(cbind, ell_meta), namesim))
    colnames(ell_meta) <- enames1
  }
  row.names(ell_meta) <- c("Method", "Level", "Volume", "Other_metadata")
  write.csv(ell_meta, namesum, row.names = TRUE)
  return(ell_meta)
}

