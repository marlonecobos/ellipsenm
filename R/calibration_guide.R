#' A guide for ellipsoid ENM calibration
#'
#' @description calibration_guide generates an Rmarkdown file with information,
#' instructions, and details, to guide users through the process of model
#' calibration.
#'
#' @param data data.frame of occurrence records. Columns must be: species,
#' longitude, and latitude. Optionally, if \code{raster_layers} is not defined,
#' \code{data} must include more columns containing the values of at least two
#' variables to be used for fitting ellipsoid* models.
#' @param species (character) name of the column with the name of the species.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @param raster_layers RasterStack of at least two environmental variables to be
#' extracted using geographic coordinates present in \code{data}. If not provided
#' data must include additional columns containing values of variables to fit
#' ellipsoid* models.
#' @param accessible_area SpatialPolygon* or RasterLayer to be used in masking
#' \code{raster_layers} according to a hypothesis of areas that have been
#' accessible to the species during a relevant period of time. See details.
#' @param output_directory name of the folder were all results will be written.
#' This helps to organize data and results and avoid confussion.
#'
#' @return
#' A folder containing an Rmarkdown file that serves as a guide to perform the
#' process of model calibration. The Rmarkdown will be automatically open after
#' running the function and all processes can be performed there.
#'
#' @export
#'
#' @details
#' \code{accessible_area}: if the user does not have a hypothesis of this area,
#' different options to create one exist. Four potential options to define this
#' area are: (1) a manual definition of accessible areas for the species; (2) the
#' use of geographic regions potentially relevant (e.g., Ecorregions, Biomes,
#' basins), generally they require further processing; (3) using other polygons
#' in the geography (e.g., buffers, convex hulls, or concave hulls); (4) finally,
#' a definition of the M can be done via simulations of dispersal based on the
#' species ecological characteristics. For the latter see (****).
#'
#' See functions \code{\link{buffer_area}}, \code{\link{convex_area}},
#' concave_area, and polygon_selection, for options to create accessible areas
#' based on geographic polygons and selection of relevant parts of existing
#' SpatialPolygon areas.
#'
#' For creation of accessible areas from simulations see function
#' \code{\link[grinnell]{M_simulation}} from the grinnell package.

calibration_guide <- function(data, species, longitude, latitude,
                              raster_layers, accessible_area,
                              output_directory = "ellipsenm_calibration") {
  # -----------
  # detecting potential errors
  if (missing(data)) {
    stop("Argument occurrences is necessary to perform the analysis.")
  }
  if (missing(species)) {
    stop("Argument species is not defined.")
  }
  if (missing(longitude)) {
    stop("Argument longitude is not defined.")
  }
  if (missing(latitude)) {
    stop("Argument latitude is not defined.")
  }
  if (missing(raster_layers)) {
    stop("Argument raster_layers is necessary to perform the analysis.")
  }
  if (missing(accessible_area)) {
    stop("Argument accessible_area is necessary to perform the analysis.")
  }

  # -----------
  # preparing directory and saving data
  dir.create(output_directory)

  save(data, species, longitude, latitude, raster_layers, accessible_area,
       file = paste0(output_directory, "/calibration_data.RData"))

  # -----------
  # producing guide
  report_format(col = "#1e84b6", name = paste0(output_directory, "guide_format"))
  guide(guide_type = "calibration", output_directory)
}
