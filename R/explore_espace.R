#' Exploration of data in enviornmental space
#'
#' @description explore_espace helps in eploring the patterns of arrangement of
#' species occurrence data and accessible areas in environmental space.
#'
#' @param data data.frame of occurrence records. Columns must be: species,
#' longitude, and latitude.
#' @param species (character) name of the column with the name of the species.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @param raster_layers RasterStack of at least two environmental variables to be
#' extracted using geographic coordinates present in \code{data}.
#' @param name (character) name (including the extention .pdf) of the PDF file to
#' be created with the exploratory figures. Default = "e_space_exploration.pdf".
#' @param open (logical) whether or not to open the file created directly.
#' Default = TRUE
#'
#' @export
#'
#' @return
#' A PDF document with indications and figures to explore the data in environmental
#' space and make desicions on which methods and variables are worth to explore
#' during model calibration.

explore_espace <- function(data, species, longitude, latitude, raster_layers,
                           name = "e_space_exploration.pdf", open = TRUE) {
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
    stop("Argument raster_layers is not defined.")
  }

  suppressPackageStartupMessages(library(GGally))
  suppressPackageStartupMessages(library(ggplot2))

  # -----------
  # data preparation
  sp_name <- as.character(data[1, species])

  ## values of variables
  v <- data.frame(code = "M", na.omit(raster::values(raster_layers)))
  o <- data.frame(code = "Species",
                  raster::extract(vars, data[, c(longitude, latitude)]))
  vo <- rbind(v, o)

  ## subsets
  subsets <- list(vo$code == "M", vo$code == "Species")

  # -----------
  # theme and functions
  my_theme <- theme(panel.background = element_rect(fill = "white", colour = "black"),
                    panel.grid = element_blank())

  f_1 <- function(data, mapping, subset = NULL) {
    ggplot(mapping = mapping) +
      geom_point(data = data[subset[[1]], ], colour = "lightblue1", alpha = 0.3) +
      geom_point(data = data[subset[[2]], ], colour = "grey25", shape = 16,
                 alpha = 0.6)
  }

  f_2 <- function(data, mapping, subset = NULL) {
    ggplot(mapping = mapping) +
      stat_density_2d(data = data[subset[[1]], ], aes(fill = (..level..)),
                      geom = "polygon", show.legend = FALSE, bins = 350) +
      scale_fill_continuous(low = "lightblue1", high = "dodgerblue4") +
      stat_density_2d(data = data[subset[[2]], ], colour = "black",
                      show.legend = FALSE, bins = 10)
  }

  # -----------
  # plot
  gp1 <- ggpairs(vo[, 2:ncol(vo)],
                 upper = list(continuous = "cor"),
                 lower = list(continuous = wrap(f_1, subset = subsets)),
                 progress = T) +
    my_theme

  gp2 <- ggpairs(vo[, 2:ncol(vo)],
                 upper = list(continuous = "cor"),
                 lower = list(continuous = wrap(f_2, subset = subsets)),
                 progress = T) +
    my_theme

  pdf(name, height = 8, width = 17)
  plot.new()
  title(main = paste0("ellipsenm: model calibration process"),
        cex.main = 2)
  legend("topleft", legend = paste0("Arrangement of ", sp_name,
                                    " occurrences records in environmental space"),
         cex = 1.4, bty = "n")
  legend("topleft", legend = "\n\nDescription", cex = 1.2, bty = "n")
  legend("topleft",
         legend = c("\n", "", "", "",
                    " The plots presented below help to visualize the arrangement of the species occurrences in the accessible environment.",
                    " With these visualizations two main aspects can be recognized: (1) how available environmental conditions are distributed and ",
                    " clustered in the environmental space (given by the kernel in blue); and (2) how the occurrences are distributed and clustered",
                    " in reference to the available environmental space.", "",
                    " Based on the distinct ways in which the data can be arranged considering distinct paris of variables, decisions about",
                    " what variables are more appropriate to use and the method to be used to create the ellipsoid that encloses the points.",
                    " For instance, if occurrence data is sorrouded "),
         cex = 1.1, bty = "n")
  cowplot::plot_grid(ggmatrix_gtable(gp1), ggmatrix_gtable(gp2))
  invisible(dev.off())

  # opening the file
  path <- paste0(getwd(), "/", name)
  if (open == TRUE) {
    invisible(system(paste0('open "', path, '"')))
  }

  cat("A PDF file named", name, "has been created, check your working directory:\n",
      getwd())
}








