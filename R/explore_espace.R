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
#' @param save (logical) whether or not to save the plot as pdf file.
#' @param name (character) name (including the extention .pdf) of the PDF file to
#' be created with the exploratory figures. Default = "e_space_exploration.pdf".
#' @param height (numeric) height of the figure; default = 8.
#' @param width (numeric) width of the figure; default = 17.
#' @param open (logical) whether or not to open the file created directly.
#' Default = TRUE
#'
#' @export
#'
#' @return
#' A plot and a PDF document with figures to explore the data in environmental
#' space, which could help to make desicions on which methods and variables are
#' worth to explore during model calibration.
#'
#' @usage
#' explore_espace(data, species, longitude, latitude, raster_layers,
#'                save = FALSE, name = "e_space_exploration.pdf",
#'                height = 8, width = 17, open = TRUE)

explore_espace <- function(data, species, longitude, latitude, raster_layers,
                           save = FALSE, name = "e_space_exploration.pdf",
                           height = 8, width = 17, open = TRUE) {
  if (missing(data)) {
    stop("Argument 'data' is necessary to perform the analysis.")
  }
  if (missing(species)) {
    stop("Argument 'species' is not defined.")
  }
  if (missing(longitude)) {
    stop("Argument 'longitude' is not defined.")
  }
  if (missing(latitude)) {
    stop("Argument 'latitude' is not defined.")
  }
  if (missing(raster_layers)) {
    stop("Argument 'raster_layers' is not defined.")
  }

  suppressPackageStartupMessages(library(GGally))
  suppressPackageStartupMessages(library(ggplot2))

  # -----------
  # data preparation
  sp_name <- as.character(data[1, species])

  ## values of variables
  v <- data.frame(code = "M", na.omit(raster::values(raster_layers)))
  if (nrow(v) > 10000) {v <- v[sample(nrow(v), 10000), ]}

  o <- data.frame(code = "Species",
                  raster::extract(raster_layers, data[, c(longitude, latitude)]))
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
      stat_density_2d(data = data[subset[[2]], ], colour = "black", lwd = 0.3,
                      show.legend = FALSE, bins = 10)
  }

  # -----------
  # plot
  ## individual pair plots
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

  ## combined pair plots
  p1 <- cowplot::plot_grid(ggmatrix_gtable(gp1), ggmatrix_gtable(gp2))

  ## plotting and exporting if needed
  print(p1)

  if (save == TRUE) {
    pdf(name, height = height, width = width)
    print(p1)
    invisible(dev.off())

    ### opening the file if needed
    if (open == TRUE) {
      path <- paste0(getwd(), "/", name)
      invisible(system(paste0('open "', path, '"')))
    }

    cat("\nA PDF file named", name, "has been created, check your working directory:\n",
        getwd(), "\n")
  }

  interpretation_espace(sp_name)
}


#' Helper function for interpretation of exploratory figures
#' @param sp_name (character) name od the species of interest. Default = "species".
#' @export
#' @return
#' A text with help for selecting variables and methods based on exploratory plots.

interpretation_espace <- function(sp_name = "species") {
  cat("\n---------------------------------------------------------------------------------------\n")
  cat("**************************ellipsenm: model calibration process*************************\n\n")
  cat("Arrangement of occurrence records in environmental space\n")
  cat("---------------------------------------------------------------------------------------\n")
  cat("\nSpecies:", sp_name, "\n\n")
  cat("Description\n")
  cat(" Plots produced with explore_espace help to visualize the arrangement of species\n",
      "occurrences in the accessible environment. With these visualizations, two main\n",
      "aspects can be recognized: (1) how available environmental conditions are distributed\n",
      "and clustered in the environmental space (given by the kernel in blue); and (2) how\n",
      "the species occurrences are distributed and clustered in reference to the available\n",
      "environmental space.\n",
      "Based on the distinct ways in which the data can be arranged considering different\n",
      "pairs of variables, decisions about what variables are more appropriate to use and\n",
      "the method to be used to create the ellipsoid that encloses the points can be taken.\n")
  cat("\nExample cases and recommendations\n")
  cat("  Selection of variables\n")
  cat("    1) A variable in which the space used by the species is narrow compared to the\n",
      "      space available may be preferable for further analyses because this fact could\n",
      "      indicate that the species tolerance limits for that variable are clear. This is,\n",
      "      if a set of variable values is not being used by the species even if they are\n",
      "      available to it, this indicates that the species may not tolerate such conditions.\n")
  cat("\n  Selection of method\n")
  cat("    Currently three methods are ready to use in ellipsenm:\n",
      "   covmat = normal ellipsoids based on a centroid and the covariance matrix.\n",
      "   mve1   = minimum volume ellipsoids (as in MASS::cov.rob() using the mve method), in\n",
      "           which the modifications to the centroid and covariance matrix are made so\n",
      "           the resulting enclosing ellipsoid has the smallest volume.\n",
      "   mve2   = moment based minimum volume ellipsoid, as with the previous method, this\n",
      "           one also creates a minimum volume ellipsoid. However, the algorithm used\n",
      "           is different. The ellipsoids created may suffer changes in the position of\n",
      "           the centroid and the values of the covariance matrix if compared with the\n",
      "           \"covmat\" method. In general,ellipsoids created with this method have smaller\n",
      "           volumes than the ones created with previous methods.\n\n")
  cat("    Selecting the method could be directed by the patterns detected in exploratory plots\n",
      "   as indicated below:\n")
  cat("    2) If occurrence data is surrounded entirely by the available conditions and areas\n",
      "      of highest density of records are close to the center of the points' clould, any\n",
      "      of the three methods could be useful. Despite the differences that they present,\n",
      "      basic assumptions related to niche centroid and niche limits may not be affected.\n")
  cat("    3) If occurrence data is located marginally and areas of highest density of records\n",
      "      are close to the center of the clould of points, the three methods may also be\n",
      "      explored during calibration.\n")
  cat("    4) If occurrence data is located marginally and areas of highest density of records\n",
      "      are close to the border of the available contiions, the limits of the species niche\n",
      "      are unknown towards such margins. In this case explorein the covmat method may be\n",
      "      better, the mve1 method may also be useful however no centainty exist.\n")
  cat("\n  The previous recommendations are only examples, other arguments to select methods and\n",
      " variables can be defined by the researcher considering the goals of the study.\n")
  cat("\n  Other recommendations will be added as other methods are made available in the future.")
  cat("\n---------------------------------------------------------------------------------------\n")
}

