explore_espace <- function(data, species, longitude, latitude, raster_layers,
                           color_palette = viridis::viridis, save = FALSE,
                           name = "e_space_exploration.pdf") {

  # -----------
  # data preparation
  sp_name <- as.character(data[1, species])

  ## all combinations of two variables
  vnames <- names(raster_layers)
  cb_vars <- combn(vnames, m = 2, simplify = FALSE)

  ## values of variables
  raster_layers <- na.omit(raster::values(raster_layers))

  # -----------
  # plot
  if (length(cb_vars) > 2) {
    init_length <- 3
    sec_length <- length(cb_vars)

    ### first page
    layout(mat = cbind(c(1, 2), c(1, 3)))
    plots <- lapply(1:init_length, function(x){
      if (x == 1) {
        par(cex = 0.6, mar = (c(0, 0.5, 3.5, 1.2) + 0.1))
        plot.new()
        title(main = paste0("Arrangement of ", sp_name," occurrences in environmental space"), cex = 1.2)
        legend("topleft", legend = "Description", cex = 1.1, bty = "n")
        legend("topleft",
               legend = c("", "",
                          "  The plots presented below help to visualize the arrangement of the species occurrences in the accessible environment.",
                          "  With these visualizations two main aspects can be recognized: (1) how the environment is distributed and clustered in",
                          "  the environmental space (given by the kernel in light blue); and (2) how the occurrences are distributed and clustered",
                          "  in reference to the available environmental space.", "", "",
                          "  Based on the distinct ways in which the data can be arranged considering distinct paris of variables, decisions about",
                          "  what variables are more appropriate to use and the method to be used to create the ellipsoid that encloses the points.",
                          "  For instance, if occurrence data is sorrouded "),
               cex = 0.9, bty = "n")

        legend("topright", legend = "Symbology        ", cex = 1.1, bty = "n")
        legend("topright", legend = sym_legend, lty = lin, pch = poi, col = colss, cex = 0.9, bty = "n")

      } else {
        par(cex = 0.6, mar = (c(4.5, 4, 3.5, 1) + 0.1))
        ### actual values
        hist(as.numeric(names(env_data[[x - 1]])), main = gsub("_", " ", sp_names[x - 1]), xlab = "Variable values")
        points(as.numeric(names(occ_data[[x - 1]])), rep(y_values[[x - 1]][1], length(occ_data[[x - 1]])), col = "darkgreen", pch = 1)
        abline(v = limits[x - 1, ], col = col1)

        ### median deviation
        hist(env_data[[x - 1]], main = gsub("_", " ", sp_names[x - 1]), xlab = "Median deviation")
        points(occ_data[[x - 1]], rep(y_values[[x - 1]][2], length(occ_data[[x - 1]])), col = "darkgreen", pch = 1)
        limit <- floor(CL_lines * length(env_data[[x - 1]]) / 100)
        abline(v = sort(env_data[[x - 1]])[limit], col = col)
      }
    })

    ### next pages
    par(mfrow = c(5, 2), cex = 0.6, mar = (c(4.5, 4, 3.5, 1) + 0.1))
    plots <- lapply(5:sec_length, function(x){
      ### actual values
      hist(as.numeric(names(env_data[[x]])), main = gsub("_", " ", sp_names[x]), xlab = "Variable values")
      points(as.numeric(names(occ_data[[x]])), rep(y_values[[x]][1], length(occ_data[[x]])), col = "darkgreen", pch = 1)
      abline(v = limits[x, ], col = col1)

      ### median deviation
      hist(env_data[[x]], main = gsub("_", " ", sp_names[x]), xlab = "Median deviation")
      points(occ_data[[x]], rep(y_values[[x]][2], length(occ_data[[x]])), col = "darkgreen", pch = 1)
      limit <- floor(CL_lines * length(env_data[[x]]) / 100)
      abline(v = sort(env_data[[x]])[limit], col = col)
    })

  } else {
    init_length <- length(cb_vars) + 1
    ## columns of layout
    c1 <- c(1, seq(2, (length(cb_vars) * 2), 2))
    c2 <- seq(1, ((length(cb_vars) * 2) + 1), 2)

    ## first page
    layout(mat = cbind(c1, c2))
    plots <- lapply(1:init_length, function(x){
      if (x == 1) {
        par(cex = 0.6, mar = (c(0, 0.5, 3.5, 1.2) + 0.1))
        plot.new()
        title(main = paste0("Values and median deviations for variable ", variable_name))
        legend("topleft", legend = "Description", cex = 1.1, bty = "n")
        legend("topleft", legend = c("", "", "  The information presented below helps to visualize the",
                                     "  distribution of values of the variable in the accessible",
                                     "  area as well as the species occurrences to facilitate the ",
                                     "  delimitation of conditions to be used in further analyses."), cex = 0.9, bty = "n")

        legend("topright", legend = "Symbology        ", cex = 1.1, bty = "n")
        legend("topright", legend = sym_legend, lty = lin, pch = poi, col = colss, cex = 0.9, bty = "n")

      } else {
        par(cex = 0.6, mar = (c(4.5, 4, 3.5, 1) + 0.1))
        ### actual values
        hist(as.numeric(names(env_data[[x - 1]])), main = gsub("_", " ", sp_names[x - 1]), xlab = "Variable values")
        points(as.numeric(names(occ_data[[x - 1]])), rep(y_values[[x - 1]][1], length(occ_data[[x - 1]])), col = "darkgreen", pch = 1)
        abline(v = limits[x - 1, ], col = col1)

        ### median deviation
        hist(env_data[[x - 1]], main = gsub("_", " ", sp_names[x - 1]), xlab = "Median deviation")
        points(occ_data[[x - 1]], rep(y_values[[x - 1]][2], length(occ_data[[x - 1]])), col = "darkgreen", pch = 1)
        limit <- floor(CL_lines * length(env_data[[x - 1]]) / 100)
        abline(v = sort(env_data[[x - 1]])[limit], col = col)
      }
    })
  }
  invisible(dev.off())
}


library(GGally)
library(ggplot2)

vo <- rbind(v, o)

my_theme <- theme(panel.background = element_rect(fill = "white", colour = "black"))

f_1 <- function(data, mapping, subset = NULL) {
  ggplot(mapping = mapping) +
    geom_point(data = data[subset[[1]], ], colour = "lightblue1", alpha = 0.3) +
    geom_point(data = data[subset[[2]], ], colour = "grey25", shape = 1) +
    my_theme
}

f_2 <- function(data, mapping, subset = NULL) {
  ggplot(mapping = mapping) +
    stat_density_2d(data = data[subset[[1]], ], geom = "polygon",
                    show.legend = FALSE, bins = 350) +
    scale_fill_continuous(low = "lightblue1", high = "dodgerblue4") +
    stat_density_2d(data = data[subset[[2]], ], colour = "black",
                    show.legend = FALSE, bins = 10) +
    my_theme
}

cors <- function(data, mapping, subset = NULL){
  data <- data[subset[[2]], ]
  x <- as.numeric(data[, 1])
  y <- as.numeric(data[, 2])
  est <- cor(x, y)
  lb_size <- 5 * abs(est)
  lbl <- paste0("r = ", round(est, 3))

  ggplot(data = data, mapping = mapping) +
    annotate("text", x = mean(x), y = mean(y), label = lbl, size = lb_size) +
    my_theme
}


subsets <- list(as.character(vo$code) == "M", as.character(vo$code) == "Species")
ggpairs(vo, columns = 2:ncol(vo),
        upper = list(continuous = wrap(cors, subset = subsets)),
        lower = list(continuous = wrap(f_1, subset = subsets))) +
  my_theme


ggpairs(vo, columns = 2:ncol(vo),
        upper = list(continuous = wrap(cors, subset = subsets)),
        lower = list(continuous = wrap(f_2, subset = subsets))) +
  my_theme


