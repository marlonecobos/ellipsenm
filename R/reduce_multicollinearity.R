reduce_multicollinearity  <- function(data, raster_layers, sample_size = 10000,
                                      correlation_limit = 0.8, VIF_limit = 10) {
  # -----------
  # detecting potential errors a
  if (missing(raster_layers)) {
    stop("Argument raster_layers is not defined.")
  }

  if (class(raster_layers)[1] == "RasterStack") {
    if (dim(raster_layers)[3]) {
      stop("raster_layers must contain at least 2 layers.")
    }
  } else {
    stop("raster_layers must be a RasterStack.")
  }

  # -----------
  # raster processing
  rvalues <- na.omit(raster::values(raster_layers))

  if (nrow(rvalues) > sample_size) {
    rvalues <- rvalues[sample(nrow(rvalues), sample_size), ] # to get the sample
  }

  # -----------
  # correlation matrix and groups
  correlations <- cor(rvalues)
  rule <- correlations > correlation_limit | correlations < -correlation_limit
  correlations <-
  correlations <-
  n_cor_groups

  # -----------
  # VIF groups
  d_r <- na.omit(raster::extract(raster_layers, data[, c(longitude, latitude)]))

  d_f <- as.data.frame(rbind(cbind(presence = 1, d_r),
                             cbind(presence = 0, rvalues)))

  form <- paste0("presence ~ ", paste0("I(", colnames(d_f)[-1], ")"))

  ## glms
  glms <- glm(form, data = d_f, family = "binomial") ########

  ## removal of variables if needed
  al_var <- attributes(alias(glms)$Complete)$dimnames #########
  al_vars <- al_var[[1]]

  if (length(al_vars) > 0) {
    d_f <- d_f[, !colnames(d_f) %in% al_vars]
    form <- paste0("presence ~ ", paste0("I(", colnames(d_f)[-1], ")"))

    glms <- glm(form, data = d_f, family = "binomial")
  }

  ## vif claculations
  vifs <- car::vif(glms)

  ## selected variables
  sel_vars <- vifs[vifs < VIF_limit]

  # -----------
  # return object
  groups <- list(
    corelation_results = n_cor_groups,
    VIF_results = list(totally_collinear = al_var,
                       selected_variables = sel_vars)
    )


  return()
}
