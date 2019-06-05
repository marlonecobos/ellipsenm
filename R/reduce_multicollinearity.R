reduce_multicollinearity  <- function(data, longitude, latitude,
                                      raster_layers, sample_size = 10000,
                                      correlation_limit = 0.8, VIF_limit = 10) {
  # -----------
  # detecting potential errors a
  if (missing(raster_layers)) {
    stop("Argument raster_layers is not defined.")
  }

  if (class(raster_layers)[1] == "RasterStack") {
    if (dim(raster_layers)[3] < 2) {
      stop("raster_layers must contain at least 2 layers.")
    }
  } else {
    stop("raster_layers must be a RasterStack.")
  }

  # -----------
  # raster processing and sampling raster data
  raster_values <- na.omit(raster::values(raster_layers))

  if (nrow(raster_values) > sample_size) {
    raster_values <- raster_values[sample(nrow(raster_values), sample_size), ]
  }

  # -----------
  # correlation groups


  # -----------
  # VIF groups


  # -----------
  # return object
  groups <- list(
    corelation_results = n_cor_groups,
    VIF_results = list(totally_collinear = al_var,
                       selected_variables = sel_vars)
    )


  return()
}

#' Helper function to find non-correlated groups of variables
#'

noncor_groups <- function(raster_values, correlation_limit) {
  correlations <- cor(raster_values)
  rule <- correlations > correlation_limit | correlations < -correlation_limit
  variables <- colnames(rule)

  cor_n <- sort(sapply(variables, function(x) {sum(rule[, x]) - 1}))

  varcor <- names(cor_n[cor_n != 0])
  cors <- lapply(varcor, function(x) {
    corel <- rule[, x]
    corvars <- names(corel[corel == TRUE & variables != x])
  })

  names(cors) <- varcor
  corl <- list(varcor, cors)

  lvar <- list()

  for (x in 1:length(cors)) {
    vec_cor <- cors[[x]]
    nvec <- length(vec_cor)

    svar <- vector()
    for (y in 1:length(vec_cor)) {
      ncor <- cors[[vec_cor[y]]]

      if (length(ncor) > nvec){
        svar[y] <- names(cors)[x]
        cors <- cors[-which(names(cors) == vec_cor[y])]
      } else {
        svar[y] <- "other"
      }
    }

    lvar[[x]] <- unique(svar)
    names(lvar[[x]]) <- vec_cor[[y]]

    if (x == length(cors)) break()
  }

  no_cor <- names(cor_n[cor_n == 0])

  no_cor_groups <- lapply(no_cor_sets, function (x) {c(no_cor, x)})
}


#' Helper function to find non-collinear groups of variables based on VIF
#'

vif_groups <- function(data, longitude, latitude, raster_layers,
                       raster_values = NULL) {
  d_r <- na.omit(raster::extract(raster_layers, data[, c(longitude, latitude)]))

  if (is.null(raster_values)) {
    raster_values <- na.omit(raster::values(raster_layers))
  }

  d_f <- as.data.frame(rbind(cbind(presence = 1, d_r),
                             cbind(presence = 0, raster_values)))

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
}
