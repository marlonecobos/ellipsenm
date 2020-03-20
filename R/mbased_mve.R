#' Fit moment-based minimum volumen ellipsoids
#'
#' @description mbased_mve helps in finding the centroid and matrix that
#' define a minimum volume ellipsoid as proposed by Moshtagh (2005).
#'
#' @param data matrix of values to be used for fitting the ellipsoid. Columns
#' represent dimensions and rows observations.
#' @param fitting_tolerance (numeric) proportion of error allowed when checking
#' if the ellipsoid incloses all values considered to fit it. Default = 0.001.
#'
#' @return
#' A named list containing the values for the centroid of th ellipsoid and a
#' matrix of covariance.
#'
#' @details
#' The algorithm used to fit the ellipsoids here was implemented by
#' Moshtagh (2005) and based on initial work of Khachiyan (1996) and
#' Rocha et al. (2002).
#'
#' Details about to the algortihm can be found in \url{https://bit.ly/2XYWlVT},
#' \url{https://bit.ly/2x6aR2s}, and \url{https://bit.ly/2KrxE1g}.
#'
#' @export
#'
#' @examples
#' # reading data
#' occurrences <- read.csv(system.file("extdata", "occurrences.csv",
#'                                     package = "ellipsenm"))
#'
#' # raster layers of environmental data
#' vars <- raster::stack(list.files(system.file("extdata", package = "ellipsenm"),
#'                                  pattern = "bio", full.names = TRUE))
#'
#' data <- raster::extract(vars, occurrences[, -1])
#'
#' mvellipsoid <- mbased_mve(data)

mbased_mve <- function(data, fitting_tolerance = 0.001) {
  # -----------
  # detecting potential errors
  if (missing(data)) {
    stop("Argument 'data' is necessary to perform the analysis")
  }

  # -----------
  # preparing data
  var_names <- colnames(data)
  data <- t(data)
  d <- nrow(data)
  n <- ncol(data)

  q <- rbind(data, rep(1, n))

  err <- 1
  u <- rep(1 / n, n)

  # -----------
  # fitting the ellipsoid according to 1 and 2 moments of data
  while (err > fitting_tolerance) {
    prod_x <- q %*% diag(u, n, n) %*% t(q)
    prod_m <- diag(t(q) %*% solve(prod_x) %*% q)
    max_m <- max(prod_m)

    step_size <- (max_m - d - 1) / ((d + 1) * (max_m - 1))
    new_u <- (1 - step_size) * u
    j <- order(prod_m)[n]
    new_u[j] <- new_u[j] + step_size

    ll <- svd(new_u - u)
    err <- max(ll[[1]])
    u <- new_u
  }

  # -----------
  # preparing results
  u1 <- diag(u, n, n)
  du <- data %*% u
  covari <- solve(solve((data %*% u1 %*% t(data)) - (du %*% t(du))) / d)

  cent <- c(du)
  names(cent) <- var_names

  return(list(centroid = cent, covariance_matrix = covari))
}
