thin_data <- function(data, longitude, latitude, thin_distance) {
  data <- data[!is.na(data[, longitude]) & !is.na(data[, latitude]), ]
  dat_sp <- sp::SpatialPointsDataFrame(data[, c(longitude, latitude)], data)
  dat_sp1 <- sp::remove.duplicates(dat_sp, zero = thin_distance)

  return(dat_sp1@data)
}
