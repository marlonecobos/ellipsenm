calibration_guide <- function(data, species, longitude, latitude, variables,
                              accessible_area,
                              output_directory = "ellipsenm_calibration") {
  # -----------
  # detecting potential errors


  # -----------
  # preparing data
  save(data, variable_names, n_var, ell_meta, mean_pred, prevalences,
       replicates, replicate_type, bootstrap_percentage, color_palette,
       file = paste0(output_directory, "/calibration_data.RData"))

  # -----------
  # producing guide
  ellipsenm::report_format(col = "#1e84b6", name = "calibration_format")
  report(report_type = "calibration", output_directory)

  file.edit(paste0(output_directory, "/enm_calibration_guide.Rmd"))
}
