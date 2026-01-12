# model_obs_acoustic_logis_trunc wrapper for model_obs_acoustic_logis_trunc_los
model_obs_acoustic_logis_trunc_los <- function(.data, .strict = TRUE) {
  .data <- copy(.data)
  # Set sensor_id column
  if (rlang::has_name(.data, "receiver_id")) {
    setnames(.data, "receiver_id", "sensor_id")
  }
  # Select columns
  cols <- c("sensor_id",
            "receiver_x", "receiver_y",
            "receiver_alpha", "receiver_beta", "receiver_gamma")
  check_names(.data, cols)
  if (.strict) {
    .data <-
      .data |>
      select(all_of(cols)) |>
      as.data.table()
  }
  # Define structure
  .data        <- list(.data)
  names(.data) <- "ModelObsAcousticLogisTruncLos"
  structure(
    .data,
    class = c("list", "ModelObs", "ModelObsAcousticLogisTrunc", "ModelObsAcousticLogisTruncLos")
  )
  
}
