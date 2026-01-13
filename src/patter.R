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

# Define line of sight 
in_line_of_sight <- function(map, x0, y0, x1, y1) {
  xy <- data.table(x0 = x0, y0 = y0, x1 = x1, y1 = y1)
  xy[, xm := (x0 + x1) * 0.5]
  xy[, ym := (y0 + y1) * 0.5]
  !is.na(terra::extract(map, cbind(xy$xm, xy$ym))[, 1])
}

# Sample initial locations for the forward filter
# * This is an R implementation of the Julia code in initalise-filters.jl
pf_filter_fwd_xinit <- function(iter, map, timeline, acoustics, model_move) {

  detections  <- acoustics[obs == 1L, ]
  detections0 <- detections[1, ]
  
  if (min(timeline) == detections0$timestamp) {
    out <- list(xinit = NULL)
    
  } else {
    
    out <- list()
    
    # Define timeline0 from timeline[1] to time of first detection
    timeline0 = seq(min(timeline), min(detections0$timestamp), by = "2 mins")
    
    # Define map0 for first detection
    container0 <- 
      cbind(detections0$receiver_x, detections0$receiver_y) |>
      terra::vect(crs = terra::crs(map)) |> 
      terra::buffer(width = detections0$receiver_gamma, quadsegs = 1000L)
    env0 <- terra::mask(map, container0)
    
    # Define initial states
    xinit0 <- 
      env0 |> 
      terra::spatSample(size = iter$n_particle_filter, xy = TRUE, na.rm = TRUE, replace = TRUE) |> 
      setDT()
    
    # Account for line of sight 
    los    <- in_line_of_sight(map, detections0$receiver_x, detections0$receiver_y, xinit0$x, xinit0$y)
    xinit0 <- xinit0[los == TRUE, ]
    if (nrow(xinit0) != iter$n_particle_filter) {
      xinit0 <- xinit0[sample.int(.N, size = iter$n_particle_filter, replace = TRUE), ]
    }
    
    # Sample headings
    xinit0 <- 
      xinit0 |> 
      lazy_dt() |> 
      mutate(map_value = 1.0, heading = runif(n(), 0, 2 * pi)) |> 
      select("map_value", "x", "y", "heading") |> 
      as.data.table()
   
    # Define dataset (non detections)
    nondetections <- acoustics[timestamp <= detections0$timestamp, ]
    nondetections[, obs := 0L]
    yobs0 <- list(ModelObsAcousticLogisTruncLos = nondetections)
    
    # Run filter backwards from first detection to start of timeline, assuming only non-detections
    bwd0 <- pf_filter(.timeline   = timeline0,
                      .state      = state_trout(),
                      .xinit      = xinit0,
                      .yobs       = yobs0,
                      .model_move = model_move,
                      .n_move     = iter$n_move,
                      .n_particle = iter$n_particle_filter,
                      .n_resample = iter$n_resample,
                      .n_record   = iter$n_particle_filter,
                      .direction  = "backward")
    out$pf_particles <- bwd0
    out$xinit        <- bwd0$states[timestep == 1L, .(map_value, x, y, heading)]
    
  }
  
  out
  
}
