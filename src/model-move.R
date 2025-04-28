state_trout <- function() {
  "StateCXY"
}

# model_move_cxy() wrapper for trout
# * .pars is a named list or 1-row data.table with model parameters
model_move_trout <- function(.pars) {
  proj.build::check_names(.pars, c("mobility", "shape", "scale", "phi"))
  mobility <- .pars$mobility 
  shape    <- .pars$shape
  scale    <- .pars$scale
  phi      <- .pars$phi
  patter::model_move_cxy(.mobility = mobility,
                         .dbn_length = glue::glue("truncated(Gamma({shape}, {scale}), upper = {mobility})"),
                         .dbn_heading_delta = glue::glue("Normal(0.0, {phi})"))
}

# Initialise the movement process
# * xinit is a data.table with map_value, x, y
model_move_xinit <- function(.xinit, .n_particle) {
  # Check inputs
  proj.build::check_names(.xinit, c("map_value", "x", "y"))
  if ((nrow(.xinit) > 1L) & .n_particle > 1L) {
    stop("If .n_particle > 1, nrow(.xinit) should be 1.")
  }
  # Duplicate data.table for .n_particle(s)
  map_value <- x <- y <- NULL
  .xinit <- copy(.xinit[, list(map_value, x, y)])
  out <- lapply(seq_len(.n_particle), function(d) {
    d <- copy(.xinit)
  }) |> rbindlist()
  # Add random variables
  heading <- NULL
  out[, heading := runif(.N, 0, 2 * pi)]
  out[, list(map_value, x, y, heading)]
  out
}
