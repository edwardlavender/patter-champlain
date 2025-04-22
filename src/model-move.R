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
if (patter:::julia_session()) {
  model_move_xinit <- function(.map, .xinit) {
    # Copy data.table (x, y coordinates)
    xinit <- copy(.xinit)
    # Add map value 
    x <- y <- heading <- NULL
    xinit[, map_value := terra::extract(.map, cbind(x, y))[, 1]]
    stopifnot(all(!is.na(xinit$map_value)))
    # Simulate heading 
    xinit[, heading := runif(.N, 0, 2 * pi)]
    # Organise
    xinit[, .(map_value, x, y, heading)]
  }
}
