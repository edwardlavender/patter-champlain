# Get parameters of shrunk/widened gamma distribution
gamma_rescale <- function(shape, scale, fact = 1.2) {
  mode      <- (shape - 1) * scale
  new_shape <- shape / fact
  new_scale <- mode / (new_shape - 1)
  c(new_shape, new_scale)
}
