# Get parameters of shrunk/widened gamma distribution
gamma_rescale <- function(shape, scale, fact = 1.2) {
  mode      <- (shape - 1) * scale
  new_shape <- shape / fact
  new_scale <- mode / (new_shape - 1)
  c(new_shape, new_scale)
}

# Mixture distribution for turning angle
# MixtureModel([truncated(Normal(0.0, 0.4), -pi, pi), Uniform(-pi, pi)], [0.99, 0.01])
dmix <- function(x, mean = 0, sd = 0.3) {
  # dnorm(x, mean = mean, sd = sd)
  d1 <- truncdist::dtrunc(x, spec = "norm", a = -pi, b =  pi, mean = mean, sd = sd)
  d2 <- dunif(x, -pi, pi)
  # w1 = 0.25 empiricaly matches the VPS data
  w1 <- 0.99
  w2 <- 1 - w1
  d1 * w1 + d2 * w2
}
