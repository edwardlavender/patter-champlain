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

# Truncated logistic detection probability function
trunclogis <- function(receiver_alpha, receiver_beta, receiver_gamma, dist) {
  p <- plogis(receiver_alpha + receiver_beta * dist)
  p[dist > receiver_gamma] <- 0
  p
}

# Compute the weighted median, matching ggplot behaviour 
# * This function handles empty groupings 
weighted.median <- function(x, w) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  } 
  as.numeric(quantreg::rq(x ~ 1, tau = 0.5, weights = w)$coefficients)
}
