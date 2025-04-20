# Add truncated gamma curves to ggplot2
stat_dtruncgamma <- function(mobility, shape, scale, col = "black", ...) {
  ggplot2::stat_function(fun = function(x) {
    truncdist::dtrunc(x, "gamma", a = 0, b = mobility, shape = shape, scale = scale)
  }, colour = col, ...)
}