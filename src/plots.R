# Plot a dbn
plot_dbn <- function(dbn, lower = 0, upper = Inf, xlim, add = FALSE, pars = list(), ...) {
  # Define x values
  x <- seq(xlim[1], xlim[2], length.out = 100)
  x <- x[x >= lower & x <= (upper + diff(x)[1])]
  # Compute (truncated) densities 
  pars$x    <- x
  pars$spec <- dbn
  pars$a    <- lower
  pars$b    <- upper
  y <- do.call(truncdist::dtrunc, pars)
  # Plot
  if (!add) {
    plot(x, y, type = "l", ...)
  } else {
    lines(x, y, ...)
  }
  invisible(NULL)
}

# Fit and plot a distribution
plot_fitted_dbn <- function(x, dbn, xlim, add = FALSE, ...) {
  pars <- fitdistrplus::fitdist(x, dbn)
  args <- list(pars$estimate[1], pars$estimate[2])
  names(args) <- names(pars$estimate)
  plot_dbn(dbn = dbn, xlim = xlim, add = add, pars = args, ... )
}

# Add polygon e.g., for density() output
add_poly <- function(dens, col = "dimgrey") {
  polygon(
    x = c(dens$x, rev(dens$x)),
    y = c(dens$y, rep(0, length(dens$y))),
    col = col,
    border = scales::alpha(col, alpha = 1), 
    lwd = 0.25
  )
}

# Add vertical arrow marking mobility
mark_mobility <- function(mobility, col = "black", ...) {
  arrows(x0 = mobility, y0 = -0.003, x1 = mobility, y1 = 0,
         length = 0.04, lwd = 1.25, col = col, ...)
}
