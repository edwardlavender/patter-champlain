# Plot a dbn
plot_dbn <- function(dbn, xlim, add = FALSE, pars = list(), ...) {
  # Define x values
  x <- seq(xlim[1], xlim[2], length.out = 100)
  # Compute densities 
  ddbn <- paste0("d", dbn)
  pars$x <- x
  y <- do.call(ddbn, pars)
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