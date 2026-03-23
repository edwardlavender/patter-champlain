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

# Quick plots of simulated occupancy distributions & simulated paths 
if (!patter:::os_linux() | (patter:::os_linux() & !patter:::julia_session())) {
  
  lapply_qplot_sim <- function(.iteration, .n_plot = 4L) {
    ind <- sample.int(nrow(.iteration), size = .n_plot)
    pp <- par(mfrow = prettyGraphics::par_mf(.n_plot))
    lapply(ind, function(i) {
      it   <- .iteration[i, ]
      map  <- terra::rast(it$file_occupancy_sim)
      path <- qs::qread(it$file_path_sim)
      # Plot occupancy distribution & add simulated path
      terra::plot(map, main = it$unit_id[i])
      patter:::add_sp_path(path$x, path$y, length = 0.01, lwd = 0.05) |> 
        suppressWarnings()
    })
    par(pp)
    nothing()
  }
  
}