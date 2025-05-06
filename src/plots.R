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


if (!patter:::os_linux() | (patter:::os_linux() & !patter:::julia_session())) {
  
  # Animate particles from the particle filter
  # * Code adapted for debugging from https://github.com/edwardlavender/patter-flapper/blob/main/src/debug-convergence.R
  # @param .sim The data.table row
  # @param .map The SpatRaster
  # @param .moorings The moorings data.table
  # @param .start,.stop Integers that define the time steps of interest
  # @param .input The named list of arguments passed to pf_filter()
  # @param .output The named list of outputs from pf_filter()
  
  ani <- function(.sim, 
                  .map, .moorings, 
                  .start = 1L, .end, 
                  .input, .output, 
                  .tnow,
                  .cl = 1L) {
    
    .moorings <- copy(.moorings)
    
    # Create directories
    frames <- here_fig("debug", .tnow, "frames")
    dir.create(frames, recursive = TRUE)
    mp4 <- here_fig("debug", .tnow)
    dir.create(mp4)
    
    # Make frames
    pf_plot_xy(.map = .map, 
               .coord = .output$states, 
               .steps = .start:.end,
               .png = list(filename = frames),
               .cl = .cl,
               .add_points = list(pch = ".", col = "red"),
               .add_layer = function(t) {
                 
                 # Add receivers
                 text(.moorings$receiver_x, .moorings$receiver_y, .moorings$receiver_id, cex = 0.5)
                 
                 # Add detection containers
                 # cbind(.moorings$receiver_x, .moorings$receiver_y) |>
                 #   terra::vect() |>
                 #   terra::buffer(width = .moorings$receiver_gamma) |> 
                 #   terra::lines()
                 
                 # Add acoustic containers
                 containers <- .input$.yobs$ModelObsContainer
                 if (!is.null(containers)) {
                   cinfo <- containers[timestamp == .input$.timeline[t], ]
                   if (nrow(cinfo) > 0L) {
                     
                     # Colour receiver(s) with next detection
                     text(cinfo$centroid_x, cinfo$centroid_y, cinfo$sensor_id, cex = 0.5, col = "blue", font = 2)
                     
                     # Add time-specific acoustic containers
                     cbind(cinfo$centroid_x, cinfo$centroid_y) |>
                       terra::vect() |>
                       terra::buffer(width = cinfo$radius) |> 
                       terra::lines(col = "royalblue") 
                   }
                 }
               })
    
    # Make animation 
    input   <- gtools::mixedsort(list.files(frames, full.names = TRUE))
    output  <- file.path(mp4, "ani.mp4")
    av::av_encode_video(input, output, framerate = 5)
    
    # Open animation (on MacOS)
    system(paste("open", shQuote(output)))
    invisible(NULL)
  }
  
}

