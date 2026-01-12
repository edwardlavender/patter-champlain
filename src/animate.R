if (!patter:::os_linux() | (patter:::os_linux() & !patter:::julia_session())) {
  
  animate_xyt <- function(.map, .coord, ..., .folder) {
    # Define folder for images
    frames <- file.path(.folder, as.numeric(Sys.time()), "frames")
    dirs.create(frames)
    # Create images
    plot_xyt(.map = .map, 
             .coord = .coord,
             .add_points = list(pch = ".", col = "red"),
             .png = list(filename = frames, 
                         height = 4, width = 2, units = "in", res = 200), ...)
    # Make video
    input   <- gtools::mixedsort(list.files(frames, full.names = TRUE))
    output  <- file.path(dirname(frames), "ani.mp4")
    av::av_encode_video(input, output, framerate = 10)
    
  }

  
  # Animate particles from the particle filter for the AC algorithm
  # * Code adapted for debugging from https://github.com/edwardlavender/patter-flapper/blob/main/src/debug-convergence.R
  # @param .iter The data.table row
  # @param .map The SpatRaster
  # @param .steps Integers that define the time steps of interest
  # @param .input The named list of arguments passed to pf_filter()
  # @param .output The named list of outputs from pf_filter()
  # @param .outdir The output directory 
  
  animate_ac <- function(.iter, 
                         .map, 
                         .steps,
                         .input, .output, 
                         .outdir,
                         .cl = 1L) {

    # Checks
    .steps <- .steps[.steps <= length(.input$.timeline)]
    
    # Create directories
    dir.create(.outdir, recursive = TRUE)
    frames <- file.path(.outdir, "frames")
    unlink(frames, recursive = TRUE)
    dir.create(frames, recursive = TRUE)
    mp4 <- .outdir
    
    # Record parameter settings
    sink(file.path(.outdir, "arguments.txt"))
    print(.input)
    sink()
    
    # Define datasets
    timeline   <- .input$.timeline
    acoustics  <- copy(.input$.yobs$ModelObsAcousticLogisTruncLos)
    acoustics[, timestep := (1:length(timeline))[match(timestamp, timeline)]]
    detections <- acoustics[obs == 1L, ]
    moorings <- 
      acoustics |> 
      lazy_dt() |> 
      group_by(sensor_id) |> 
      slice(1L) |> 
      as.data.table()
    
    # Precompute likelihood surfaces for each receiver
    likelihoods <- lapply(split(moorings, seq_len(nrow(moorings))), function(m) {
      # m <- moorings[1, ]
      # Define receiver position 
      v <- 
        cbind(m$receiver_x, m$receiver_y) |> 
        terra::vect(crs = terra::crs(.map))
      # Define probability of detection & non detection
      dist <- terra::distance(.map, v)
      p1 <- terra::app(dist, function(x) plogis(m$receiver_alpha + m$receiver_beta * x))
      p0 <- 1 - p1
      # Mask by land
      p1 <- terra::mask(p1, .map) |> readAll()
      p0 <- terra::mask(p0, .map) |> readAll()
      list(nondetection = p0, detection = p1)
    })
    names(likelihoods) <- as.character(moorings$sensor_id)
    
    # Iteratively create plots
    cl_lapply(.steps, .cl = .cl, .fun = function(t) {
      
      # Define time-specific datasets
      # t = 1
      # print(t)
      tstamp <- timeline[t]
      acc    <- acoustics[timestamp == tstamp, ]
      det    <- detections[timestamp == tstamp, ]
      cinfo  <- .input$.yobs$ModelObsContainer[timestamp == tstamp, ]
      states <- .output$states[timestamp == tstamp, ]
      
      # Set up plot
      png(file.path(frames, paste0(t, ".png")), 
          height = 5, width = 6, res = 200, units = "in")
      pp <- par(mfrow = c(1, 2))
      
      
      #### Plot detection time series ------------------------------------------
      
      plot(detections$timestep, detections$sensor_id, 
           xlab = "Time (steps)", ylab = "Receiver ID", 
           xlim = c(0, length(timeline)))
      mtext(side = 3, t, line = 0, font = 2)
      if (nrow(det) > 0L) {
        points(det$timestep, det$sensor_id, col = "green", lwd = 2)
      }
      abline(v = t, col = "red", lty = 3)
    
      
      #### Plot map ------------------------------------------------------------
      
      # (optional) Compute likelihood surface
      likelihood <- 
        lapply(split(acc, seq_len(nrow(acc))), function(m) {
        if (m$obs == 1) {
          likelihoods[[as.character(m$sensor_id)]]$detection
        } else {
          likelihoods[[as.character(m$sensor_id)]]$nondetection
        }
      }) |> 
        terra::rast() |> 
        terra::app("prod", na.rm = TRUE)
      
      # Plot base map
      terra::plot(likelihood,
                  col = scales::alpha(terra::map.pal("viridis", 100), 0.8))
      
      # Add states 
      points(states$x, states$y, pch = ".", col = "blue")
      
      # Add operational receivers
      # * Detection = green; non detection = grey 
      points(acc$receiver_x, acc$receiver_y, 
             col = ifelse(acc$obs == 1L, "green", "black"), 
             lwd = ifelse(acc$obs == 1L, 2, 1))
      
      # Add detection containers 
      if (nrow(cinfo) > 0L) {
        points(cinfo$centroid_x, cinfo$centroid_y, lwd = 1.5, col = "red")
        cbind(cinfo$centroid_x, cinfo$centroid_y) |>
          terra::vect() |>
          terra::buffer(width = cinfo$radius) |> 
          terra::lines(col = "red", lwd = 2) 
      }

      
      #### Global updates ------------------------------------------------------
      
      # Title
      det_lbl  <- length(which(acc$obs == 1L))
      move_lbl <- paste0(strwrap(.input$.model_move, width = 50), collapse = "\n")
      main_lbl <- paste0("Time: ", t, "; Detection(s):", det_lbl, "; \n", move_lbl)
      mtext(side = 3, main_lbl, font = 2, adj = 0.1, line = -3, outer = TRUE, cex = 0.75)
      
      # Close plot
      dev.off()
      
    })
    
    # Make animation 
    tictoc::tic()
    input   <- gtools::mixedsort(list.files(frames, full.names = TRUE))
    output  <- file.path(mp4, "ani.mp4")
    av::av_encode_video(input, output, framerate = 150)
    unlink(frames, recursive = TRUE)
    tictoc::toc()
    
    # Open animation (on MacOS)
    system(paste("open", shQuote(output)))
    invisible(NULL)
  }
  
}

