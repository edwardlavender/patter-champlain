###########################
###########################
#### Core constructor function for analyses

# `.sim` must contain:
# * shape, scale, mobility, phi for model_move_trout()
# * file_detections, receiver_alpha, receiver_beta, receiver_gamma for model_obs_champlain()
# * file_output for output files

constructor_ac_core <- function(.sim, .datasets, .verbose, ...) {
    
    # Checks 
    stopifnot(length(list(...)) == 0L)
    proj.build::check_names(.sim, c( "mobility", "file_detections", "file_output"))
    
    # Enable testing & define tuning settings 
    test <- FALSE
    n_particle_filter <- 5e4
    n_particle_smo    <- 2e3L
    n_sim_smo         <- 100L
    if (test) {
      warning("test = TRUE!", immediate. = TRUE)
      n_particle_filter <- 5e3L
      n_particle_smo    <- 100L
      n_sim_smo         <- 30L
    }
    n_record <- max(c(1000L, n_particle_smo))
    
    # Read datasets
    proj.build::check_names(.datasets, c("cap_recap", "moorings"))
    cap_recap  <- .datasets$cap_recap
    moorings   <- .datasets$moorings
    detections <- qs::qread(.sim$file_detections)
    map_bbox   <- qs::qread(here_input("map-bbox.qs"))
    
    # Define timeline
    timeline <- get_dataset_timeline(.sim = .sim)
    if (test) {
      timeline <- timeline[1:200L]
    }
    
    # Define movement model
    state      <- state_trout()
    model_move <- model_move_trout(.sim)
    
    # (optional) Define initial states
    # * If provided, cap_recap should be a two-row data.table with map_value, x, y 
    # * This contains the [1] capture and [2] recapture locations (used in simulations)
    if (!is.null(cap_recap)) {
      stopifnot(nrow(cap_recap) == 2L)
      proj.build::check_names(cap_recap, c("map_value", "x", "y"))
      map_value <- x <- y <- NULL
      cap       <- cap_recap[1, list(map_value, x, y)]
      recap     <- cap_recap[2, list(map_value, x, y)]
      xinit_fwd <- model_move_xinit(.xinit = cap, .n_particle = n_particle_filter)
      xinit_bwd <- model_move_xinit(.xinit = recap, .n_particle = n_particle_filter)
      
      # Visual checks of initial locations (for debugging)
      # map  <- terra::rast(here_input("map.tif"))
      # path <- qs::qread(here_input_sim("paths.qs"))[path_id == .sim$unit_id, ]
      # terra::plot(map)
      # points(xinit_fwd$x, xinit_fwd$y)
      # points(path$x[1], path$y[1], pch = ".")
      # points(xinit_bwd$x, xinit_bwd$y)
      # points(path$x[nrow(path)], path$y[nrow(path)], pch = ".")
      
    } else {
      xinit_fwd <- NULL
      xinit_bwd <- NULL
    }
    
    # (optional) Assemble capture/recapture containers
    containers <- xinit_containers <- NULL 
    if (!is.null(cap_recap)) {
      # Set capture/recapture locations
      # * If test = TRUE, we use a timeline[1:n]
      # * So we can set x1 but we must set xT to NULL
      x1 <- cap
      xT <- recap
      if (test) {
        xT <- NULL
      }
      # Assemble containers 
      containers <- xinit_containers <- 
        assemble_xinit_containers(.timeline = timeline, 
                                  .xinit    = list(forward = x1, backward = xT), 
                                  .radius   = .sim$mobility, 
                                  .mobility = .sim$mobility,
                                  .map      = map_bbox)
    } 
    
    # Assemble acoustic observations
    moorings  <- model_obs_champlain(.moorings     = moorings,
                                      .pars        = .sim, 
                                      .as_ModelObs = FALSE)
    acoustics <- assemble_acoustics(.timeline    = timeline, 
                                     .detections = detections, 
                                     .moorings   = moorings)
    yobs_fwd <- yobs_bwd <- list(ModelObsAcousticLogisTrunc = acoustics)
    # terra::plot(map)
    # points(acoustics$receiver_x, acoustics$receiver_y)
    
    # Assemble acoustic containers
    if (length(which(acoustics$obs == 1L)) > 2L) {
      acoustic_containers <- assemble_acoustics_containers(.timeline = timeline, 
                                                           .acoustics = acoustics, 
                                                           .mobility = .sim$mobility, 
                                                           .map = map_bbox)
      if (!is.null(xinit_containers)) {
        containers <- assemble_containers(xinit_containers, acoustic_containers)
      } else {
        containers <- acoustic_containers
      }
    }
    
    # Collate containers
    if (use_containers(.containers = containers, .direction = "forward")) {
      yobs_fwd$ModelObsContainer <- containers$forward
    } 
    if (use_containers(.containers = containers, .direction = "backward")) {
      yobs_bwd$ModelObsContainer <- containers$backward
    }
    
    # Define arguments for forward filter run
    # * With batching, we can keep .collect = TRUE to record diagnostics/callstats
    args_fwd <- list(.timeline   = timeline, 
                     .state      = state,
                     .model_move = model_move, 
                     .yobs       = yobs_fwd, 
                     .n_particle = n_particle_filter, 
                     .n_record   = n_record,
                     .direction  = "forward", 
                     .batch      = particle_batch(.sim = .sim, .type = "fwd"),
                     .collect    = TRUE,
                     .verbose    = .verbose, 
                     .progress   = julia_progress(enabled = test))
    stopifnot(all(names(args_fwd) %in% names(formals(pf_filter))))
    
    # (optional) Interactive debugging of filter convergence failures
    if (FALSE) {
      
      # Load additional packages
      library(ggplot2)
      
      # Run filter for troublesome individual/time
      if (FALSE) {
        
        # Run filter
        args_fwd$.n_particle <- 5e4L
        args_fwd$.n_record   <- 1000L
        args_fwd$.batch      <- NULL
        args_fwd$.progress   <- julia_progress(enabled = TRUE)
        fwd <- do.call(pf_filter, args_fwd)
        # > unit_id  individual_id    time_id 
        # > 16       24321         2016-01-01
        # > Weights from filter (1 -> 18762) are zero at time 12806:returning outputs from 1:12806. Note that all (log) weights at 12806 are -Inf.
        # > @ 2014-11-22 17:20:00 
        
        # Record inputs/outputs
        qs::qsave(.sim, here_data("debug", "sim.qs"))
        qs::qsave(moorings, here_data("debug", "moorings.qs"))
        qs::qsave(args_fwd, here_data("debug", "input.qs"))
        qs::qsave(fwd, here_data("debug", "output.qs"))

      } else {
        fwd <- qs::qread(here_data("debug", "output.qs"))
      }
      
      # Identify problematic time step/stamps
      tdt    <- data.table(timestep = 1:length(timeline), timestamp = timeline)
      tstep  <- 12806
      tstamp <- timeline[tstep]
      
      # Visualise detection time series
      debug_detections <-
        args_fwd$.yobs$ModelObsAcousticLogisTrunc |> 
        lazy_dt() |> 
        mutate(timestep = tdt$timestep[match(timestamp, tdt$timestamp)]) |>
        filter(obs == 1L) |> 
        as.data.table() 
      debug_gg <- 
        debug_detections |>
        ggplot() + 
        geom_point(aes(timestep, sensor_id)) +
        geom_vline(aes(xintercept = tstep)) 
      plotly::ggplotly(debug_gg)
      
      # Check detections before/after convergence failure
      # 12770, sensor id 56
      # 12811, sensor id 43
      
      # Identify detection containers
      debug_container <- 
        args_fwd$.yobs$ModelObsContainer |> 
        lazy_dt() |> 
        mutate(timestep = tdt$timestep[match(timestamp, tdt$timestamp)]) |>
        as.data.table()
      
      # Check detection containers correctly shrink in the 5 time steps before failure: ok.
      debug_container[timestamp >= (tstamp - 5 * 60), ] |> head(10)
      # timestamp   obs sensor_id centroid_x centroid_y radius
      # <POSc> <int>     <int>      <num>      <num>  <num>
      #   1: 2014-11-22 17:16:00     1        43    1793697   927281.8   9512
      # 2: 2014-11-22 17:18:00     1        43    1793697   927281.8   9296
      
      # Container size in moment before detection should be .sim$receiver_gamma + .sim$mobility: ok. 
      debug_container[timestep == 12810, ]
      
      # Plot receivers (not in Julia session on Linux)
      map      <- terra::rast(here_input("map.tif"))
      moorings <- qs::qread(here_data("debug", "moorings.qs"))
      terra::plot(map)
      r0 <- moorings[receiver_id == 56, .(receiver_x, receiver_y)]
      points(r0)
      r1 <- moorings[receiver_id == 43, .(receiver_x, receiver_y)]
      points(r1)
      
      # Compute distance/movement speed between receivers: 
      # 12770, sensor id 56
      # 12811, sensor id 43
      debug_dist <- patter:::dist_2d(as.matrix(r0), as.matrix(r1)) 
      # average travel distance required per time step 
      # > 78 m per time step is not unreasonable
      debug_dist / (12811 - 12770) 
      
      # Create animation (not in Julia session on Linux)
      # unlink(here_fig("debug", qs::qread(here_data("debug", "sim.qs"))$index, "frames"), recursive = TRUE)
      ani(.sim      = qs::qread(here_data("debug", "sim.qs")), 
          .map      = terra::rast(here_input("map.tif")), 
          .moorings = qs::qread(here_data("debug", "moorings.qs")), 
          .start    = 12500L, 
          .end      = 12811, 
          .input    = qs::qread(here_data("debug", "input.qs")), 
          .output   = qs::qread(here_data("debug", "output.qs")), 
          .cl = 50L
          )
      
      # Take home message
      # > Particle animation suggests insufficient directness (Normal(0, 1.8)) in movement model

    }
    
    # Prepare smoothing outputs unless .sim$smooth = FALSE explicitly specified 
    if (!rlang::has_name(.sim, "smooth") || (rlang::has_name(.sim, "smooth") & .sim$smooth)) {
      
      # Define arguments for backward filter run
      args_bwd            <- args_fwd
      args_bwd$.yobs      <- yobs_bwd
      args_bwd$.direction <- "backward"
      args_bwd$.batch     <- particle_batch(.sim = .sim, .type = "bwd")
      
      # Define smoother arguments
      # * Note .collect = TRUE is required for particle_success()
      args_smo <- list(.n_particle = n_particle_smo, 
                       .n_sim = n_sim_smo, 
                       .cache = TRUE, 
                       .batch = particle_batch(.sim = .sim, .type = "smo"),
                       .progress = julia_progress(enabled = test), 
                       .collect = TRUE, 
                       .verbose = .verbose)
      stopifnot(args_smo$.collect)
      
      # Checks
      stopifnot(all(names(args_bwd) %in% names(formals(pf_filter))))
      stopifnot(all(names(args_smo) %in% names(formals(pf_smoother_two_filter))))
      
    } else {
      # Set args_bwd and args_smo to NULL to suppress backward filter/smoother
      args_bwd <- args_smo <- NULL
    }
    
    # Collate arguments
    # * estimate_coord_particle() requires the following arguments:
    # - `forward`
    # - `backward`  (if smoothing desired, NULL otherwise)
    # - `smooth`    (if smoothing desired, NULL otherwise)
    list(forward = args_fwd, backward = args_bwd, smooth = args_smo, verbose = .verbose)
    
  }
  

# Call constructor_ac_core()
call_constructor_ac_core <- function(.sim, 
                                     .timeline, .cap_recap, .moorings, .detections, .map_bbox,
                                     .verbose, ...) {
  datasets <- list(cap_recap  = .cap_recap, moorings   = .moorings)
  args <- list(.sim = .sim, .datasets = datasets, .verbose = .verbose, ...)
  do.call(constructor_ac_core, args)
}


###########################
###########################
#### User-facing constructors 

# These functions wrap the core routine, constructor_ac_core()

constructor_ac_sim <- function(.sim, .datasets, .verbose, ...) {
  
  # Checks
  stopifnot(length(.datasets) == 0L)
  
  # Read datasets
  # * For convenience, cap_recap for all individuals is stored in one list
  cap_recap  <- qs::qread(here_input_sim("xinits.qs"))[[.sim$unit_id]]
  moorings   <- qs::qread(here_input_sim("moorings.qs"))
  
  # Call constructor_ac_core()
  call_constructor_ac_core(.sim        = .sim, 
                           .cap_recap  = cap_recap, 
                           .moorings   = moorings,
                           .verbose    = .verbose, ...)

}

constructor_ac_real <- function(.sim, .datasets, .verbose, ...) {
  
  # Checks
  stopifnot(length(.datasets) == 0L)
  
  # Read type-specific datasets
  cap_recap  <- NULL
  moorings   <- qs::qread(here_input_real("moorings.qs"))
  
  # Run constructor_ac_core()
  call_constructor_ac_core(.sim        = .sim, 
                           .cap_recap  = cap_recap, 
                           .moorings   = moorings,
                           .verbose    = .verbose, ...)
  
}


###########################
###########################
#### pf_filter_loglik_optim() constructor

constructor_pf_filter_loglik_optim <- function(.sim, .datasets, .verbose, ...) {
  # pf_filter_loglik_optim() only requires a .sim argument
  list(.sim = copy(.sim))
}
