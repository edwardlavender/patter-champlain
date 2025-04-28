#### Define custom estimate_coord constructor function 
# .sim must contain:
# * file_detections, 
# * shape, scale, mobility, phi for model_move_trout()
# * file_detections, receiver_alpha, receiver_beta, receiver_gamma for model_obs_champlain()
# * file_output for output files
constructor_ac_sim <- function(.sim, .datasets, .verbose, ...) {
  
  # Checks 
  stopifnot(length(list(...)) == 0L)
  proj.build::check_names(.sim, c("file_detections", "mobility", "file_output"))
  stopifnot(file.exists(.sim$file_detections))
  
  # Enable testing & define tuning settings 
  test              <- TRUE
  n_particle_filter <- 2.5e4
  n_particle_smo    <- 1e3L
  n_sim_smo         <- 100L
  if (test) {
    warning("test = TRUE!", immediate. = TRUE)
    n_particle_filter <- 5e3L
    n_particle_smo    <- 100L
    n_sim_smo         <- 30L
  }
  
  # Define timeline
  timeline <- qs::qread(here_input_sim("timeline.qs"))
  if (test) {
    timeline <- timeline[1:250L]
  }
  
  # Define movement model
  state      <- state_trout()
  model_move <- model_move_trout(.sim)
  
  # TO DO (optional) Define initial states
  #
  
  # TO DO (optional) Assemble capture/recapture containers
  #
  
  # Assemble acoustic observations
  # (We could also read datasets for .sim$unit_id from file)
  moorings   <- qs::qread(here_input_sim("moorings.qs"))
  moorings   <- model_obs_champlain(.moorings = moorings,
                                    .pars = .sim, 
                                    .as_ModelObs = FALSE)
  detections <- qs::qread(.sim$file_detections)
  acoustics  <- assemble_acoustics(.timeline   = timeline, 
                                   .detections = detections, 
                                   .moorings   = moorings)
  yobs_fwd <- yobs_bwd <- list(ModelObsAcousticLogisTrunc = acoustics)
  
  # Assemble acoustic containers
  if (length(which(acoustics$obs == 1L) > 2L)) {
    map_bbox   <- qs::qread(here_input("map-bbox.qs"))
    containers <- assemble_acoustics_containers(.timeline = timeline, 
                                                .acoustics = acoustics, 
                                                .mobility = .sim$mobility, 
                                                .map = map_bbox)
    yobs_fwd$ModelObsContainer <- containers$forward
    yobs_bwd$ModelObsContainer <- containers$backward
  }
  
  # Define arguments for forward filter run
  # * With batching, we can keep .collect = TRUE to record diagnostics/callstats
  args_fwd <- list(.timeline   = timeline, 
                   .state      = state,
                   .model_move = model_move, 
                   .yobs       = yobs_fwd, 
                   .n_particle = n_particle_filter, 
                   .direction  = "forward", 
                   .batch      = particle_batch(.sim = .sim, .type = "fwd"),
                   .collect    = TRUE,
                   .verbose    = .verbose, 
                   .progress   = julia_progress(enabled = test))
  
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
  stopifnot(all(names(args_fwd) %in% names(formals(pf_filter))))
  stopifnot(all(names(args_bwd) %in% names(formals(pf_filter))))
  stopifnot(all(names(args_smo) %in% names(formals(pf_smoother_two_filter))))
  
  # Collate filter arguments
  # * particle_algorithm() requires the following arguments:
  # - `forward`
  # - `backward`  (if smoothing desired, NULL otherwise)
  # - `smooth`    (if smoothing desired, NULL otherwise)
  list(forward = args_fwd, backward = args_bwd, smooth = args_smo, verbose = .verbose)
  
}

# pf_filter_loglik_optim() constructor
constructor_pf_filter_loglik_optim <- function(.sim, .datasets, .verbose, ...) {
  # pf_filter_loglik_optim() only requires a .sim argument
  list(.sim = copy(.sim))
}