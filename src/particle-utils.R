# Use containers
use_containers <- function(.containers, .direction) {
  !is.null(.containers) && 
    rlang::has_name(.containers, .direction) && 
    !is.null(.containers[[.direction]]) && 
    nrow(.containers[[.direction]]) > 0L
}


# Define batches for particle algorithm
particle_batch <- function(.sim, .type = c("fwd", "bwd", "smo")) {
  .type <- match.arg(.type)
  proj.build::check_names(.sim, "file_output")
  file.path(dirname(.sim$file_output), paste0(.type, "-", 1:3, ".jld2"))
}

# Determine success of particle algorithms
# * x is the output of particle algorithms
particle_success <- function(x) {
  # Set success = FALSE if the smoothing element is empty
  # * This is due to a convergence failure during smoothing 
  # * (since .collect = TRUE for the smoother)
  if (is.null(x$smooth)) {
    FALSE
  } else {
    # Otherwise, extract x$smooth$callstats$convergence
    x$smooth$callstats$convergence
  }
}

# Cleanup after particle algorithms
# * We delete fwd/bwd batch files after use
#   to minimise storage requirements (on each iteration)
particle_cleanup <- function(.sim, .cl) {
  # Delete fwd files
  fwd <- particle_batch(.sim = .sim, .type = "fwd")
  sapply(fwd, function(f) {
    if (file.exists(f)) {
      unlink(f)
    }
  })
  # Delete bwd files
  bwd <- particle_batch(.sim = .sim, .type = "bwd")
  sapply(bwd, function(f) {
    if (file.exists(f)) {
      unlink(f)
    }
  })
  invisible(NULL)
}

# Collate batches in R
particle_collate <- function(.sim, .timeline) {
  
  # Check names
  # * We run particle algorithms with batching
  # * Each output file contains the output of estimate_coord_particle():
  # * list(forward = fwd, backward = bwd, smooth = smo)
  # * Each element is list(states = NULL, diagnostics, callstats)
  # * We will create a file_coord files
  check_names(.sim, c("file_diag", "file_coord"))
  
  # (optional) Skip .sim$file_diag files that don't exist
  # * This is relevant during testing
  
  # Open .sim$file_diag
  out <- qs::qread(.sim$file_diag)
  if (!out$smooth$callstats$convergence) {
    return(FALSE)
  }
  
  # Define timeline 
  julia_assign("timeline", .timeline)
  
  # Define batch files (for smoother)
  batch <- particle_batch(.sim = .sim, .type = "smo")
  julia_assign("bbb", batch)
  
  # Collate states in Julia 
  julia_command('
    smo_states = hcat([f["xsmo"] for f in map(jldopen, bbb)]...);')
  
  # Collect states in R
  out$smooth$states <- julia_eval('
    Patter.r_get_states(smo_states, collect(1:length(timeline)), timeline);
    ')
  
  # Write coord.qs file
  qs::qsave(out$smooth, .sim$file_coord)
  TRUE
}

