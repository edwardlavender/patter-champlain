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

# particle_placeholder
# * This function is a placeholder for estimate_coord_particle()
# * It is used to check the constructor functions work for all datasets
particle_placeholder <- function(...) {
  list(forward = NULL, backward = NULL, smooth = NULL)
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

# particle convergence (post-hoc check)
# * This function expects a list from estimate_coord_particle() 
#   with $forward, $backward and $smooth elements
particle_convergence <- function(l) {
  # Iterate over each element & check $callstats$convergence
  convergences <- sapply(c("forward", "backward", "smooth"), function(direction) {
    convergence <- FALSE
    if (rlang::has_name(l, direction) && !is.null(l[[direction]]$callstats)) {
      convergence <- l[[direction]]$callstats$convergence
    }
    convergence
  }) 
  data.table(forward = convergences[1], backward = convergences[2], smooth = convergences[3])
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
  if (!rlang::has_name(out$smooth, "callstats") || !out$smooth$callstats$convergence) {
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

# Iterative workflow
# * constructor_ac_sim() implements batching
# * For each iteration, we should collate the estimated coordinates across batches
# * This function expects the .iteration and .timeline
# * A subsetted timeline may be required if used for testing
cl_lapply_particle_collate <- function(.iteration) {
  
  # Check inputs
  # * file_output required by particle_collate -> particle_batch()
  proj.build::check_names(.iteration, c("folder_coord", "file_output"))
  stopifnot(all(file.exists(.iteration$file_output)))
  .iteration <- copy(.iteration)
  
  # For each simulation, collate particles across batches & cleanup smo-{i}.jld2 files
  convergence <- cl_lapply(split(.iteration, seq_len(nrow(.iteration))), function(.sim) {
    
    # Define timeline for simulation
    timeline <- get_dataset_timeline(.sim = .sim)
    # timeline <- timeline[1:500L]
    
    # Collate particles across batches and write file_coord
    convergence <- particle_collate(.sim = .sim,
                                    .timeline = timeline)
    
    # (optional) Clean up smo-{i}.jld2 files to save space
    if (FALSE) {
      batches <- list.files(.sim$folder_coord, 
                            pattern = "^smo-\\d+\\.jld2$", 
                            full.names = TRUE)
      unlink(batches)
    } else {
      warning("Clean up of old ^smo-\\d+\\.jld2$ files is currently suppressed.", immediate. = TRUE)
    }
    
    # Return convergence
    convergence
  })
  
  unlist(convergence)
}

if (!patter:::os_linux() | (patter:::os_linux() & !patter:::julia_session())) {
  
  # Compute area spanned by 95 % of particles via 2D histogram
  # (Computation via ks::kde() and ks::contourSizes() is too slow)
  particle_hr <- function(.map, .coord) {
    
    # Compute area of grid cell (assuming UTM grid)
    A <- prod(terra::res(.map))
    
    # Compute average area spanned by 95 % of particles
    .coord |>
      lazy_dt() |>
      select(timestep, x, y) |>
      # filter(timestep %in% seq(1, max(timestep), by = 100)) |> 
      # Discretise coordinates
      mutate(id = terra::cellFromXY(.map, cbind(.data$x, .data$y)),
             x = terra::xFromCell(.map, .data$id),
             y = terra::yFromCell(.map, .data$id)) |>
      # Assign equal weights (marks)
      group_by(.data$timestep) |>
      mutate(mark = 1 / n()) |>
      ungroup() |>
      # Calculate the total weight of each location within time steps (2D histogram)
      group_by(.data$timestep, .data$id) |>
      summarise(mark = sum(.data$mark)) |>
      ungroup() |>
      # Calculate for each timestep the area containing 95 % of probability mass
      # * Sum weights
      # * Identify the number of cells required to put us above 95 % probability mass * A
      group_by(.data$timestep) |> 
      arrange(desc(mark), .by_group = TRUE) |> 
      mutate(cmark = cumsum(mark)) |>
      summarise(area = which(cmark >= 0.95)[1] * A ) |>
      ungroup() |> 
      # Compute mean area over all time steps
      summarise(area_mean = mean(area)) |>
      ungroup() |>
      pull(area_mean)
    
  }
  
}
