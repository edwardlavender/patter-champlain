###########################
###########################
#### analysis-patter.R

#### Aims
# 1) This script provides a generic workflow for analysing simulated/real-world datasets with patter

#### Prerequisites
# 1) Run simulations


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Set global options
patter::julia_connect()

#### Load essential packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(glue)
library(JuliaCall)
library(patter)
library(patter.workflows)
library(proj.verse)
library(tictoc)
files_source_r(here_src())
expect_no_geospatial()


###########################
###########################
#### Select analysis type

# Select analysis type ("sim", "real")
# analysis <- "sim"
# analysis <- "real"

# Select iterations by mobility (162, 216, 270)
analysis_mobility <- 216

# Set development mode
# * Use one core
# * Set maps 
dev <- TRUE

# (optional) commandArgs() override for server deployments
cmd_args <- commandArgs(trailingOnly = TRUE)
if (length(cmd_args) > 0L) {
  stopifnot(length(cmd_args) == 3L)
  analysis          <- cmd_args[1]
  analysis_mobility <- as.numeric(cmd_args[2])
  dev               <- as.logical(cmd_args[3])
}

# Check input settings
stopifnot(analysis %in% c("sim", "real"))
print(glue("Arguments: analysis = '{analysis}'; analysis_mobility = {analysis_mobility}; dev = {dev}."))

#### Define analysis-specific routines
here_input_analysis       <- switch_here_input_analysis(analysis)
here_output_analysis_main <- switch_here_output_analysis_main(analysis)
constructor_ac_analysis   <- switch_constructor_ac_analysis(analysis)

#### Define analysis-specific data
iteration <- qs::qread(here_input_analysis("iteration-patter.qs"))
nrow(iteration)


###########################
###########################
#### Estimate coordinates

#### (optional) Reset directories
if (FALSE) {
  # unlink(dirname(iteration$folder_coord), recursive = TRUE)
  unlink(iteration$folder_coord, recursive = TRUE)
  dirs.create(iteration$folder_coord)
}

#### Select iterations (by mobility)
stopifnot(!any(duplicated(iteration$index)))
table(iteration$mobility)
iteration <- iteration[mobility == analysis_mobility, ]
iteration[, file_diag := file.path(folder_coord, "diagnostics.qs")]
iteration[, file_output := file_diag]
# (optional) Further subset for testing
if (TRUE) {
  iteration <- iteration[sensitivity == "best", ]
  iteration <- iteration[1:min(c(.N, 100L)), ]
}
nrow(iteration)
stopifnot(nrow(iteration) > 0L)

#### Set maps
if (dev) {
  # This is implemented below on each node. 
  set_map(here_input("map.tif"))
  set_vmap(.vmap = here_input("vmap", iteration$mobility[1], "vmap.tif"))
}

#### Set logs
if (dev) {
  log.txt <- TRUE
} else {
  dir.create(here_output_analysis_main("logs", "R"))
  log.txt <- here_output_analysis_main("logs", "R", paste0("log-", iteration$mobility[1], ".txt"))
  # unlink(log.txt)
}

#### Setup cluster
# Quick check on Linux
if (FALSE) {
  # On siam-linux20, julia_connect() may fail on individual nodes
  # Updating package versions can solve errors 
  # This code is a quick check to see if it works with current package versions
  cl <- parallel::makeCluster(2L)
  parallel::clusterEvalQ(cl, {
    patter::julia_connect(.socket = TRUE)
  })
  pbapply::pblapply(1:2, function(i) {
    JuliaCall::julia_eval(glue::glue('{i} + {i}'))
  }, cl = cl)
  parallel::stopCluster(cl)
}
# Setup cluster
if (!dev) {
  # Define number of workers
  ncl <- min(c(nrow(iteration), 100L, parallel::detectCores() - 1L))
  # Define required memory for export
  lobstr::mem_used() * ncl
  # Initialise cluster
  cl  <- parallel::makeCluster(ncl)
  cl_init(iteration = iteration, cl = cl, varlist = ls())
} else {
  ncl <- 1L
  cl  <- NULL
}

#### Estimate coordinates: time trials (simulations)
# iteration[1, ] on SIA-LAVENDED-M
# * iteration[1, ], 1 thread, 2.5e4 filter particles, 1e3 smoothing particles, 1e2 smoothing sims
# * 26.8 min XX min: 6.95 min (filter) + 6.82 min (filter) + 11.98 min (smoother)
# * Warning: All smoothing weights (from xbwd[k, t] to xfwd[j, t - 1]) are zero at 8 time step(s) (0.04 %).
# iteration[1:90, ] on siam-linux20
# * 1 hour, 90 cl

#### Estimate coordinates
# TO DO In patter.workflows, update .verbose for parallelisation
# debug(constructor_ac_core)
# iteration <- iteration[1:1L, ]
print(glue("Using {ncl} core(s) for {nrow(iteration)} iteration row(s) (mobility = {iteration$mobility[1]})."))
coord_list <- 
  cl_lapply_workflow(.iteration   = iteration,
                     .datasets    = list(),
                     .constructor = constructor_ac_analysis, 
                     .algorithm   = estimate_coord_particle, 
                     .success     = particle_success, 
                     .cleanup     = particle_cleanup,
                     .cl          = cl,
                     .verbose     = log.txt)

#### Progress:
# sim-1   : TO DO
# sim-2   : TO DO
# sim-3   : TO DO
# real-1  : TO DO
# real-2  : TO DO
# real-3  : TO DO

#### Review file sizes (storage requirements) 
# Each iteration records callstats, diagnostics, smo (particles)
# We record 1000 particles (4 state dimensions) for _up to_ one month
# This is ~714 MB per iteration
nr <- 1000L
ns <- 4L
nt <- length(qs::qread(here_input_sim("timeline.qs")))
(nr * ns * nt * 8) / 1e6
# Check output file sizes for iteration[1, ]: 
# * "callstats.qs"                           : < 0.1 MB
# * "diagnostics.qs"                         : 0.7 MB
# * "smo-1.jld2", "smo-2.jld2", "smo-3.jld2" : 238 MB each 
files <- list.files(iteration$folder_coord[1], full.names = TRUE)
sapply(files, file.size) / 1e6
# As an approximation of the total storage requirements (GB) is:
# > 107 GB for 150 simulations
# > 2.35 TB for 3290 real-world analyses 
nrow(qs::qread(here_input_analysis("iteration-patter.qs"))) * 714 / 1000
# TO DO Consider reducing storage requirements:
# * Batch particles appropriately & output summary statistics rather than particles
# * Periodically summarise and cleanup particles (run iteration in batches)

#### Collate coordinates across batches
if (FALSE) {
  
  #### An explanation of callstats.qs
  # All iterations should output a callstats.qs file
  # * This is derived from cl_lapply:::workflow()
  list.files(dirname(iteration$file_output)[1])
  iteration[, file_callstats := file.path(dirname(file_output), "callstats.qs")]
  stopifnot(all(file.exists(iteration$file_callstats)))
  qs::qread(iteration$file_callstats[1])

  #### An explanation of file_output
  # All iterations should also output a iteration$file_output (diagnostics.qs) file
  # * This contains the output of estimate_coord_particle()
  # * $forward$states = NULL, forward$diagnostics, foward$callstats
  # * $backward$states = NULL, backward$diagnostics, backward$callstats
  # * $smooth$states = NULL, smooth$diagnostics, smooth$callstats
  # > The states elements are NULL because we implemented batching
  # > Filter/smoother files are written to file
  # > We clean up fwd/bwd files on the fly & only retain smo states
  # qs::qread(iteration$file_output[1])
  table(file.exists(iteration$file_output))

  #### Check convergence 
  # Compute convergence of filters and smoother
  convergence_dt <- cl_lapply_iteration_file(
    iteration, 
    .file = "file_output", 
    .fun = function(.sim, .input) {
      particle_convergence(.input)
    }) |> rbindlist()
  # Check the number of forward/backward filter runs that converged
  head(convergence_dt)
  table(convergence_dt$forward & convergence_dt$backward)
  # Check the number of smoothing runs that converged
  # * This relies on a high threshold of > 95 % time steps with 'proper' smoothing
  table(convergence_dt$smooth)
  # Examine proportion of time steps with 'proper' smoothing for each run
  smooth_ess_prop <- cl_lapply_iteration_file(
    iteration, 
    .file = "file_output", 
    .fun = function(.sim, .input) {
      out <- NULL
      if (rlang::has_name(.input$smooth, "diagnostics")) {
        n0 <- nrow(.input$smooth$diagnostics)
        n1 <- length(which(!is.na(.input$smooth$diagnostics$ess)))
        out <- n1 / n0
      }
      out 
    }) |> unlist()
  hist(smooth_ess_prop, breaks = 100)
  table(smooth_ess_prop > 0.9)
  utils.add::basic_stats(smooth_ess_prop)
  # Summarise smoother ESS for successful smoothing runs
  ess_mean <- cl_lapply_iteration_file(
    iteration,
    .file = "file_output", 
    .fun = function(.sim, .input) {
      if (rlang::has_name(.input$smooth, "diagnostics")) {
        return(mean(.input$smooth$diagnostics$ess, na.rm = TRUE))
      }
      return(NA)
    }, .combine = unlist) 
  hist(ess_mean, xlim = c(0, 2000))
  # Results (100 real iterations)
  # * With 2.5e4 filter particles & 1,000 smoothing particles:
  # - 1 hour (100 cl)
  # - 73/100 forward/backward filter successes
  # - 52/73 forward/backward & smoothing successes (95 % threshold)
  # - 58/73 forward/backward & smoothing successes (90 % 'patter-flapper' threshold)
  # * With 5e4 filter particles & 2,000 smoothing particles:
  # - 3 hours (100 cl)
  # - 76/100 forward/backward filter successes
  # - 54/76 forward/backward & smoothing successes (95 % threshold)
  # - 60/76 forward/backward & smoothing successes (90 % 'patter-flapper' threshold)
  # * For convergence, there is little benefit in boosting the number of particles
  
  # Collate smoothed states 
  # * ETA: ~5 s per row on 1 cl (07m 03s for 84 rows)
  iteration[, file_coord := file.path(folder_coord, "coord.qs")]
  success <- cl_lapply_particle_collate(.iteration = iteration[file.exists(file_output), ])
  
  # Check all files were successfully created
  # * For file structure, see below
  table(unlist(convergence))
  
}

#### Examine coordinates
if (FALSE) {
  
  #### Examine example file
  # Each file_coord file is an standard pf_particles output from the smoother
  # $states, $diagnostics, $callstats (for smoothing)
  list.files(iteration$folder_coord)
  eg <- qs::qread(iteration$file_coord[1])
  summary(eg)
  # View(eg)

  #### Visual checks
  
  # Visually check 'best' maps for a few individuals
  map      <- terra::rast(here_input("map.tif"))
  unit_ids <- unique(iteration$unit_id)[1:9L]
  lapply_qplot_coord(iteration[unit_id %in% unit_ids & sensitivity == "best", ], 
                     .datasets = list(coordinates = function(smo) smo$states),
                     .map = map, 
                     .n_plot = 9L)

  # Visually check 'sensitivity' for a few individuals
  lapply_qplot_coord(iteration[unit_id == 1L, ], 
                     .datasets = list(coordinates = function(smo) smo$states),
                     .map = map)
  lapply_qplot_coord(iteration[unit_id == 2L, ], 
                     .datasets = list(coordinates = function(smo) smo$states),
                     .map = map)
  
  #### Quantitative checks
  # (A) Check distances between individual/receiver @ moment of detections
  cl_lapply_iteration_file(
    iteration[1:5], 
    .file = "file_coord", 
    .fun = function(.sim, .input) {
      # Define states/detections
      states     <- .input$states
      detections <- qs::qread(.sim$file_detections)
      # Compute distances 
      distances  <- 
        right_join(states, detections, by = "timestamp", relationship = "many-to-many") |> 
        select(x, y, receiver_x, receiver_y) |> 
        mutate(dist = terra::distance(cbind(x, y), cbind(receiver_x, receiver_y), 
                                      lonlat = FALSE, pairwise = TRUE)) |>
        as.data.table()
      # Validate distances
      stopifnot(all(distances$dist < .sim$receiver_gamma))
    })
  
}

#### Examine diagnostics
if (FALSE) {
  
  #### Summarise convergence 
  # (see above)
  
  #### Summarise ESS
  # Compute the average mean ess for selected algorithm runs 
  ess_mean <- cl_lapply_iteration_file(
    iteration[sensitivity == "best", ], 
    .file = "file_coord", 
    .fun = function(.sim, .input) {
      mean(.input$diagnostics$ess, na.rm = TRUE)
    }, .cl = 10L, .combine = unlist) 
  mean(ess_mean)
  
  #### Summarise the area spanned by 95 % of the distribution
  # This provides information on how well an animal has been localised 
  areas_mean <- cl_lapply_iteration_file(
    iteration[sensitivity == "best", ], 
    .file = "file_coord", 
    .fun = function(.sim, .input) {
      # Compute average area spanned by 95 % probability mass via 2d hist (~6 s)
      # (averaged over all time steps)
      particle_hr(.map = map, .coord = .input$states)
    }, .cl = 10L, .combine = unlist)
  # Compute average mean area (km2) across iterations
  mean(areas_mean / 1e6)
  
}


#### End of code. 
###########################
###########################