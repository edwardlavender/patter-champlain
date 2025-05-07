###########################
###########################
#### debug-patter.R

#### Aims
# 1) This script develops the patter analysis

#### Prerequisites
# < This code is designed for MacOS/Windows >
# 1) Develop best model
# 2) Implement model for a subset of real-world datasets
# 3) Identify convergence failures
# 4) Use this script to dig into the causes of convergence failures
# 5) Revise best model, as designed in develop-model-move.R & develop-model-detection.R

#### Results
# * A weakly correlated random walk (phi = 1.8) works for ~73/100 test runs (analysis-patter.R):
# - (2e4 particles, default n_resample, default n_move)
# - The problem with this model appears to be movement between more distant receivers
# - Particle animations show particles spread out too slowly to reach necessary receivers 
# - In contrast, high correlations make it easier for particles to move between disparate receivers
# - High correlations can be expensive with a truncated movement model (if lots of moves)
# - But can work well with an un-truncated model and few particles (fast)
# * A more correlated random walk (phi = 0.3) works for ~93/100 individuals (analysis-patter.R):
# - (2e4 particles, default n_resample, default n_move)
# - For the 7 remaining individuals, the challenge lies in the south of the study area
# - Particles can die easily in the narrow channels and receiver gates
#   can block particles from where they need to get to. 
# - Strategies that help include:
#   - Less resampling
#   - More particles
#   - More moves (expensive)
# - The challenge remains unsolved for at least two individuals


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Set global options
patter::julia_connect()
if (JuliaCall::julia_eval("Threads.nthreads()") == 1L) {
  warning("This script should leverage multiple Julia threads.")
}

#### Load essential packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(glue)
library(JuliaCall)
library(patter)
library(patter.workflows)
library(proj.verse)
library(tictoc)
files_source_r(here_src())
expect_no_geospatial()

#### Load data
map <- terra::rast(here_input("map.tif")) |> readAll()


###########################
###########################
#### Setup

#### Define analysis type
analysis <- "real"

#### Define analysis-specific routines
here_input_analysis       <- switch_here_input_analysis(analysis)
here_output_analysis_main <- switch_here_output_analysis_main(analysis)
constructor_ac_analysis   <- switch_constructor_ac_analysis(analysis)

#### Define analysis-specific data
iteration <- qs::qread(here_input_analysis("iteration-patter.qs"))
iteration[, file_diag := file.path(folder_coord, "diagnostics.qs")]
iteration[, file_output := file_diag]
nrow(iteration)

#### Select iterations (by mobility)
iteration <- iteration[mobility == 216, ]
if (TRUE) {
  iteration <- iteration[sensitivity == "best", ]
  iteration <- iteration[1:min(c(.N, 100L)), ]
}

#### Set maps 
set_map(here_input("map.tif"))
set_vmap(.vmap = here_input("vmap", iteration$mobility[1], "vmap.tif"))


###########################
###########################
#### Run filter

#### Select individual
sim <- iteration[index == 617, ]
# > unit_id  individual_id    time_id    index 
# > 16       24321         2016-01-01
# > 28       24322         2014-11-01    64
# > 102      24325         2016-04-01    358
# > 151      24328         2015-10-01    400
# > 231      24331         2016-03-01    610 -> unsolved
# > 232      24331         2016-04-01    617 -> unsolved 

#### (optional) Tweak selected parameters
sim[, phi := 0.3]

#### Define baseline filter args 
# (We may further customise inputs below)
args <- constructor_ac_real(.sim = sim,
                            .datasets = list(), 
                            .verbose = TRUE)
args_fwd <- args$forward
args_fwd$.n_particle <- 2e4L
args_fwd$.n_record   <- 2000L
args_fwd$.batch      <- NULL
args_fwd$.progress   <- julia_progress(enabled = TRUE)
# Collect relevant arguments
timeline <- args_fwd$.timeline
moorings <- 
  args_fwd$.yobs$ModelObsAcousticLogisTrunc |> 
  lazy_dt() |>
  group_by(sensor_id) |> 
  slice(1L) |> 
  select(receiver_id = sensor_id, receiver_x, receiver_y, receiver_alpha, receiver_beta, receiver_gamma) |> 
  as.data.table()

#### (optional) Visualise movement model
# You need quite small turning angles to generate correlated looking paths
# Paths only start to 'explore' study area when phi < 1
if (FALSE) {
  sim_path_walk(.map = map, 
                .timeline = args_fwd$.timeline, 
                .state = args_fwd$.state, 
                .model_move = "ModelMoveCXY(env, 216, truncated(Gamma(3.25, 25), upper = 216), Normal(0.0, 0.5));", 
                .n_path = 4L, 
                .one_page = TRUE)
}

#### (optional) Visualise particles without observations
if (FALSE) {
  # Run filter without observations
  args_nodata             <- args_fwd
  args_nodata$.model_move <- "ModelMoveCXY(env, 216, truncated(Gamma(3.25, 25), upper = 216), Normal(0.0, 1.0))"
  args_nodata$.yobs       <- list()
  args_nodata$.xinit      <- data.table(x = 1778217, y = 947554.8, map_value = 1, heading = 0.3)
  args_nodata$.n_particle <- 5e3L
  args_nodata$.n_record   <- 1e3L
  args_nodata$.n_resample <- 1
  fwd_nodata              <- do.call(pf_filter, args_nodata)
  # Examine output
  animate_xyt(.map = map, 
              .coord = fwd_nodata$states, 
              .steps = seq(1, length(args_nodata$.timeline), by = 200),
              .folder = here_fig("model-move"), 
              .cl = 12L)
  # Results
  # > With phi = 1.0, particles spread slowly over landscape
  # > With phi = 0.5, particles spread more quickly
  # - The filter is slower
  # - It is still difficult to travel down narrow channels
  # - Larger time steps may help but cause difficulty with land
}

#### (optional) Tweak detection probability model
if (FALSE) {
  # Plot detection probability model 
  args_fwd$.yobs$ModelObsAcousticLogisTrunc |> 
    model_obs_acoustic_logis_trunc() |>
    plot()
  # Visualise more restrictive model 
  a <- 1.414281; b <- -0.002016435
  lines(1:8000, plogis(a + 1:8000 * b), col = "red")
  # (optional) Update .yobs with more restrictive model
  args_fwd$.yobs$ModelObsAcousticLogisTrunc[, receiver_alpha := a]
  args_fwd$.yobs$ModelObsAcousticLogisTrunc[, receiver_beta := -b]
}

#### (optional) Tweak resampling settings
# With resampling at each time step, populations of particles in some channels ultimately die out
# We'll try resampling near to acoustic detections or when ESS < 500
tdet <- args_fwd$.yobs$ModelObsAcousticLogisTrunc[obs == 1, .(timestamp)]
tdet[, timestep := (1:length(timeline))[match(timestamp, timeline)]]
tres <- sapply(tdet$timestep, function(t) (t - 5):t, USE.NAMES = FALSE) |> 
  unlist() |> sort() |> unique()
tres <- tres[tres > 0]
length(tres)

#### Run filter
# ~2 mins with 2e4L particles
# ~10 mins with 5e4L particles, phi = 0.3, .n_move = 5000
# ~16 mins with 2e4L particles, phi = 0.4, .n_move = 100_000
args_fwd$.n_resample <- 500
args_fwd$.t_resample <- tres
args_fwd$.n_move     <- 10000
args_fwd$.n_particle <- 5e4
args_fwd$.model_move <- "ModelMoveCXY(env, 216, truncated(Gamma(3.25, 25), upper = 216), Normal(0.0, 0.3))";
fwd <- do.call(pf_filter, args_fwd)

#### Results: on the causes of convergence failures
# (Results derived from filter runs + examination of particle animations (see below))
# 
# iteration[2, ], 2e4L particles:
# * phi = c(1.8, 1.3, 0.5) stuck around time step ~12806 
# * phi = 0.5 _can_ work but is somewhat 'forced' to by containers
# * phi = 0.4 works:
#   - With regular resampling this looks like 'luck' though
#   - (acoustic containers have a strong effect)
#   - less regular sampling seems to help particles to spread out
# * phi = 0.3, 5e3 particles, .n_move = 1 works, in 1 min
# * phi = 0.3-0.4 looks like the right setting for this individual
#
# iteration[index == 64, ]:
# * this individual gets stuck when it moves south into the narrow channel (@ time 5395)
# * the particles behave sensibly
# * but they 'die out' too soon in the narrow channel (before the detection)
# * It is very hard to move back into the channel given a receiver in the way, given:
#   - High movement correlation 
#   - Current detection probability model
# * Possible solutions:
#   - More particles 
#   - Truncated movement model
#   - Weaker movement correlation 
# * Evaluation: 
#   - Fails with 1e5 particles and .n_move = 1000
#   - With 5e4L particles and .n_move = 1000 or .n_move = 5000, we can get this working
#   - (But only just, looking at the animation!)
#   - Reducing correlation in turning angle (N(0, 1.0), N(0, 1.3), N(0, 1.8)) does not help
#   - (N(0, 1.3), N(0, 1.8) fail even earlier)
#   - This suggests it is just hard for particles to survive in some of the narrow channels
#
# iteration[index == 358, ]:
# * This seems to work here (failed on server)
#
# iteration[index == 400, ]:
# * The particles get stuck in the narrow channel between two receivers
# * With the current detection probability model, it is hard to escape
#   given receivers at both ends of the channel. 
# * This happens even with 1e5 particles & n_move = 1
# * With 5e4 particles & n_move = 5000L, the problem is solved
#
# iteration[index == 610, ]:
# * Similar issue
# * Die out in the narrow channel in the south (long gap)
# * Hard to move through receiver barrier back into channel 
# * It is hard to get this individual to converge
#   even with reduced resampling, 1e5 particles, 5000 moves
#
# iteration[index = 617, ]
# * Similar issue as index 610
# * Unsolved with 5e4 particles and 100,000 moves
# * A more restrictive detection probability model 
#   (which would reduce the barrier effect) did not help

#### Record problematic time step/stamps
timeline <- args_fwd$.timeline
tdt      <- data.table(timestep = 1:length(timeline), timestamp = timeline)
tstep    <- max(fwd$diagnostics$timestep)
tstamp   <- timeline[tstep]

#### Record inputs/outputs
if (TRUE) {
  tnow <- as.numeric(Sys.time())
  dir.create(here_debug(tnow))
  qs::qsave(sim, here_debug(tnow, "sim.qs"))
  qs::qsave(moorings, here_debug(tnow, "moorings.qs"))
  qs::qsave(args_fwd, here_debug(tnow, "input.qs"))
  qs::qsave(fwd, here_debug(tnow, "output.qs"))
} else {
  # Read latest outputs
  tnows <- gtools::mixedsort(list.files(here_debug()))
  tnow <- tnows[length(tnows)]
  fwd  <- qs::qread(here_debug(tnow, "output.qs"))
}


###########################
###########################
#### Animation

# Create animation
# * For 12,000 steps: This takes 6 min (12 cl) plus >2 min for 12,000 steps
# * NB: Running this for a few steps with .cl = 1L seems to suppress a segmentation
#   fault when it then run in parallel for a larger time series below.
#   If you jump to the parallel version, it can throw a segmentation fault
animate_ac(.sim      = sim,
           .map      = map,
           .steps    = 1:10L,
           .input    = args_fwd, 
           .output   = fwd,
           .tnow     = tnow, 
           .cl       = 1L)
animate_ac(.sim      = sim,
           .map      = map,
           .steps    = 10000:11000,
           .input    = args_fwd, 
           .output   = fwd,
           .tnow     = tnow, 
           .cl       = 10L)

# Map
map_pou(.map = map, .coord = fwd$states)
map_dens(.map = map, .coord = fwd$states, .discretise = TRUE)


###########################
###########################
#### Interactive debugging

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


#### End of code. 
###########################
###########################