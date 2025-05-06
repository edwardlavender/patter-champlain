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
map <- terra::rast(here_input("map.tif"))


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
sim <- iteration[2, ]
# > unit_id  individual_id    time_id 
# > 16       24321         2016-01-01

#### (optional) Tweak selected parameters
sim[, phi := 1.3]

#### Define filter args 
args <- constructor_ac_real(.sim = sim,
                            .datasets = list(), 
                            .verbose = TRUE)
args_fwd <- args$forward
args_fwd$.n_particle <- 2.5e4L
args_fwd$.n_record   <- 1000L
args_fwd$.batch      <- NULL
args_fwd$.progress   <- julia_progress(enabled = TRUE)
# Collect moorings
moorings <- 
  args_fwd$.yobs$ModelObsAcousticLogisTrunc |> 
  lazy_dt() |>
  group_by(sensor_id) |> 
  slice(1L) |> 
  select(receiver_id = sensor_id, receiver_x, receiver_y, receiver_alpha, receiver_beta, receiver_gamma) |> 
  as.data.table()

#### Run filter
fwd <- do.call(pf_filter, args_fwd)

# iteration[2, ]
# > Weights from filter (1 -> 18762) are zero at time 12806:returning outputs from 1:12806. Note that all (log) weights at 12806 are -Inf.
# > @ 2014-11-22 17:20:00 

#### Record inputs/outputs
tnow <- as.numeric(Sys.time())
qs::qsave(sim, here_debug(tnow, "sim.qs"))
qs::qsave(moorings, here_debug(tnow, "moorings.qs"))
qs::qsave(args_fwd, here_debug(tnow, "input.qs"))
qs::qsave(fwd, here_debug(tnow, "output.qs"))


###########################
###########################
#### Animation

# Create animation (not in Julia session on Linux)
ani(.sim      = sim,
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


###########################
###########################
#### Interactive debugging

# Identify problematic time step/stamps
timeline <- fwd_args$timeline
tdt      <- data.table(timestep = 1:length(timeline), timestamp = timeline)
tstep    <- 12804
tstamp   <- timeline[tstep]

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