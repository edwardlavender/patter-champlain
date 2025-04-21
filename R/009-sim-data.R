###########################
###########################
#### sim-data.R

#### Aims
# 1) Simulates trajectories and observations

#### Prerequisites
# 1) Develop movement and observation models
# 2) This code should be run on MacOS or Windows


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
library(patter)
library(proj.verse)
library(tictoc)
files_source_r(here_src())

#### Load data
map             <- terra::rast(here_input("map.tif"))
fish            <- qs::qread(here_input("fish.qs"))
moorings        <- qs::qread(here_input("moorings.qs"))
pars_model_move <- qs::qread(here_input("pars-model-move.qs"))
pars_model_obs  <- qs::qread(here_input("pars-model-obs.qs"))


###########################
###########################
#### Simulate data 

#### Setup Julia
set_seed()
set_map(map)

#### Define n_sim
n_sim <- 30L

#### Define timeline (one-month)
timeline <- seq(as.POSIXct("2025-01-01 00:00:00", tz = "UTC"), 
                as.POSIXct("2025-01-31 23:58:00", tz = "UTC"), 
                by = "2 mins")

#### Define movement model
model_move <- model_move_trout(pars_model_move)
plot(model_move)

#### Define observation model
# We assume all receivers were active over the simulated study period
# This is a 'best-case' scenario! 
model_obs <- model_obs_champlain(moorings, pars_model_obs)
plot(model_obs)

#### Define tagging locations
xinit <- model_move_xinit(.map = map, 
                          .xinit = fish[sample.int(n_sim, replace = TRUE), ])

#### Simulate movement paths (~7 s)
# This returns a data.table with trajectories
tic()
paths <- sim_path_walk(.map        = map, 
                       .timeline   = timeline, 
                       .state      = state_trout(), 
                       .model_move = model_move,
                       .xinit      = xinit,
                       .n_path     = n_sim)
toc()
# (optional) Visualise moorings on plot
points(model_obs$ModelObsAcousticLogisTrunc$receiver_x, 
       model_obs$ModelObsAcousticLogisTrunc$receiver_y)
# Validate that each path starts with the simulated xinit
stopifnot(dplyr::all_equal(
  xinit[, .(x, y, heading)],
  paths |> 
    group_by(path_id) |> 
    slice(1L) |> 
    ungroup() |>
    select("x", "y", "heading") |> 
    as.data.table()
))

#### Simulate observations for each path
# ETA: 30 s x n_sim = 15 mins!
# TO DO: Improve speed of Patter.jl.sim_observations()
# A) Simulate acoustic observations 
tic()
acoustics_by_path <- sim_observations(.timeline = timeline, 
                                      .model_obs = model_obs)
acoustics_by_path <- acoustics_by_path$ModelObsAcousticLogisTrunc
toc()
# B) Focus on detections
detections_by_path <- acoustics_by_path
for (i in seq_len(n_sim)) {
  detections <- detections_by_path[[i]]
  detections_by_path[[i]] <- detections[obs == 1L, ]
}
# C) Validation
# * We should only record detections within receiver_gamma of receiver
lapply(seq_len(n_sim), function(i) {
  # Join path and detection data.tables
  path       <- paths[path_id == i, ][, .(timestamp, x, y)]
  detections <- detections_by_path[[i]][, .(timestamp, receiver_x, receiver_y, receiver_gamma)]
  positions  <- right_join(path, detections, by = "timestamp") |> as.data.table()
  # Compute distances between individual and receiver @ which it was detected
  positions[, dist := terra::distance(cbind(x, y), 
                                      cbind(receiver_x, receiver_y),
                                      lonlat = FALSE,
                                      pairwise = TRUE)]
  # Verify all distances are less than the detection threshold
  stopifnot(all(positions$dist <= positions$receiver_gamma))
})

#### Write datasets to file
qs::qsave(paths, here_input_sim("paths.qs"))
qs::qsave(acoustics_by_path, here_input_sim("acoustics-by-path.qs"))
qs::qsave(detections_by_path, here_input_sim("detections-by-path.qs"))


#### End of code. 
###########################
###########################