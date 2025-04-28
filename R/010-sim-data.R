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
moorings        <- qs::qread(here_input_sim("moorings-xy.qs"))
pars_model_move <- qs::qread(here_input("pars-model-move-best.qs"))
pars_model_obs  <- qs::qread(here_input("pars-model-obs-best.qs"))
unitsets        <- qs::qread(here_input_sim("unitsets.qs"))


###########################
###########################
#### Simulate data 

#### Setup Julia
set_seed()
set_map(map)

#### Define n_sim
n_sim <- nrow(unitsets)

#### Define timeline (one-month)
timeline <- seq(as.POSIXct("2025-01-01 00:00:00", tz = "UTC"), 
                as.POSIXct("2025-01-31 23:58:00", tz = "UTC"), 
                by = "2 mins")

#### Define movement model
model_move <- model_move_trout(pars_model_move)
plot(model_move)

#### Define observation model
# * We consider the average receiver positions in each StationName
#   following Futia et al. (2024). 
# * This is a better representation of the study design than 
#   assuming all 153 rows in moorings represent available receivers.
# * We update receiver_start and receiver_end for the simulation timeline. 
stopifnot(nrow(moorings) == 31L)
moorings[, receiver_start := min(timeline) -  24 * 60 * 60]
moorings[, receiver_end := max(timeline) + 24 * 60 * 60]
model_obs <- model_obs_champlain(moorings, pars_model_obs)
plot(model_obs)

#### Define tagging locations
xinit <- fish[sample.int(n_sim, replace = TRUE), ]
xinit[, map_value := terra::extract(map, cbind(x, y))]
xinit <- xinit[, .(map_value, x, y)]
xinit <- model_move_xinit(.xinit = xinit, .n_particle = NULL)

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

#### Collate capture/recapture locations for each unit_id
# Select capture/recapture locations
xinits <- 
  paths |> 
  group_by(path_id) |> 
  slice(c(1, n())) |> 
  select(path_id, timestep, map_value, x, y) |> 
  as.data.table()
# Convert to list
xinits <- split(xinits, xinits$path_id)
lobstr::obj_size(xinits)

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
( ndet <- sapply(detections_by_path, nrow) )
stopifnot(all(ndet > 0))
# C) Validation
# (i) We should only record detections within receiver_gamma of receiver
cl_lapply(seq_len(n_sim), function(i) {
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
# (ii) Check the number of detections per simulation
sapply(detections_by_path, nrow) |> sort()


###########################
###########################
#### Filter simulations

#### Check which simulated datasets meet criteria for modelling
sapply(detections_by_path, function(d) {
  d <- detections_by_path[[1]]
  d |>
  summarise(
    duration = as.numeric(difftime(max(timestamp), min(timestamp), units = "days")), 
    ndays = length(unique(lubridate::floor_date(timestamp, "days"))),
    pdays = ndays / duration
  ) |> 
  pull(pdays) > 0.75
})


#### TO DO Focus on time series that pass criteria for modelling
# TO DO
# * Define criteria & implement as function
# * Implement function in both sim-data.R (here) and prepare-real.R


###########################
###########################
#### Write datasets to file

qs::qsave(timeline, here_input_sim("timeline.qs"))
qs::qsave(xinits, here_input_sim("xinits.qs"))
qs::qsave(paths, here_input_sim("paths.qs"))
qs::qsave(moorings, here_input_sim("moorings.qs"))
qs::qsave(acoustics_by_path, here_input_sim("acoustics-by-path.qs"))
qs::qsave(detections_by_path, here_input_sim("detections-by-path.qs"))
for (i in 1:nrow(unitsets)) {
  detections <- detections_by_path[[i]]
  qs::qsave(detections, unitsets$file_detection[i])
}


#### End of code. 
###########################
###########################