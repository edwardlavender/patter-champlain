###########################
###########################
#### sim-data.R

#### Aims
# 1) Simulates trajectories and observations

#### Prerequisites
# 1) Develop movement and observation models


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
library(JuliaCall)
library(patter)
library(proj.verse)
library(tictoc)
files_source_r(here_src())
expect_no_geospatial()

#### Load data
fish            <- qs::qread(here_input("fish.qs"))
moorings        <- qs::qread(here_input_sim("moorings-xy.qs"))
pars_model_move <- qs::qread(here_input("pars-model-move-best.qs"))
pars_model_obs  <- qs::qread(here_input("pars-model-obs-best.qs"))


###########################
###########################
#### Simulate data 

#### Setup Julia
set_seed()
set_map(here_input("map.tif"))

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
# (Extract map_value via Patter for linux handling)
xinit <- fish[sample.int(n_sim, replace = TRUE), ]
julia_assign("x0", xinit$x)
julia_assign("y0", xinit$y)
xinit[, map_value := julia_eval('[Patter.extract(env, x0[i], y0[i]) for i in eachindex(x0)]')]
xinit <- model_move_xinit(.xinit = xinit, .n_particle = NULL)
stopifnot(all(xinit$map_value == 1L))

#### Simulate movement paths (~7 s)
# This returns a data.table with trajectories
tic()
paths <- sim_path_walk(.timeline   = timeline, 
                       .state      = state_trout(), 
                       .model_move = model_move,
                       .xinit      = xinit,
                       .n_path     = n_sim)
toc()
# (optional) Visualise moorings on plot, if .map specified
# points(model_obs$ModelObsAcousticLogisTrunc$receiver_x, 
#        model_obs$ModelObsAcousticLogisTrunc$receiver_y)
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

#### Simulate acoustic observations for each path
# ETA:
# * ~2.6 mins (SIA-LAVENDED-M)
# * ~6.5 mins (siam-linux20)
# * TO DO: Improve speed of Patter.jl.sim_observations()
tic()
acoustics_by_path <- sim_observations(.timeline = timeline, 
                                      .model_obs = model_obs)
acoustics_by_path <- acoustics_by_path$ModelObsAcousticLogisTrunc
toc()

#### Collate detections data.table (to match real-world data structure)
detections <- lapply(seq_len(n_sim), function(i) {
  acoustics_by_path[[i]] |> 
    filter(obs == 1L) |> 
    mutate(individual_id = i, time_id = lubridate::floor_date(timestamp, "months")) |> 
    select(individual_id, time_id, timestamp, receiver_id = sensor_id, 
           receiver_x, receiver_y, receiver_alpha, receiver_beta, receiver_gamma) |> 
    as.data.table()
}) |> rbindlist()

#### Validation
# We should only record detections within receiver_gamma of receiver
positions <- 
  paths |> 
  right_join(detections, by = c("path_id" = "individual_id", "timestamp")) |> 
  mutate(dist = patter:::dist_2d(cbind(x, y), 
                                 cbind(receiver_x, receiver_y),
                                 pairwise = TRUE)) |> 
  as.data.table()
stopifnot(all(positions$dist <= positions$receiver_gamma))


###########################
###########################
#### Write datasets to file

qs::qsave(timeline, here_input_sim("timeline.qs"))
qs::qsave(xinits, here_input_sim("xinits.qs"))
qs::qsave(paths, here_input_sim("paths.qs"))
qs::qsave(moorings, here_input_sim("moorings.qs"))
qs::qsave(acoustics_by_path, here_input_sim("acoustics-by-path.qs"))
qs::qsave(detections, here_input_sim("detections.qs"))


#### End of code. 
###########################
###########################