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

#### Load essential packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(JuliaCall)
library(patter)
library(proj.verse)
library(tictoc)
files_source_r(here_src())

#### Load data
fish            <- qs::qread(here_input("fish.qs"))
moorings        <- qs::qread(here_input_sim("main", "moorings-xy.qs"))
pars_model_move <- qs::qread(here_input("pars-model-move-best.qs"))
pars_model_obs  <- qs::qread(here_input("pars-model-obs-best.qs"))


###########################
###########################
#### Simulate data 

#### Setup Julia
julia_connect()
set_seed()
julia_source(file.path("Julia", "src", "observation-model.jl"))
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

#### Simulate trajectories and observations (~179 s)
# This is implemented iteratively to generate exactly n_sim time series for modelling
# I.e., that pass filter_detections() criteria 
tic()
count              <- 1L
total              <- 1L
paths_by_path      <- list()
acoustics_by_path  <- list()
detections_by_path <- list()
while (count <= n_sim & total < 100L) {
  
  cli::cat_rule()
  print(paste(count, ", ", total))
  
  #### Define tagging locations
  # (Extract map_value via Patter for linux handling)
  xinit <- fish[sample.int(1L, replace = TRUE), ]
  julia_assign("x0", xinit$x)
  julia_assign("y0", xinit$y)
  xinit[, map_value := julia_eval('[Patter.extract(env, x0[i], y0[i]) for i in eachindex(x0)]')]
  xinit <- model_move_xinit(.xinit = xinit, .n_particle = NULL)
  stopifnot(all(xinit$map_value == 1L))
  
  #### Simulate movement paths (~1.6 s)
  # This returns a data.table with trajectories
  tic()
  path <- sim_path_walk(.timeline   = timeline, 
                        .state      = state_trout(), 
                        .model_move = model_move,
                        .xinit      = xinit,
                        .n_path     = 1L)
  path[, path_id := count]
  toc()
  # Validate that each path starts with the simulated xinit
  stopifnot(dplyr::all_equal(
    xinit[, .(x, y, heading)],
    path |> 
      group_by(path_id) |> 
      slice(1L) |> 
      ungroup() |>
      select("x", "y", "heading") |> 
      as.data.table()
  ))
  
  #### Collate capture/recapture locations
  # xinit <- 
  #   path |> 
  #   group_by(path_id) |> 
  #   slice(c(1, n())) |> 
  #   select(path_id, timestep, map_value, x, y) |> 
  #   as.data.table()
  
  #### Simulate acoustic observations (~5.5 s)
  tic()
  acoustics <- sim_observations(.timeline = timeline, 
                                .model_obs = model_obs)
  acoustics <- acoustics$ModelObsAcousticLogisTruncLos[[1]]
  toc()
  
  #### Collate detections data.table (to match real-world data structure)
  detections <- 
    acoustics |> 
    filter(obs == 1L) |> 
    mutate(individual_id = count, time_id = lubridate::floor_date(timestamp, "months")) |> 
    select(individual_id, time_id, timestamp, receiver_id = sensor_id, 
           receiver_x, receiver_y, receiver_alpha, receiver_beta, receiver_gamma) |> 
    as.data.table()
  
  #### Filter detections
  n <- nrow(filter_detections(detections))
  
  #### Record outputs, if successful
  if (n > 0L) {
    paths_by_path[[count]]      <- copy(path)
    acoustics_by_path[[count]]  <- copy(acoustics)
    detections_by_path[[count]] <- copy(detections)
    count <- count + 1L
  }
  total <- total + 1L
  
}
toc()

#### Collate datasets
# Check lists
stopifnot(length(paths_by_path) == n_sim)
stopifnot(length(acoustics_by_path) == n_sim)
stopifnot(length(detections_by_path) == n_sim)
# Collate data
paths      <- rbindlist(paths_by_path)
detections <- rbindlist(detections_by_path)
# Additional checks
stopifnot(identical(unique(paths$path_id), 1:n_sim))
stopifnot(identical(unique(detections$individual_id), 1:n_sim))

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

qs::qsave(timeline, here_input_sim("main", "timeline.qs"))
qs::qsave(paths, here_input_sim("main", "paths.qs"))
qs::qsave(moorings, here_input_sim("main", "moorings.qs"))
qs::qsave(acoustics_by_path, here_input_sim("main", "acoustics-by-path.qs"))
qs::qsave(detections, here_input_sim("main", "detections.qs"))


#### End of code. 
###########################
###########################