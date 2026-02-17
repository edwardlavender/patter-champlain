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
library(ggplot2)
library(patter)
library(proj.verse)
library(tictoc)
files_source_r(here_src())

#### Load data
map             <- terra::rast(here_input("map.tif"))
champlain_utm   <- qs::qread(here_input("champlain-utm.qs")) |> sf::st_geometry()
pars_model_move <- qs::qread(here_input("pars-model-move-best.qs"))
pars_model_obs  <- qs::qread(here_input("pars-model-obs-best.qs"))
moorings        <- qs::qread(here_input_sim("main", "moorings-xy.qs"))
xinit           <- fread(here_data_raw_mf("champlain_patter_startLocations.csv"))

#### Set local variables
set.seed(123L)


###########################
###########################
#### Simulate data 

#### Define n_sim
n_sim <- 100L 

#### Define simulation settings
iter <- 
  data.table(file_timeline = here_input_sim("main", "timeline.feather"), 
             file_xinit    = here_input_sim("main", "xinit.feather"),
             file_moorings = here_input_sim("main", "moorings.feather"), 
             file_paths_raw = here_input_sim("main", "paths-raw.feather"), 
             file_acoustics_raw = here_input_sim("main", "acoustics-raw.feather")) |> 
  cbind(as.data.table(pars_model_move)) |> 
  cbind(as.data.table(pars_model_obs))
arrow::write_feather(iter, here_input_sim("main", "iter.feather"))

#### Define timeline (one-month)
timeline <- seq(as.POSIXct("2025-01-01 00:00:00", tz = "UTC"), 
                as.POSIXct("2025-01-31 23:58:00", tz = "UTC"), 
                by = "2 mins")
timeline <- data.table(timestamp = timeline)
arrow::write_feather(timeline,  iter$file_timeline)

#### Simulate initial locations
# (1) Sample starting points as in Futia et al. (2024)
# Define starting points
if (FALSE) {
  xinit <- 
    cbind(xinit$start_lon, xinit$start_lat) |> 
    terra::vect(crs = "WGS84") |> 
    terra::project(terra::crs(map)) |> 
    terra::geom() |> 
    as.data.table() |> 
    select("x", "y") 
  # Resample 
  xinit <- 
    xinit[sample.int(.N, size = n_sim, replace = TRUE)] |> 
    mutate(map_value = terra::extract(map, cbind(x, y))$map_value, 
           heading = runif(n(), 0, 2 * pi)) |> 
    select("map_value", "x", "y", "heading") |>
    as.data.table()
}
# (2) Sample starting points uniformaly across study area 
if (TRUE) {
  xinit <-
    map |>
    terra::spatSample(size = 100L, method = "random",
                      replace = FALSE, na.rm = TRUE,
                      xy = TRUE) |>
    mutate(heading = runif(n(), 0, 2 * pi)) |>
    select("map_value", "x", "y", "heading") |>
    as.data.table()
  stopifnot(nrow(xinit) == 100L)
}
stopifnot(!all(is.na(xinit$map_value)))
arrow::write_feather(xinit, iter$file_xinit)
terra::plot(map)
points(xinit$x, xinit$y)
# text(xinit$x, xinit$y, 1:nrow(xinit))

#### Define model_move
# The settings are defined above

#### Define model_obs_pars
moorings <-
  moorings |> 
  cbind(as.data.table(pars_model_obs)) |> 
  # mutate(receiver_start = lubridate::floor_date(min(timeline$timestamp)), 
  #        receiver_end = lubridate::floor_date(max(timeline$timestamp))) |>
  select(sensor_id = "receiver_id", 
         "receiver_x", "receiver_y", 
         "receiver_alpha", "receiver_beta", "receiver_gamma") |> 
  as.data.table()
arrow::write_feather(moorings, iter$file_moorings)

#### Simulate observations (~500 s)
# Run 001-sim-data.jl 
if (!file.exists(iter$file_acoustics_raw)) {
  tic()
  julia <- Sys.which("julia")
  stopifnot(nzchar(julia))
  system2("/Users/lavended/.juliaup/bin/julia",
          args = "./Julia/001-sim-data.jl",
          env = c(LD_LIBRARY_PATH = ""),
          stdout = TRUE,
          stderr = TRUE)
  toc()
}

#### Collate outputs
paths     <- arrow::read_feather(iter$file_paths_raw)
acoustics <- arrow::read_feather(iter$file_acoustics_raw)
moorings  <- 
  moorings |>
  mutate(receiver_start = lubridate::floor_date(min(timeline$timestamp)),
         receiver_end = lubridate::floor_date(max(timeline$timestamp))) |> 
  select("receiver_id" = "sensor_id", 
         "receiver_x", "receiver_y", 
         "receiver_start", "receiver_end",
         "receiver_alpha", "receiver_beta", "receiver_gamma") |> 
  as.data.table()

#### Process datasets
# Process paths
head(paths)
setDT(paths)
# Process acoustics
acoustics <- 
  acoustics |> 
  select("path_id", "timestamp", 
         "obs", "sensor_id", 
         "receiver_x", "receiver_y", 
         "receiver_alpha", "receiver_beta", "receiver_gamma") |> 
  as.data.table()
# Process detections data.table (to match real-world data structure)
detections <- 
  acoustics |> 
  filter(obs == 1L) |> 
  mutate(individual_id = path_id) |> 
  select("individual_id", "timestamp", receiver_id = "sensor_id", 
         "receiver_x", "receiver_y", "receiver_alpha", "receiver_beta", "receiver_gamma") |> 
  as.data.table()

#### Validate datasets
# Check path_ids 
stopifnot(identical(unique(paths$path_id), 1:n_sim))
stopifnot(identical(unique(acoustics$path_id), 1:n_sim))
# stopifnot(identical(unique(detections$individual_id), 1:n_sim))
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
# We should only record detections within receiver_gamma of receiver
positions <- 
  paths |> 
  right_join(detections, by = c("path_id" = "individual_id", "timestamp")) |> 
  mutate(dist = patter:::dist_2d(cbind(x, y), 
                                 cbind(receiver_x, receiver_y),
                                 pairwise = TRUE)) |> 
  as.data.table()
stopifnot(all(positions$dist <= positions$receiver_gamma))

#### Examine simulated datasets
# Examine paths 
head(paths)
# Plot paths (12 mins)
if (FALSE) {
  tic()
  png(here_fig("sim", "main", "qplot-paths.png"), 
      height = 20, width = 20, units = "in", res = 800)
  p <- 
    ggplot() +
    geom_sf(data = champlain_utm) + 
    geom_path(aes(x, y, colour = timestep), data = paths, linewidth = 0.2) + 
    facet_wrap(~path_id, nrow = 20, ncol = 5) 
  print(p)
  dev.off()
  toc()
}
# Review the number of detections per individual
acoustics |> 
  group_by(path_id) |> 
  summarise(n_detections = length(which(obs == 1L))) |> 
  ungroup() |>
  summarise(utils.add::basic_stats(n_detections))
# Review number of individuals with zero detections
acoustics |> 
  group_by(path_id) |> 
  summarise(zero_detections = all(obs == 0L)) |> 
  filter(zero_detections)
# Plot detections (~5 s)
tic()
png(here_fig("sim", "main", "qplot-detections.png"), 
    height = 20, width = 20, units = "in", res = 800)
p <- 
  acoustics[obs == 1L, ] |> 
  ggplot() + 
  geom_point(aes(timestamp, path_id)) + 
  theme_bw()
print(p)
dev.off()
toc()


###########################
###########################
#### Write datasets to file

qs::qsave(timeline, here_input_sim("main", "timeline.qs"))
qs::qsave(paths, here_input_sim("main", "paths.qs"))
qs::qsave(moorings, here_input_sim("main", "moorings.qs"))
qs::qsave(acoustics, here_input_sim("main", "acoustics.qs"))
qs::qsave(detections, here_input_sim("main", "detections.qs"))


#### End of code. 
###########################
###########################