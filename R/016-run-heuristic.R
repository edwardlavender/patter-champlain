###########################
###########################
#### run-heuristic.R

#### Aims
# 1) Run the heuristic (Int) algorithm for the detection dataset 
#    (We will compare residency estimates derived by this approach to patter)

#### Prerequisites
# 1) NA


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Load packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(proj.verse)
library(tictoc)
files_source_r()

#### Load data
analysis    <- "sim"
subanalysis <- "main"
detections  <- qs::qread(here_input(analysis, subanalysis, "detections.qs"))
regions     <- terra::rast(here_input("regions.tif"))


###########################
###########################
#### Run interpolation algorithm

#### Format detections data.table for glatos 
# Define receiver coordinates (lon, lat)
rll <- 
  cbind(detections$receiver_x, detections$receiver_y) |> 
  terra::vect(crs = terra::crs(regions)) |> 
  terra::project("WGS84") |> 
  terra::geom() |> 
  as.data.table()
# Define detections 
detections_glatos <- 
  detections |> 
  rename(animal_id               = individual_id,
         detection_timestamp_utc = timestamp) |> 
  mutate(deploy_long = rll$x, 
         deploy_lat  = rll$y) |> 
  select("animal_id", "detection_timestamp_utc", "deploy_long", "deploy_lat") |> 
  as.data.frame()

#### Define transition matrix
# Define map (lon, lat) 
map_ll <-
  regions |> 
  terra::clamp(lower = 1, upper = 1) |> 
  terra::app(as.numeric) |> 
  terra::project("WGS84")
terra::plot(map_ll)
# Make transition matrix (~6 s)
tic()
tr <- gdistance::transition(raster::raster(map_ll), 
                            transitionFunction = mean, 
                            directions = 8)
tr <- gdistance::geoCorrection(tr, type = "c")
toc()

#### Run interpolation (~193 s)
# Interpolate data with 60 min bin timestamp (3600 seconds)
set.seed(111)
tic()
positions <- 
  glatos::interpolate_path(detections_glatos,
                           trans          = tr,
                           int_time_stamp = 3600, 
                           lnl_thresh     = 0.999, 
                           show_progress  = TRUE)
toc()
# Confirm that all coordinates are valid on the map
stopifnot(all(!is.na(terra::extract(map_ll, cbind(positions$longitude, positions$latitude))[, 1])))

#### Average positions for timesteps that include multiple detection locations
# Define x, y locations for averaging
pxy <- 
  cbind(positions$longitude, positions$latitude) |>
  terra::vect(crs = "WGS84") |> 
  terra::project(terra::crs(regions)) |> 
  terra::geom()|> 
  as.data.table()
# Implement averaging 
positions_centroid <- 
  positions |>
  mutate(x = pxy$x, y = pxy$y) |> 
  group_by(animal_id, bin_timestamp) |> 
  reframe(x = mean(x),
          y = mean(y),
          record_type = unique(record_type)) |> 
  mutate(region = terra::extract(regions,
                                 terra::vect(cbind(x, y), crs = terra::crs(regions)),
                                 search_radius = 2000)$map_value)
# (optional) Recompute lon, lat coordinates
pll <- 
  cbind(positions_centroid$x, positions_centroid$y) |>
  terra::vect(crs = terra::crs(regions)) |>
  terra::project("WGS84") |> 
  terra::geom() |>
  as.data.frame()
positions_centroid$longitude <- pll$x
positions_centroid$latitude  <- pll$y

#### Validate coordinates
# When we move to UTM coordinates, we produce some NAs 
v <- is.na(terra::extract(regions, cbind(pxy$x, pxy$y))[, 1])
table(v)
v <- is.na(terra::extract(regions, cbind(positions_centroid$x, positions_centroid$y))[, 1])
table(v)
# There are some NAs even if we use average lon, lat coordinates
v <- is.na(terra::extract(map_ll, cbind(pll$x, pll$y))[, 1])
table(v)
terra::plot(regions)
points(positions_centroid$x[which(v)], positions_centroid$y[which(v)])
# Fix NAs
# * NAs are handled using terra::extract(, ... search_radius above)
# * (optional) TO DO Re-design glatos::interpolate_path() using UTM grid
unique(positions_centroid$region)
table(is.na(positions_centroid$region))
stopifnot(all(!is.na(positions_centroid$region)))

#### Compute residency for each region relative to total number of positions
# head(qs::qread(here_output_sim_main("synthesis", "residency-skill.qs")))
residency <- 
  positions_centroid |> 
  group_by(animal_id) |> 
  mutate(n_positions_total = n_distinct(bin_timestamp)) |> 
  ungroup() |> 
  group_by(animal_id, region) |> 
  summarise(n_positions_in_region = n_distinct(bin_timestamp),
            estimate = n_positions_in_region / n_positions_total[1]) |> 
  ungroup() |> 
  group_by(animal_id) |> 
  tidyr::complete(region, fill = list(estimate = 0)) |> 
  ungroup() |> 
  mutate(algorithm = "Int",
         sensitivity = "Int", 
         sensitivity_label = "Int") |> 
  select(individual_id = "animal_id", "algorithm", "sensitivity", "sensitivity_label",  "region", "estimate") |> 
  arrange(individual_id, algorithm, sensitivity, region) |> 
  as.data.table()

# Confirm every individual has one residency estimate for every region
stopifnot(all(table(residency$individual_id, residency$region) == 1))

# Record outputs
qs::qsave(residency, here_output(analysis, subanalysis, "synthesis", "residency-int.qs"))


#### End of code