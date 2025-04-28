###########################
###########################
#### setup-data-detection.R

#### Aims
# 1) Sets up detection data for analysis with patter

#### Prerequisites
# 1) Raw data provided by M. Futia 


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Set global options
Sys.setenv("JULIA_SESSION" = FALSE)
set.seed(1L)

#### Load essential packages
library(proj.verse)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(lubridate)
library(tictoc)
files_source_r(here_src())

#### Load data 
epsg_utm   <- qs::qread(here_input("epsg-utm.qs"))
map        <- terra::rast(here_input("map.tif"))
map_bbox   <- qs::qread(here_input("map-bbox.qs"))
moorings   <- readRDS(here_data_raw_mf("OriginalReceiverSummary_2013-2017.rds"))
detections <- readRDS(here_data_raw_mf("lkt_detections_2013-2017.rds"))
surgery    <- fread(here_data_raw("mfutia", "model_comparison", "surgery_log.csv"))


###########################
###########################
#### Identify fish 

#### Define fish (id, size, tagging location)
# Define fish 
fish <- 
  detections |> 
  group_by(animal_id) |> 
  summarise(individual_id = animal_id[1], 
            len = length[1] / 1000, 
            lat = deploy_lat[1], 
            lon = deploy_long[1],
            # Checks 
            nlen = n_distinct(length), 
            nlat = n_distinct(lat), 
            nlon = n_distinct(lon)
            ) |> 
  as.data.table()
# Define tagging locations (UTM)
# * This code requires internet
xy <- 
  cbind(fish$lon, fish$lat) |> 
  terra::vect(crs = "EPSG:4326") |> 
  terra::project(epsg_utm) |> 
  terra::crds()
stopifnot(nrow(xy) > 0L)
fish[, x := xy[, 1]]
fish[, y := xy[, 2]]
# Check tagging locations
stopifnot(all(!is.na(terra::extract(map, xy)[, 1])))
terra::plot(map)
points(xy)
# Check tagging dates
# * Note that fish were tagged at different times
# * If we analyse the data in blocks e.g., months, we need to account for this
# * TO DO Confirm date format %m/%d/%Y
# * We could check that each fish is only associated with detections after tagging
range(detections$detection_timestamp_utc)
range(as.Date(surgery$cap_date, format = "%m/%d/%Y"))

#### Checks
# Each individual is associated with one length (presumably length @ tagging)
# deploy_lat and deploy_long seem to refer to fish tagging locations
stopifnot(all(fish$nlen == 1L))
stopifnot(all(fish$nlat == 1L))
stopifnot(all(fish$nlon == 1L))

#### Clean up
fish |> 
  select(individual_id, len, x, y, lon, lat) |> 
  as.data.table()

#### Comments
# length is total length (mm)
# Hansen et al. 2022 provide an equation to convert TL to FL for lake trout 
# (FL = 0.9143*TL-8.2772)


###########################
###########################
#### Prepare moorings

#### Define receiver coordinates (UTM)
# Note that re-projection requires an internet connection! 
head(moorings)
rxy <- 
  cbind(moorings$deploy_lon, moorings$deploy_lat) |> 
  terra::vect(crs = "EPSG:4326") |> 
  terra::project("EPSG:3175") |>
  terra::crds()
stopifnot(nrow(rxy) > 0L)
stopifnot(all(!is.na(terra::extract(map, rxy)$map_value)))

#### Receiver depths 
utils.add::basic_stats(moorings$depth, na.rm = TRUE)
# min  mean median  max   sd   IQR  MAD
# 1 3.4 13.56  11.85 45.7 9.06 11.85 8.82

#### Process moorings
moorings <- 
  moorings |> 
  mutate(receiver_station = StationName, 
         receiver_id = row_number(),
         receiver_sn = as.integer(as.character(receiver_sn)),
         receiver_start = as.POSIXct(paste0(deploy_date_time, "00:00:00"), tz = "UTC"), 
         receiver_end = as.POSIXct(paste0(recover_date_time, "00:00:00"), tz = "UTC"), 
         receiver_int = lubridate::interval(receiver_start, receiver_end),
         receiver_x = rxy[, 1],
         receiver_y = rxy[, 2]) |> 
  select(receiver_station, 
         receiver_id, receiver_sn, receiver_start, receiver_end, 
         receiver_int, receiver_x, receiver_y) |>
  as.data.frame()

#### Check deployment periods
ggplot(moorings) +
  geom_segment(aes(
    x    = receiver_start,
    xend = receiver_end,
    y    = factor(receiver_id),
    yend = factor(receiver_id)
  ), size = 2)

#### Add detection probability parameters
# This is implemented later


###########################
###########################
#### Prepare detections

#### Examine selected columns
nrow(detections)
head(detections)
table(detections$passFilter)
hist(detections$length)
max(detections$length)

#### Process detections
detections <- 
  detections |> 
  mutate(individual_id = as.integer(as.character(animal_id)), 
         timestamp = as.POSIXct(detection_timestamp_utc, tz = "UTC"),
         receiver_id = NA_integer_, 
         receiver_sn = as.integer(as.character(receiver_sn))) |>
  select(individual_id, 
         timestamp,
         receiver_sn) |> 
  as.data.table()

#### Order detections by duration
# This improves speed during algorithm testing 
# (NB: the code below works because duration is unique to each individual)
detections <- 
  detections |>
  group_by(individual_id) |> 
  mutate(duration = as.numeric(difftime(max(timestamp), min(timestamp), units = "days"))) |> 
  arrange(duration, timestamp) |> 
  mutate(-duration) |>
  as.data.table()

#### Define receiver_id (~4 s)
# Match using receiver_sn and the time stamps
tic()
for (i in 1:nrow(moorings)) {
  detections[receiver_sn == moorings$receiver_sn[i] & 
               timestamp %within% moorings$receiver_int[i], receiver_id := moorings$receiver_id[i]]
}
toc()
table(is.na(detections$receiver_id))

#### Drop 'extra' detections
detections <- detections[!is.na(receiver_id), ]


###########################
###########################
#### Clean up

# Define study period
study_start <- min(detections$timestamp)
study_end   <- max(detections$timestamp)
study_int   <- lubridate::interval(study_start, study_end)

#### Clean up moorings 
head(moorings)
nrow(moorings)
moorings <- 
  moorings |> 
  as.data.frame() |>
  mutate(int = lubridate::interval(receiver_start, receiver_end)) |> 
  filter(int_overlaps(int, study_int)) |> 
  select(receiver_station, receiver_id, receiver_start, receiver_end, receiver_x, receiver_y) |> 
  as.data.table()
nrow(moorings)
# Define moorings for simulation analyses
# * We average the locations of the receivers in each Station
# * (Receivers were redeployed in the same area (station) after servicing)
# * The receiver_start and receiver_end columns will be replaced later
#   in line with the simulation timeline (see sim-data.R)
moorings_sim <- 
  moorings |> 
  group_by(receiver_station) |> 
  mutate(receiver_x = mean(receiver_x), 
         receiver_y = mean(receiver_y)) |> 
  slice(1L) |> 
  mutate(receiver_id = row_number()) |> 
  select(-receiver_station) |> 
  as.data.table()
nrow(moorings_sim)
# Define moorings for real-world analyses
moorings_real <- 
  moorings |> 
  select(-receiver_station) |> 
  as.data.table()
rm(moorings)
# Validate positions
terra::plot(map)
points(moorings_sim$receiver_x, moorings_sim$receiver_y, pch = ".")
terra::plot(map)
points(moorings_real$receiver_x, moorings_real$receiver_y, pch = ".")

#### Clean up detections
head(detections)
detections <-
  detections |> 
  select(individual_id, timestamp, receiver_id) |> 
  as.data.table()


###########################
###########################
#### Write outputs

qs::qsave(fish, here_input("fish.qs"))
qs::qsave(moorings_sim, here_input_sim("moorings-xy.qs"))
qs::qsave(moorings_real, here_input_real("moorings.qs"))
qs::qsave(detections, here_input_real("detections.qs"))


#### End of code. 
###########################
###########################