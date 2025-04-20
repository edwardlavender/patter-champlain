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
library(lubridate)
library(tictoc)
files_source_r(here_src())

#### Load data 
epsg_utm   <- qs::qread(here_input("epsg_utm.qs"))
map        <- terra::rast(here_input("map.tif"))
map_bbox   <- qs::qread(here_input("map-bbox.qs"))
moorings   <- readRDS(here_data_raw_mf("OriginalReceiverSummary_2013-2017.rds"))
detections <- readRDS(here_data_raw_mf("lkt_detections_2013-2017.rds"))


###########################
###########################
#### Identify fish 

#### Define fish & lengths 
fish <- 
  detections |> 
  group_by(animal_id) |> 
  summarise(len = length[1] / 1000, 
            nlen = n_distinct(length)) |> 
  as.data.table()

#### Comments
# length is total length (mm)
# Hansen et al. 2022 provide an equation to convert TL to FL for lake trout 
# (FL = 0.9143*TL-8.2772)


###########################
###########################
#### Prepare moorings

#### Define receiver coordinates (UTM)
head(moorings)
rxy <- 
  cbind(moorings$deploy_lon, moorings$deploy_lat) |> 
  terra::vect(crs = "EPSG:4326") |> 
  terra::project("EPSG:3175") |>
  terra::crds()
stopifnot(all(!is.na(terra::extract(map, rxy)$map_value)))

#### Process moorings
moorings <- 
  moorings |> 
  mutate(receiver_id = row_number(),
         receiver_sn = as.integer(as.character(receiver_sn)),
         receiver_start = as.POSIXct(paste0(deploy_date_time, "00:00:00"), tz = "UTC"), 
         receiver_end = as.POSIXct(paste0(recover_date_time, "00:00:00"), tz = "UTC"), 
         receiver_int = lubridate::interval(receiver_start, receiver_end),
         receiver_x = rxy[, 1],
         receiver_y = rxy[, 2]) |> 
  select(receiver_id, receiver_sn, receiver_start, receiver_end, receiver_int, receiver_x, receiver_y) |>
  as.data.frame()

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

#### Clean up moorings 
head(moorings)
moorings <- 
  moorings |> 
  select(receiver_id, receiver_start, receiver_end, receiver_x, receiver_y) |> 
  as.data.table()

#### Clean up detections
head(detections)
detections <-
  detections |> 
  select(individual_id, timestamp, receiver_id) |> 
  as.data.table()


###########################
###########################
#### Write outputs

qs::qsave(fish, here_input_real("fish.qs"))
qs::qsave(moorings, here_input("moorings.qs"))
qs::qsave(detections, here_input_real("detections.qs"))


#### End of code. 
###########################
###########################