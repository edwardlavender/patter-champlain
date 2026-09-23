###########################
###########################
#### develop-data-validation.R

#### Aims
# 1) This script supports in-field validation of patter outputs

#### Prerequisites
# 1) setup-data-detection-pr-futia.R


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Set global options
Sys.setenv("JULIA_SESSION" = FALSE)
set.seed(123L)

#### Load essential packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(leaflet)
library(proj.verse)
files_source_r()

#### Load data
map <- terra::rast(here_input("map.tif"))
# Futia & Marsden (2025) dataset
futia_detections <- 
  qs::qread(here_data("supp", "model-obs", "futia-raw-validation.qs"))
futia_moorings <- 
  qs::qread(here_data_raw("model-obs","futia-et-al-2025",
                          "range_test_metadata_2021-2022.qs"))
# Pinheiro (unpublished dataset)
pinheiro_detections <-
  qs2::qs_read(
    here_data_raw("mfutia", 
                  "validation", 
                  "Champlain_RangeTestDetections_Summer2016_mhf.qs2"))


###########################
###########################
#### Process datasets

#### Process Futia & Marsden (2025) moorings (non-independent)
futia_moorings <- 
  futia_moorings |> 
  mutate(dataset = "F") |> 
  rename(receiver_lon = rec_lon, 
         receiver_lat = rec_lat) |> 
  distinct(dataset, receiver_sn, receiver_lon, receiver_lat) |> 
  select("dataset", "receiver_sn", "receiver_lon", "receiver_lat") |> 
  as.data.table()

n_futia_moorings <- nrow(futia_moorings)

#### Process associated detections 
futia_detections <-
  futia_detections |> 
  mutate(dataset = "F") |> 
  rename(receiver_sn = receiver_id) |> 
  filter(receiver_sn %in% futia_moorings$receiver_sn) |> 
  select("dataset", "transmitter_id", "tag_lon", "tag_lat", "start", "end",
         "timestamp", "receiver_sn",  "receiver_lon", "receiver_lat"
         ) |> 
  as.data.table()

#### Update moorings
# Receiver deployment periods are unrecorded, by some receivers were re-deployed 
table(futia_moorings$receiver_sn)
# We assume receivers were active only over time period of detections in each location
futia_moorings <- 
  futia_moorings |> 
  left_join(futia_detections |> 
              group_by(receiver_sn, receiver_lon, receiver_lat) |>
              mutate(receiver_start = min(timestamp), 
                     receiver_end = max(timestamp)) |>
              slice(1L) |> 
              select("receiver_sn", "receiver_lon", "receiver_lat", 
                     "receiver_start", "receiver_end") |> 
              ungroup() |> 
              as.data.table()) |>  
  select("dataset", 
         "receiver_sn",
         "receiver_start", "receiver_end", 
         "receiver_lon", "receiver_lat") |>
  filter(!is.na(receiver_start)) |>
  as.data.table()

#### Process Pinheiro (unpublished) moorings
# For this dataset, surrounding receivers are unknown so we focus on the (one)
# receiver that recorded detections to define moorings. 
pinheiro_moorings <- 
  pinheiro_detections |>  
  janitor::clean_names() |> 
  mutate(dataset = "P") |> 
  rename(receiver_lon = rec_deploy_lon, 
         receiver_lat = rec_deploy_lat) |> 
  distinct(dataset, receiver_sn, receiver_lon, receiver_lat) |> 
  mutate(receiver_start = min(pinheiro_detections$detection_timestamp_utc),
         receiver_end = max(pinheiro_detections$detection_timestamp_utc)) |> 
  select("dataset", 
         "receiver_sn", 
         "receiver_lon", "receiver_lat", 
         "receiver_start", "receiver_end") |> 
  as.data.table()

#### Process Pinheiro (unpublished) detections
pinheiro_detections <- 
  pinheiro_detections |>  
  janitor::clean_names() |>
  mutate(dataset = "P") |> 
  rename(transmitter_id = transmitter, 
         timestamp = detection_timestamp_utc, 
         tag_lon = tag_deploy_lon, 
         tag_lat = tag_deploy_lat, 
         receiver_lon = rec_deploy_lon, 
         receiver_lat = rec_deploy_lat) |> 
  group_by(transmitter_id, tag_lon, tag_lat) |>
  # Assign start/end times based on time of range test(s) using detections
  # This is necessary because unlike the Futia dataset test times are unknown
  mutate(start = min(timestamp), 
         end = max(timestamp)) |> 
  select("dataset", "transmitter_id", "tag_lon", "tag_lat", "start", "end",
         "timestamp", "receiver_sn",  "receiver_lon", "receiver_lat") |> 
  as.data.table()

#### Collate moorings

# Filter Futia & Marsden moorings to receivers within 
# (Note that all futia moorings are within 8000 m of at least one range testing tag)
rll <- 
  cbind(futia_moorings$receiver_lon, futia_moorings$receiver_lat)
tll <- 
  futia_detections |> 
  distinct(tag_lon, tag_lat) |> 
  as.matrix()
terra::distance(rll, tll, lonlat = TRUE) |> apply(1, function(x) any(x < 8000))

# Collate moorings
moorings <- 
  rbind(futia_moorings, pinheiro_moorings) |> 
  mutate(receiver_id = row_number(), .before = 1) |>
  as.data.table()

# Define receiver_x and receiver_y coordinates on map
rxy <- 
  cbind(moorings$receiver_lon, moorings$receiver_lat) |>
  terra::vect(crs = "WGS84") |>
  terra::project(terra::crs(map)) |> 
  terra::geom(df = TRUE)
moorings[, receiver_x := rxy$x]
moorings[, receiver_y := rxy$y]

#### Collate detections
detections <- 
  rbind(futia_detections, pinheiro_detections) |> 
  # Define receiver_id in detections
  left_join(moorings |> 
              select("dataset", 
                     "receiver_id", "receiver_sn", 
                     "receiver_lon", "receiver_lat") |> 
              as.data.table(), 
            by = c("dataset", "receiver_sn", "receiver_lon", "receiver_lat")) |> 
  filter(!is.na(receiver_id)) |>
  # Define an 'individual_id' column used to distinguish separate range tests
  arrange(dataset, transmitter_id, tag_lon, tag_lat, start, timestamp, receiver_id) |> 
  group_by(dataset, transmitter_id, tag_lon, tag_lat, start) |> 
  mutate(individual_id = cur_group_id(), .before = 1L) |>
  as.data.table()

#### Collate 'tagging' (range test) information
# This defines test IDs, the location of the range testing tagt, and the start/end time
tests <- 
  detections |> 
  distinct(individual_id, dataset, tag_lon, tag_lat, start, end) |>
  as.data.table()
         
#### Clean up
moorings |> 
  select("receiver_id", 
         "receiver_x", "receiver_y",
         "receiver_start", "receiver_end") |> 
  arrange(receiver_id) |> 
  as.data.table()

detections |> 
  select("individual_id", "timestamp", "receiver_id") |> 
  arrange(individual_id, timestamp, receiver_id) |> 
  as.data.table()

head(moorings)
head(moorings)

#### End of code. 
###########################
###########################
