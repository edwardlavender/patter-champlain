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
library(lubridate)
library(proj.verse)
files_source_r()

#### Load data
map <- terra::rast(here_input("map.tif"))
# Futia & Marsden (2025) dataset
futia_detections <- 
  qs::qread(here_data("supp", "model-obs", "futia-raw-validation.qs"))
futia_moorings <- 
  fread(here_data_raw("mfutia", "validation", "dissertation_receiver_log.csv"))
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
# The futia_moorings contains all receiver deployments
# We use this over the range_test_metadata_2021-2022.qs
# used in the setup-detection-pr-futia.R script b/c that is pre-processed
# data from MF that does not include receiver start/end times
# which we need for the in-field validation. 

futia_moorings <- 
  futia_moorings |> 
  mutate(
    dataset = "F", 
    # Define receiver_start and receiver_end
    # Use date_deploy/date_recover columns as deploy_date_time not provided for all receivers
    # Assume deployment dates defined for America/New York time zone
    receiver_start = 
      date_deploy |> 
      as.Date(format = "%m/%d/%Y") |> 
      as.POSIXct(tz = "America/New York") |> 
      lubridate::with_tz("UTC"),
    receiver_end = 
      date_recover |> 
      as.Date(format = "%m/%d/%Y") |> 
      as.POSIXct(tz = "America/New York") |> 
      lubridate::with_tz("UTC")) |> 
  filter(!is.na(receiver_start) & !is.na(receiver_end)) |>
  # Focus on receivers 
  filter(int_overlaps(interval(receiver_start, receiver_end), 
                      interval(min(futia_detections$detection_timestamp_utc), 
                               max(futia_detections$detection_timestamp_utc)))) |> 
  rename(receiver_lon = deploy_lon, 
         receiver_lat = deploy_lat) |> 
  select("dataset", "receiver_sn",
         "receiver_start", "receiver_end",
         "receiver_lon", "receiver_lat") |> 
  as.data.table()

# Filter Futia & Marsden moorings to receivers within ~10000 m of a tag 
# * This improves efficiency for range tests
rll <- 
  cbind(futia_moorings$receiver_lon, futia_moorings$receiver_lat)
tll <- 
  futia_detections |> 
  distinct(tag_lon, tag_lat) |> 
  as.matrix()
futia_moorings <-
  futia_moorings |>
  filter(terra::distance(rll, tll, lonlat = TRUE) |> 
           apply(1, function(x) any(x < 10000))) |> 
  as.data.table()

# Record number of moorings
n_futia_moorings <- nrow(futia_moorings)

# Note that none of the receiver deployment/retrieval dates overlap with range test dates (good)
unique(futia_moorings$receiver_start) %in%
  unique(as.Date(lubridate::with_tz(futia_detections$start, "UTC"))) |> any()
unique(futia_moorings$receiver_end) %in%
  unique(as.Date(lubridate::with_tz(futia_detections$start, "UTC"))) |> any()

#### Process associated detections 
# detection time stamps defined by detection_timestamp_utc
# test start/end times defined by start/end in America/New York and converted to UTC
futia_detections <-
  futia_detections |> 
  mutate(dataset = "F", 
         timestamp = detection_timestamp_utc, 
         start     = lubridate::with_tz(start, "UTC"), 
         end       = lubridate::with_tz(end, "UTC")) |> 
  rename(receiver_sn = receiver_id) |> 
  filter(receiver_sn %in% futia_moorings$receiver_sn) |> 
  select("dataset", "transmitter_id", "tag_lon", "tag_lat", "start", "end",
         "timestamp", "receiver_sn",  "receiver_lon", "receiver_lat"
  ) |> 
  as.data.table()

#### Process Pinheiro (unpublished) moorings
# For this dataset, surrounding receivers are unknown so we focus on the (one)
# receiver that recorded detections to define moorings. 
# We assume a deployment period over the duration of detections 
# as the full deployment period & timing of the range test is also unknown. 
pinheiro_moorings <- 
  pinheiro_detections |>  
  janitor::clean_names() |> 
  mutate(dataset = "P") |> 
  rename(receiver_lon = rec_deploy_lon, 
         receiver_lat = rec_deploy_lat) |> 
  distinct(dataset, receiver_sn, receiver_lon, receiver_lat) |> 
  mutate(receiver_start = min(pinheiro_detections$detection_timestamp_utc),
         receiver_end   = max(pinheiro_detections$detection_timestamp_utc), 
         receiver_start = lubridate::force_tz(receiver_start, "UTC"), 
         receiver_end   = lubridate::force_tz(receiver_end, "UTC")
  ) |> 
  select("dataset", 
         "receiver_sn", 
         "receiver_lon", "receiver_lat", 
         "receiver_start", "receiver_end") |> 
  as.data.table()

#### Process Pinheiro (unpublished) detections
pinheiro_detections <- 
  pinheiro_detections |>  
  janitor::clean_names() |>
  mutate(dataset   = "P", 
         timestamp = lubridate::force_tz(detection_timestamp_utc, "UTC")) |> 
  rename(transmitter_id = transmitter, 
         tag_lon        = tag_deploy_lon, 
         tag_lat        = tag_deploy_lat, 
         receiver_lon   = rec_deploy_lon, 
         receiver_lat   = rec_deploy_lat) |> 
  group_by(transmitter_id, tag_lon, tag_lat) |>
  # Assign start/end times based on time of range test(s) using detections
  # This is necessary because unlike the Futia dataset test times are unknown
  mutate(start = min(timestamp), 
         end = max(timestamp)) |> 
  select("dataset", "transmitter_id", "tag_lon", "tag_lat", "start", "end",
         "timestamp", "receiver_sn",  "receiver_lon", "receiver_lat") |> 
  as.data.table()

#### Collate moorings

# Collate moorings
moorings <- 
  rbind(futia_moorings, pinheiro_moorings) |> 
  as.data.table()

# Identify receivers recorded in two places at once
# * This is the case for one receiver: 110541
receivers_in_two_places_at_once <-
  moorings |>
  arrange(receiver_sn, receiver_start) |>
  group_by(receiver_sn) |>
  filter(any(receiver_start < lag(receiver_end), na.rm = TRUE)) |>
  pull(receiver_sn) |>
  unique()

# Assign receiver_id
moorings <- 
  moorings |>
  # Drop receivers in two places at once before defining receiver_id
  filter(!receiver_sn %in% receivers_in_two_places_at_once) |>
  arrange(dataset, receiver_sn, receiver_start) |> 
  mutate(receiver_id = row_number(), .before = 1L) |> 
  as.data.table()

# Check for overlapping deployment periods of the same receiver
# * This validates that no remaining receivers were deployed in 2 places at once
# * NB: Some receivers were deployed/redeployed on same day
# * We could fix this by shifting deployment times by two min
# * But no adjustments occurred during range tests (above) so this shouldn't be necessary
stopifnot(all(moorings$receiver_start < moorings$receiver_end))
lapply(split(moorings, moorings$receiver_sn), function(d) {
  # d <- moorings[receiver_sn == 110541, ]
  print(d$receiver_sn[1])
  d <- d |> arrange(receiver_start)
  if (nrow(d) > 1L) {
    for (i in 1:(nrow(d) - 1)) {
      print(i)
      # The start of the first deployment must be before start of next deployment
      stopifnot(d$receiver_start[i] <= d$receiver_start[i + 1])
      # The end of the first deployment must be before start of next deployment
      stopifnot(d$receiver_end[i] <= d$receiver_start[i + 1])
    }
  }
}) |> invisible()

# Define receiver_x and receiver_y coordinates on map
rxy <- 
  cbind(moorings$receiver_lon, moorings$receiver_lat) |>
  terra::vect(crs = "WGS84") |>
  terra::project(terra::crs(map)) |> 
  terra::geom(df = TRUE)
moorings[, receiver_x := rxy$x]
moorings[, receiver_y := rxy$y]
terra::plot(map)
points(rxy$x, rxy$y, col = "red")

# Note that some receivers were redeployed in the same locations
# * So to assign moorings$receiver_id to detections, 
#   we need to account for receiver_sn, location and time stamp (below)
moorings |> 
  group_by(receiver_sn, receiver_lat, receiver_lon) |> 
  summarise(n = n()) |> 
  ungroup() |>
  pull(n)

#### Collate detections
detections <- 
  rbind(futia_detections, pinheiro_detections) |> 
  # Define receiver_id in detections
  left_join(
    moorings |>
      select(dataset, receiver_id, receiver_sn,
             receiver_lon, receiver_lat, receiver_start, receiver_end),
    by = join_by(
      dataset,
      receiver_sn,
      receiver_lon,
      receiver_lat,
      between(timestamp, receiver_start, receiver_end)
    )
  ) |>
  filter(!is.na(receiver_id)) |> 
  # Define an 'individual_id' column used to distinguish separate range tests
  arrange(dataset, transmitter_id, tag_lon, tag_lat, start, timestamp, receiver_id) |> 
  group_by(dataset, transmitter_id, tag_lon, tag_lat, start) |> 
  mutate(individual_id = cur_group_id(), .before = 1L) |>
  ungroup() |> 
  as.data.table()

#### Collate 'tagging' (range test) information
# This defines test IDs, the location of the range testing tag, and the start/end time
tests <- 
  detections |> 
  distinct(individual_id, dataset, tag_lon, tag_lat, start, end) |>
  as.data.table()
txy <- 
  tests |> 
  select("tag_lon", "tag_lat") |> 
  as.matrix() |> 
  terra::vect(crs = "WGS84") |> 
  terra::project(terra::crs(map)) |> 
  terra::geom(df = TRUE)
tests <- 
  tests |> 
  mutate(tag_x = txy$x, 
         tag_y = txy$y, 
         .after = tag_lat) |> 
  as.data.table()

#### Clean up
moorings <- 
  moorings |> 
  select("receiver_id", 
         "receiver_x", "receiver_y",
         "receiver_start", "receiver_end") |> 
  arrange(receiver_id) |> 
  as.data.table()

detections <- 
  detections |> 
  select("individual_id", "timestamp", "receiver_id") |> 
  arrange(individual_id, timestamp, receiver_id) |> 
  as.data.table()

#### Review datasets
moorings
detections

#### Automated checks
# Validate time zones (UTC)
stopifnot(all(sapply(list(
  futia_moorings$receiver_start,
  futia_moorings$receiver_end,
  pinheiro_moorings$receiver_start,
  pinheiro_moorings$receiver_end,
  moorings$receiver_start, 
  moorings$receiver_end, 
  
  futia_detections$timestamp,
  pinheiro_detections$timestamp, 
  detections$timestamp, 
  
  tests$start, 
  tests$end
), lubridate::tz) == "UTC"))

#### Write to file
qs::qsave(tests, here_input_validation("main", "tests.qs"))
qs::qsave(moorings, here_input_validation("main", "moorings.qs"))
qs::qsave(detections, here_input_validation("main", "detections.qs"))


#### End of code. 
###########################
###########################
