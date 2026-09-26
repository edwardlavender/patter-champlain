###########################
###########################
#### setup-data-detection-pr-futia.R

#### Aims
# 1) Setup range testing data from Lake Champlain

#### Prerequisites
# 1) Detection probability data provided by Futia & colleagues


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Set global options
Sys.setenv("JULIA_SESSION" = FALSE)

#### Load essential packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = TRUE)
library(ggplot2)
library(lubridate)
library(proj.verse)
files_source_r(here_src())

#### Load data
map        <- terra::rast(here_input("map.tif"))
map_ll     <- terra::project(map, "WGS84")
moorings   <- readRDS(here_data_raw_mf("OriginalReceiverSummary_2013-2017.rds"))
detections <- qs::qread(here_data_raw("model-obs","futia-et-al-2025",
                                      "range_test_detections_2021-2022.qs"))
metadata   <- qs::qread(here_data_raw("model-obs","futia-et-al-2025",
                                      "range_test_metadata_2021-2022.qs"))


###########################
###########################
#### Define detections/non-detections

#### Review time zones
# detections$detection_timestamp_utc is in UTC
# local_time, set_dt and pull_dt are in America/New_York
# For the purpose of the range test analyses, we use local_time
# We keep track of detection_timestamp_utc (below) for the validation analysis
head(detections[, c("detection_timestamp_utc", "local_time")])
lubridate::tz(detections$detection_timestamp_utc)
lubridate::tz(detections$detection_timestamp_utc)
lubridate::tz(metadata$set_dt)
lubridate::tz(metadata$pull_dt)

#### Review range tests
# There are 95 transmitter_id/set_dt combinations
#   (of which 93 are associated with detections, below)
# Metadata contains 172 rows
#   The extra rows are b/c metadata includes rows for multiple receivers
metadata <- as.data.table(metadata)
metadata[, metadata_row_id := seq_len(.N)]
metadata[, metadata_test_id := .GRP, by = c("transmitter_id", "set_dt")]
nrow(metadata)
uniqueN(metadata$metadata_row_id)
uniqueN(metadata$metadata_test_id)
# This code shows that the same test ID includes multiple receivers 
metadata |> 
  group_by(transmitter_id, set_dt) |> 
  mutate(n = n(), group = cur_group_id()) |> 
  filter(n > 1) |> 
  select(group, location, receiver_sn, rec_lat, rec_lon, transmitter_id, set_dt) |> 
  arrange(group, receiver_sn) |> 
  as.data.table()
# Check metadata for example tags 
metadata |> 
  filter(transmitter_id %in% c(4331, 4632)) |> 
  distinct(transmitter_id, tag_type, transmitter_dB)

#### Merge range test detections with associated metadata
# The set and pull times are exact to the second that the range test tags 
# were set and deployed, so there should not be any more detections than 
# what is possible in that window. I added the extension to the recover times 
# as some receivers had minor clock drift that was not completely accounted 
# for in the time corrections, so there were some deployments where the detections 
# did not line up perfectly within the recorded deployment period. Adding 75 
# seconds allowed me to capture all detections but was a short enough interval 
# that it did not overlap with the following deployment. 
detections <- 
  detections |> 
  mutate(receiver_sn = as.integer(receiver_sn)) |> 
  left_join(metadata, relationship = "many-to-many") |>  
  filter(local_time %within% interval(set_dt, pull_dt + 75))
# 131 rows in metadata are associated with detections 
uniqueN(detections$metadata_row_id)
# 93/95 range tests associated with detections
uniqueN(detections$metadata_test_id)

#### Examine range test locations (compared to moorings)
# Define range test locations 
range_test_locs <- 
  detections |> 
  distinct(location, .keep_all = TRUE) |> 
  select(location, lon = rec_lon, lat = rec_lat) 
# Define mooring locations
moorings_locs <- 
  moorings |> 
  distinct(StationName, .keep_all = TRUE) |> 
  select(location = StationName, lon = deploy_long, lat = deploy_lat)
# Map
png(here_fig("model", "model-obs", "champlain-range-tests.png"), 
    height = 10, width = 10, units = "in", res = 600)
pp <- par(mfrow = c(1, 2))
terra::plot(map_ll)
points(range_test_locs$lon, range_test_locs$lat)
basicPlotteR::addTextLabels(range_test_locs$lon, 
                            range_test_locs$lat, 
                            range_test_locs$location)
terra::plot(map_ll)
points(moorings_locs$lon, moorings$locs$lat)
basicPlotteR::addTextLabels(moorings_locs$lon, 
                            moorings_locs$lat, 
                            moorings_locs$location)
dev.off()
# Summary stats
# > Range tests were conducted for 9 receivers
# > in 9 broad locations (for 17 unique coordinates)
# > In the same broad area as many of the receivers
# > but locations only overlap for 2/9 locations
length(unique(detections$receiver))
detections |> 
  distinct(rec_lon, rec_lat, .keep_all = TRUE) |> 
  nrow()
table(range_test_locs$location %in% moorings_locs$location)

#### Reformat data
detections <- 
  validation <- 
  detections |> 
  select(transmitter_id, 
         receiver_id = receiver_sn,
         # Keep track of detection_timestamp_utc for validation analysis
         detection_timestamp_utc,
         # For range tests, use local_time 
         timestamp = local_time,
         tag_lon = deploy_lon,
         tag_lat = deploy_lat,
         receiver_lon = rec_lon, 
         receiver_lat = rec_lat,
         dB = transmitter_dB, 
         delay = tag_delay, 
         start = set_dt, 
         end = pull_dt
  ) |>
  # Add variables required for modelling  
  mutate(dist = terra::distance(cbind(tag_lon, tag_lat), 
                                cbind(receiver_lon, receiver_lat), 
                                lonlat = TRUE, pairwise = TRUE)
  ) |> 
  as.data.table()

# Cleanup
detections <-
  copy(detections) |>
  arrange(transmitter_id, timestamp) |>
  select(transmitter_id, start, end, timestamp, receiver_id, dB, delay, dist) |> 
  as.data.table()

#### Summarise raw dataset
nrow(detections)
length(unique(detections$receiver_id))
range(detections$timestamp)


###########################
###########################
#### Define detection counts

#### Create summarises of observed/expected number of detections for modelling
# Tags were deployed for approximately 30 mins with a high ping frequency (7 or 15 s)
# Each tag deployment is called an experiment
# For each experiment, we compute the observed/expected number of detections 
# (by receiver)
dcounts <- 
  detections |> 
  # Define experiments
  mutate(experiment = paste(transmitter_id, start), 
         experiment = factor(experiment, levels = unique(experiment)), 
         experiment = as.integer(experiment)) |>
  # Compute observed number of detections at each receiver per experiment
  group_by(experiment, receiver_id) |> 
  mutate(observed = n()) |> 
  slice(1L) |>
  ungroup() |> 
  # Compute expected number of transmissions per experiment
  # * Note the round() to ensure the expected number of transmissions is Int
  mutate(expected = round(as.numeric(difftime(end, start, units = "secs")) / delay)) |>
  # Use success/failure for glm
  mutate(
    observed = as.integer(observed), 
    expected = as.integer(expected),
    success = observed, 
    failure = expected - observed, 
    prop = success / (success + failure)) |>
  # Cleanup
  select(transmitter_id, timestamp, prop, success, failure, dB, dist) |> 
  arrange(transmitter_id, timestamp) |>
  as.data.table()

#### Quick validation
# The number of successes and failures should be integers
str(dcounts)
# The expected number of transmissions should be >= observed number
table(dcounts$failure < 0)
# Plot detection efficacy
dcounts |> 
  ggplot(aes(dist, prop)) + 
  geom_point() + 
  geom_smooth()

#### Write to file
qs::qsave(dcounts, here_data("supp", "model-obs", "futia-raw.qs"))
qs::qsave(validation, here_data("supp", "model-obs", "futia-raw-validation.qs"))

#### End of code. 
###########################
###########################