###########################
###########################
#### setup-data-detection-pr.R

#### Aims
# 1) Setup detection probability data

#### Prerequisites
# 1) Detection probability data provided by Klinard et al. (2019)


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
library(proj.verse)
files_source_r(here_src())

#### Load data
detections <- fread(here_data_raw("model-obs", "klinard-et-al-2019", 
                                  "Aug2015-Jun2016_Range_Detections.csv"))


###########################
###########################
#### Klinard et al. (2019) Methods

# Study period: 
# The complete receiver array, including the 85 receivers
# from the bloater telemetry project and the five receivers
# for range testing, was deployed from 22 October, 2015
# to 23 May, 2016 (215 days). To ensure consistency across
# detection distances and probabilities, only detections for
# these dates were used in analyses. 

# Transmitters:
# * 3x V9-2x 69-kHz range tags 
#   (power output 145 dB, nominal delay 1800 s, random interval 1750–1850 s), 
# * 1x V13-169-kHz range tag 
#   (power output 153 dB, nominal delay 1800 s)
# * four V16-6X 69-kHz range tags 
#   (power output 158 dB, nominal delay 1800 s)

# Depth zones:
# To examine spatial variability in DE across tag types and
# depths, detection data were separated into five _categories_: 
# deep V9, shallow V9, deep V13, deep V16, and shallow V16. 

# Tag deployment depths:
# * The deep group of tags (two** V9, one V13, one V16) was situated
#   below the thermocline at a depth of 50 m
#   ** There is a typo in the manuscript, two V9s were deployed at 50 m
#      See Table 1
# * The shallow group (one V9, one V16) was above the thermocline
#   at a depth of 11 m to evaluate the impact of tag depth and
#   thermal stratification on DE. 

# cf. Lake Champlain system:
# * Receivers deployed in depths of 3 - 50 m (see setup-data-detection.R)
# * V13 147 dB tags used

# We will select tags in the following categories:
# * V9 11 m (one tag)  
# * V9 50 m (two tags tag)
# * V13 50 m (one tag)
# (These tags bound the characteristics of the tags used in our system)

# Transmitter IDs
# * See LakeOntario_2015-16_range_tag_specs
# * This defines transmitter models (V9, V13) and IDs (but not depths)
# * This list is the subset of (a) relevant and (b) properly functional tags
# * We identified the depth category for each tag by:
#   - Compute max detection range (below)
#   - Compare to Table 1 in Klinard et al. (2019) 

# 30838: V13 (deep)
# 57347: V9 (deep)
# 57349: V9 (deep)
# 57350: V9 (shallow)

# Analyses were performed separately for each tag category. For each tag and
# receiver combination (n=720), DE was calculated for
# each day of deployment by dividing the number of detections 
# by the expected number of transmissions per day
# (48 for a nominal transmission interval of 1800 s). 
# Daily DE was used to estimate DE for the entire study period
# using generalized additive mixed models (GAMMs)


###########################
###########################
#### Setup data

#### Define transmitter IDs (see above)
transmitter_ids <- c("57347", "57349", "57350", "30838")
transmitter_dB <- c(145, 145, 145, 153)
transmitters    <- data.table(id = transmitter_ids,
                              dB = transmitter_dB)
stopifnot(all(transmitter_ids %in% detections$transmitter_id))

#### Identify transmitter depth categories (see above)
detections |> 
  filter(transmitter_id %in% transmitter_ids) |> 
  mutate(dist = terra::distance(cbind(tag_lon, tag_lat), cbind(longitude, latitude), 
                                lonlat = TRUE, pairwise = TRUE)) |> 
  group_by(transmitter_id) |>
  summarise(max(dist))

#### Process detections 
detections <- 
  detections |> 
  select(transmitter_id, 
         receiver_id = receiver_sn,
         timestamp = datetime_UTC,
         tag_lon, tag_lat,
         receiver_lon = longitude, receiver_lat = latitude) |> 
  # Focus on the relevant time window of detections
  filter(timestamp >= as.POSIXct("2015-10-22 00:00:00", tz = "UTC")) |> 
  filter(timestamp <= as.POSIXct("2016-05-23 00:00:00", tz = "UTC")) |> 
  # Focus on relevant transmitters
  filter(transmitter_id %in% transmitter_ids) |> 
  # Add covariates for models
  mutate(dist = terra::distance(cbind(tag_lon, tag_lat), cbind(receiver_lon, receiver_lat), 
                                lonlat = TRUE, pairwise = TRUE)) |> 
  # Cleanup
  arrange(transmitter_id, timestamp) |>
  as.data.table()

#### Create daily summarises of observed/expected number of detections for modelling
klinard <- 
  detections |> 
  # Compute observed number of detections per transmitter/receiver/day
  mutate(timestamp = lubridate::floor_date(timestamp, "days")) |> 
  group_by(transmitter_id, receiver_id, timestamp) |> 
  mutate(observed = n()) |> 
  slice(1L) |>
  ungroup() |> 
  # Compute expected number of transmissions per day
  mutate(expected = (24 * 60 * 60) / 1800) |>
  # Use success/failure for glm
  mutate(success = observed, failure = expected - observed, 
         prop = success / (success + failure)) |>
  # Cleanup
  select(transmitter_id, timestamp, prop, success, failure, dist) |> 
  arrange(transmitter_id, timestamp) |>
  as.data.table()

#### Compute weights
# We assign weights so that on average, a GLM of detection probability
# behaves as though all detections came from 147 dB tags (Lake Champlain tags)
# Thus, we upweight 147 dB tags & downweight 153 tags, 
# accounting for the number of observations. 
# Since we have more measurements from 147 dB tags, this utimately involves
# a downweighting of those measurements and an upweighting of the 153 measurements. 
klinard <- 
  klinard |> 
  mutate(dB = transmitters$dB[match(transmitter_id, transmitters$id)]) |>
  group_by(dB) %>%
  mutate(
    n_obs = n(),
    # Compute weight, accounting for dB scaling and number of observations
    w = 10^((147 - first(dB)) / 10) / n_obs
  ) %>%
  ungroup() %>%
  # Normalise weights so average weight is one 
  mutate(w = w / mean(w)) |> 
  as.data.table()

# Examine weights
table(klinard$w)

#### Checks
# The expected number of transmissions should be >= observed number
# There are a few cases where that is not the case
table(sort(klinard$failure))
klinard[failure < 0, c("success", "failure") := .(48, 0)]
table(sort(klinard$failure))
table(sort(klinard$success))

#### Write to file
qs::qsave(klinard, here_data("supp", "model-obs", "klinard.qs"))


#### End of code. 
###########################
###########################