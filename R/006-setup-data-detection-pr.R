###########################
###########################
#### setup-data-detection-pr.R

#### Aims
# 1) Collates detection probability datasets & prepares data for modelling 

#### Prerequisites
# 1) Process detection probability datasets


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
# TO DO Include Futia data
klinard <- qs::qread(here_data("supp", "model-obs", "klinard-raw.qs"))


###########################
###########################
#### Collate datasets for modelling 

#### Collate datasets
# TO DO Collate datasets
head(klinard)

#### Create daily summarises of observed/expected number of detections for modelling
klinard <- 
  klinard |> 
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
  select(transmitter_id, timestamp, prop, success, failure, dB, dist) |> 
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