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
library(lubridate)
library(proj.verse)
files_source_r(here_src())

#### Load data
detections_raw <- qread(here_data_raw("model-obs","futia-et-al-2025",
                                  "range_test_detections_2021-2022.qs"))
metadata <- qread(here_data_raw("model-obs","futia-et-al-2025",
                                "range_test_metadata_2021-2022.qs"))


###########################
###########################
#### Process data
# merge range test detections with associated metadata
detections <- 
  detections_raw |> 
  mutate(receiver_sn = as.integer(receiver_sn)) |> 
  left_join(metadata, relationship = "many-to-many") |>  
  filter(local_time %within% interval(set_dt, pull_dt+75))

# reformat data
detections <- 
  detections |> 
  select(transmitter_id, 
         receiver_id = receiver_sn,
         timestamp = local_time,
         tag_lon = deploy_lon,
         tag_lat = deploy_lat,
         receiver_lon = rec_lon, 
         receiver_lat = rec_lat,
         dB = transmitter_dB
  ) |>
  # Add variables required for modelling  
  mutate(dist = terra::distance(cbind(tag_lon, tag_lat), cbind(receiver_lon, receiver_lat), 
                                lonlat = TRUE, pairwise = TRUE)
  )|> 
  # Cleanup
  select(transmitter_id, timestamp, receiver_id, dB, dist) |> 
  arrange(transmitter_id, timestamp) |>
  as.data.table()

#### Write to file
qs::qsave(detections, here_data("supp", "model-obs", "futia-raw.qs"))


#### End of code. 
###########################
###########################