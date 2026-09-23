###########################
###########################
#### validation.R

#### Aims
# 1) This script supports in-field validation of patter outputs

#### Prerequisites
# 1) NA


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

#### Load data
val <-
  qs2::qs_read(
    here_data_raw("mfutia", 
                  "validation", 
                  "Champlain_RangeTestDetections_Summer2016_mhf.qs2"))


###########################
###########################
#### Examine data

#### Quick checks
val <- as.data.table(val)
head(val)
table(val$transmitter)
table(val$receiver_sn)

#### Examine the time frame of the range test
val |> 
  group_by(transmitter, tag_deployLat, tag_deployLon) |> 
  reframe(range(detection_timestamp_utc)) |> 
  as.data.table()

#### Plot detection time series 
val |> 
  ggplot(aes(detection_timestamp_utc, 
             factor(transmitter), 
             colour = factor(receiver_sn))) + 
  geom_point()

#### Plot locations of range test, including tags and receivers
leaflet() |>
  addTiles() |>
  addCircleMarkers(
    data = distinct(val, transmitter, tag_deployLat, tag_deployLon),
    lng = ~tag_deployLon,
    lat = ~tag_deployLat,
    label = ~paste("Tag:", transmitter),
    popup = ~paste("Tag:", transmitter),
    radius = 5,
    color = "blue"
  ) |>
  addCircleMarkers(
    data = distinct(val, receiver_sn, rec_deployLat, rec_deployLon),
    lng = ~rec_deployLon,
    lat = ~rec_deployLat,
    label = ~paste("Receiver:", receiver_sn),
    popup = ~paste("Receiver:", receiver_sn),
    radius = 5,
    color = "red"
  )

#### Examine distances between tags and receivers (m)
pts <- 
  val |>
  summarise(
    tag1_lon = first(tag_deployLon),
    tag1_lat = first(tag_deployLat),
    rec_lon = first(rec_deployLon),
    rec_lat = first(rec_deployLat),
    .by = transmitter
  )
xy <- rbind(
  c(pts$tag1_lon[1], pts$tag1_lat[1]),
  c(pts$tag1_lon[2], pts$tag1_lat[2]),
  c(pts$rec_lon[1], pts$rec_lat[1])
)
p <- terra::vect(xy, crs = "EPSG:4326")
terra::distance(p)


#### End of code. 
###########################
###########################
