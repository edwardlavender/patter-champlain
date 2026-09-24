###########################
###########################
#### analysis-patter-validation.R

#### Aims
# 1) This script analyses in-situ validation outputs

#### Prerequisites
# 1) Run previous scripts


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
library(spatial.extensions)
library(tictoc)
files_source_r(here_src())

#### Load data
map           <- terra::rast(here_input("map.tif"))
iteration     <- qs::qread(here_input_validation("main", "iteration.qs"))
champlain_utm <- qreadvect(here_input("champlain-utm.qs"))
tests         <- qs::qread(here_input_validation("main", "tests.qs"))
moorings      <- qs::qread(here_input_validation("main", "moorings.qs"))


###########################
###########################
#### Analysis

#### Process coastline
land <- terra::erase(
  terra::as.polygons(terra::ext(champlain_utm), crs = terra::crs(champlain_utm)),
  champlain_utm
)

#### Plot spatial distribution of tags/receivers
# TO DO

#### Plot histogram of distances between tags/receivers
# TO DO

#### Plot example tag location & occurrence distribution
# Define tag location 
it   <- iteration[1, ]
test <- tests[individual_id == it$individual_id, ]
# Define active receivers
m <- 
  moorings |> 
  filter(int_overlaps(interval(receiver_start, receiver_end),
                      interval(it$block_start, it$block_end))) |> 
  as.data.table() 
# Define occurrence distribution, zoomed in a bit
r  <- terra::rast(it$file_occupancy)
r  <- terra::classify(r, cbind(0, NA))
r0 <- r
r  <- terra::trim(r)
e  <- terra::ext(r)
e  <- e + 20000
r  <- terra::crop(r0, e)
# Make map
terra::plot(r, legend = FALSE)
terra::plot(land, col = "lightgrey", add = TRUE)
points(test$tag_x, test$tag_y, col = "red3", lwd = 3)
points(m$receiver_x, m$receiver_y, col = "black", pch = 4, cex = 0.5, lwd = 2)
terra::sbar()
# patter::map_hr_home(r, .add = TRUE)


#### End of code.
###########################
###########################