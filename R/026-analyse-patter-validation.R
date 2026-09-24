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
detections    <- qs::qread(here_input_validation("main", "detections.qs"))


###########################
###########################
#### Process data

#### Process iteration
iteration <- 
  iteration |> 
  # Focus on successful runs
  filter(file.exists(file_occupancy)) |> 
  # Add test information i.e., tag_x and tag_y
  left_join(tests, by = "individual_id") |> 
  as.data.table()

#### Process observations
detections <- 
  detections |> 
  # Add test information i.e., tag_x and tag_y
  left_join(tests, by = "individual_id") |> 
  as.data.table()

#### Define land (for clean plots)
land <- terra::erase(
  terra::as.polygons(terra::ext(champlain_utm), crs = terra::crs(champlain_utm)),
  champlain_utm
)


###########################
###########################
#### Analysis

#### Plot spatial distribution of tags/receivers
# TO DO

#### Plot histogram of distances between tags/receivers
# TO DO

#### Plot detection time series
if (FALSE) {
  # Plot detections by dataset
  detections |> 
    ggplot() +
    geom_point(aes(timestamp, factor(individual_id))) +
    facet_wrap(~dataset, scales = "free")
  # Plot detections for each individual
  detections |> 
    ggplot() +
    geom_point(aes(timestamp, receiver_id)) +
    facet_wrap(~individual_id, scales = "free")
}


#### Plot example occurrence distribution with tag location
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
terra::plot(land, col = scales::alpha("lightgrey", 0.8), add = TRUE)
points(test$tag_x, test$tag_y, col = "red3", lwd = 3)
points(m$receiver_x, m$receiver_y, col = "black", pch = 4, cex = 0.5, lwd = 2)
terra::sbar()
# patter::map_hr_home(r, .add = TRUE)
# cf. Distance between distribution centre and tag location
p <- terra::as.points(r, na.rm = TRUE)
terra::crds(p) |>
  apply(2, weighted.mean, w = terra::values(p)[, 1]) |>
  matrix(ncol = 2) |>
  terra::distance(cbind(test$tag_x, test$tag_y), lonlat = FALSE)

#### Plot best occurrence distributions with tag location
table(iteration$sensitivity)
iteration_best <- iteration[sensitivity == "best", ]
png(here_fig("validation", "main", "maps-best.png"), 
    height = 6, width = 10, units = "in", res = 800)
pp <- par(mfrow = c(5, 20), 
          mar = c(0, 0, 0, 0),
          oma = c(0, 0, 0, 0))
pbapply::pblapply(split(iteration_best, iteration_best$index), function(it) {
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
  e  <- terra::ext(r)
  xlim <- as.numeric(e[1:2])
  ylim <- as.numeric(e[3:4])
  # Make map
  terra::plot(r, 
              xlim = xlim, ylim = ylim,
              axes = FALSE, box = FALSE, legend = FALSE,
              mar = NA, buffer = FALSE, 
              font = 2)
  terra::plot(land, col = scales::alpha("lightgrey", 0.5), add = TRUE, lwd = 0.5)
  points(it$tag_x, it$tag_y, col = "red3", lwd = 1.5)
  points(m$receiver_x, m$receiver_y, col = "black", pch = 4, cex = 0.5, lwd = 1)
  # terra::sbar()
  # mtext(side = 3, 
  #       text = it$individual_id, 
  #       line = -1.75, adj = 0.03, font = 2, cex = 1)
  legend(
    "topleft",
    legend = it$individual_id,
    bty = "o",
    bg = scales::alpha("white", 0.9),
    box.col = NA,
    text.font = 2,
    text.width = strwidth(it$individual_id),
    y.intersp = 0.1,
    cex = 1,
    inset = 0.01
  )
  box(lwd = 1)
}) |> invisible()
par(pp)
dev.off()

#### (optional) Compute the distance between the tag location and the distribution centre
if (FALSE) {
  # Compute distances
  distances <- 
    pbapply::pbsapply(split(iteration, iteration$index), function(it) {
      r  <- terra::rast(it$file_occupancy)
      r  <- terra::classify(r, cbind(0, NA))
      r  <- terra::trim(r)
      p  <- terra::as.points(r, na.rm = TRUE)
      centre <- apply(
        terra::crds(p),
        2,
        weighted.mean,
        w = terra::values(p)[, 1]) |> 
        matrix(ncol = 2, byrow = FALSE)
      terra::distance(centre, cbind(test$tag_x, test$tag_y), lonlat = FALSE)
    })
  # Make histogram
  hist(distances)
}


#### End of code.
###########################
###########################