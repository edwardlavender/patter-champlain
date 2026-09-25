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
library(proj.verse)
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
  # Add detection summary statistics (for plot)
  left_join(
    detections |> 
      left_join(
        tests |> 
          select(individual_id, start, end),
        by = "individual_id"
      ) |> 
      group_by(individual_id) |> 
      arrange(timestamp, .by_group = TRUE) |> 
      summarise(
        # Number of detections for each individual
        detection_count = n(),
        # Longest period without detection (mins)
        detection_gap_max = max(
          c(
            difftime(first(timestamp), first(start), units = "mins"),
            Tools4ETS::serial_difference(timestamp, units = "mins"),
            difftime(first(end), last(timestamp), units = "mins")
          ), na.rm = TRUE
        ),
        detection_gap_max = as.numeric(round(detection_gap_max)),
        # Duration of the test 
        test_duration = as.numeric(round(difftime(max(end), min(start), units = "mins")))
      ) |> 
      as.data.table(),
    by = "individual_id"
  )
# Review test durations
iteration |> 
  distinct(individual_id, start, end, test_duration) |> 
  arrange(test_duration)

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

#### Test summary statistics
# Number of range-testing tags
tests |> 
  group_by(dataset) |> 
  summarise(n = n())
# Number of unique tag deployment locations
tests |>
  group_by(dataset) |>
  summarise(uniqueN(paste(tag_lon, tag_lat)))
# Time period of range tests 
tests |> 
  group_by(dataset) |> 
  summarise(min(start), max(end))
# Duration of tests 
tests |> 
  group_by(dataset) |> 
  reframe(utils.add::basic_stats(as.numeric(difftime(end, start, units = "days"))))
# Distance from each test to the nearest active receiver (m)
tests <- 
  tests |> 
  left_join(
    moorings |> 
      select(receiver_x, receiver_y, receiver_start, receiver_end),
    by = join_by(
      start <= receiver_end,
      end >= receiver_start
    )
  ) |> 
  mutate(
    distance_nearest_receiver = sqrt(
      (tag_x - receiver_x)^2 + 
        (tag_y - receiver_y)^2
    )
  ) |> 
  summarise(
    distance_nearest_receiver = min(distance_nearest_receiver),
    .by = all_of(names(tests))
  ) |> 
  as.data.table()
# Summarise distances
tests |> 
  group_by(dataset) |>
  reframe(utils.add::basic_stats(distance_nearest_receiver))
# Summarise distances for P
tests |> 
  filter(dataset == "P") |> 
  group_by(individual_id) |> 
  summarise(unique(distance_nearest_receiver))

#### Detection summary statistics
# Detection period duration
detections |>
  group_by(individual_id) |> 
  mutate(length = difftime(max(timestamp), min(timestamp), units = "mins")) |>
  slice(1L) |> 
  group_by(dataset) |> 
  reframe(utils.add::basic_stats(length))
# Detection period duration for P
detections |>
  filter(dataset == "P") |> 
  group_by(individual_id) |> 
  summarise(length = difftime(max(timestamp), min(timestamp), units = "mins"))
# Number of detections per individual
detections |>
  group_by(individual_id) |> 
  mutate(n = n()) |>
  slice(1L) |> 
  group_by(dataset) |> 
  reframe(utils.add::basic_stats(n))
# Number of detections per individual for P dataset
detections |>
  filter(dataset == "P") |> 
  group_by(individual_id) |> 
  summarise(n = n()) 
# Number of detections overall
detections |>
  group_by(dataset) |> 
  summarise(n = n())
# Detection gap duration 
detections |>
  group_by(individual_id) |> 
  arrange(timestamp, .by_group = TRUE) |> 
  mutate(gap = Tools4ETS::serial_difference(timestamp, units = "mins")) |>
  ungroup() |> 
  filter(!is.na(gap)) |> 
  group_by(dataset) |> 
  reframe(utils.add::basic_stats(gap))
# Detection gap duration for P 
detections |>
  filter(dataset == "P") |> 
  group_by(individual_id) |> 
  arrange(timestamp, .by_group = TRUE) |> 
  mutate(gap = Tools4ETS::serial_difference(timestamp, units = "mins")) |>
  filter(!is.na(gap)) |> 
  reframe(utils.add::basic_stats(gap))
# Relationship between gap duration & distance from nearest receiver
iteration |> 
  filter(sensitivity == "best") |> 
  distinct(individual_id, dataset, detection_gap_max) |> 
  left_join(
    tests |> select(individual_id, distance_nearest_receiver),
    by = "individual_id"
  ) |> 
  ggplot(aes(distance_nearest_receiver, detection_gap_max)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~dataset, scales = "free") +
  labs(
    x = "Distance to nearest receiver (m)",
    y = "Maximum detection gap (min)"
  )

#### Plot spatial distribution of tags/receivers
png(here_fig("validation", "main", "champlain-range-tests.png"), 
    height = 10, width = 3, units = "in", res = 800)
terra::plot(land, col = scales::alpha("dimgrey", 0.3))
# points(moorings$receiver_x, moorings$receiver_y, 
#        pch = 4, cex = 0.35)
points(tests$tag_x[tests$dataset == "F"], 
       tests$tag_y[tests$dataset == "F"], 
       pch = 4, col = "red", cex = 0.2)
dev.off()

#### Plot histogram of distances between tags/receivers (~1 s)
distances <- 
  pbapply::pblapply(split(tests, tests$individual_id), function(test) {
    # Define active receivers
    m <- 
      moorings |> 
      filter(int_overlaps(interval(receiver_start, receiver_end),
                          interval(test$start, test$end))) |> 
      as.data.table() 
    # Compute distances between active receivers and tag location
    distances <- terra::distance(cbind(m$receiver_x, m$receiver_y), 
                                 cbind(test$tag_x, test$tag_y),
                                 lonlat = FALSE)
    # distances <- distances[distances < 10000]
    data.frame(individual_id = test$individual_id, 
               distance = distances)
  }) |> rbindlist()
# Summarise the average distance of a range testing tag from an active receiver
utils.add::basic_stats(distances$distance)
# As above but focusing on receivers within ~10 km
utils.add::basic_stats(distances$distance[distances$distance < 10000])
# Visualise histogram 
hist(distances$distance)

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
terra::plot(land, col = scales::alpha("dimgrey", 0.3), add = TRUE)
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

#### Plot occurrence distributions for independent range test (including sensitivity)
if (TRUE) {
  # Define iterations 
  iteration_selected <- iteration[dataset == "P", ]
  nr <- uniqueN(iteration_selected$individual_id)
  nc <- uniqueN(iteration_selected$sensitivity_label)
  # Define map aspect ratio
  it <- iteration_selected[1, ]
  buffer <- 
    cbind(it$tag_x, it$tag_y) |> 
    terra::vect(crs = terra::crs(map)) |> 
    terra::buffer(width = 15e3)
  r <- 
    terra::rast(it$file_occupancy) |> 
    terra::crop(buffer)
  e <- terra::ext(r)
  map_asp <- (e$xmax - e$xmin) / (e$ymax - e$ymin)
  # Define device dimensions to match map aspect ratio
  height <- 2.5
  width <- height * nc / nr * map_asp
  # Make plot 
  png(here_fig("validation", "main", "maps-pinheiro-full.png"), 
      height = height, width = width, units = "in", res = 800)
  pp <- par(mfrow = c(nr, nc), 
            mar = c(0, 0, 0, 0),
            oma = c(0, 0, 0, 0))
  # Iterate over individuals and make plot
  pbapply::pblapply(
    split(
      iteration_selected,
      interaction(
        iteration_selected$individual_id,
        iteration_selected$sensitivity_label,
        drop = TRUE,
        lex.order = TRUE
      )
    ), 
    function(it) {
      # Define active receivers
      m <- 
        moorings |> 
        filter(int_overlaps(interval(receiver_start, receiver_end),
                            interval(it$block_start, it$block_end))) |> 
        as.data.table() 
      # Define tag region
      buffer <- 
        cbind(it$tag_x, it$tag_y) |> 
        terra::vect(crs = terra::crs(map)) |> 
        terra::buffer(width = 15e3)
      # Define occurrence distribution around tag
      r0 <- terra::rast(it$file_occupancy)
      r  <- terra::classify(r0, cbind(0, NA))
      r  <- terra::crop(r, buffer)
      # Make map
      terra::plot(
        r, 
        axes = FALSE, box = FALSE, legend = FALSE,
        mar = NA, buffer = FALSE, 
        font = 2
      )
      # patter::map_hr_home(r0, .add = TRUE)
      terra::plot(land, col = "white", add = TRUE, border = NA)
      terra::plot(
        land, 
        col = scales::alpha("lightgrey", 0.5), 
        add = TRUE, lwd = 0.5
      )
      points(it$tag_x, it$tag_y, col = "red3", lwd = 1.5)
      points(
        m$receiver_x, m$receiver_y, 
        col = "black", pch = 4, cex = 0.5, lwd = 1
      )
      terra::sbar(
        d = 5000, xy = "bottomleft", labels = "", 
        lonlat = FALSE, halo = FALSE
      )
      # Add panel label
      usr <- par("usr")
      yadj <- 0.2
      rect(
        xleft   = usr[1],
        ybottom = usr[4] - yadj * diff(usr[3:4]),
        xright  = usr[2],
        ytop    = usr[4],
        col     = scales::alpha("dimgrey", 0.7),
        border  = NA
      )
      text(
        x = mean(usr[1:2]),
        y = usr[4] - yadj / 2 * diff(usr[3:4]),
        labels = paste0(
          match(
            it$individual_id,
            sort(unique(iteration_selected$individual_id))
          ), " | ",
          it$sensitivity_label
        ),
        font = 2,
        cex = 1.2
      )
      
      box(lwd = 1)
    }
  ) |> invisible()
  par(pp)
  dev.off()
}

#### Plot occurrence distributions for non-independent range tests (best only)
if (TRUE) {
  png(here_fig("validation", "main", "maps-futia-best.png"), 
      height = 10, width = 10, units = "in", res = 800)
  pp <- par(mfrow = c(10, 10), 
            mar = c(0, 0, 0, 0),
            oma = c(0, 0, 0, 0))
  iteration_selected <- iteration[sensitivity == "best" &  dataset == "F", ]
  nrow(iteration_selected)
  pbapply::pblapply(split(iteration_selected, iteration_selected$index), function(it) {
    # Define active receivers
    m <- 
      moorings |> 
      filter(int_overlaps(interval(receiver_start, receiver_end),
                          interval(it$block_start, it$block_end))) |> 
      as.data.table() 
    # Define tag region
    buffer <- 
      cbind(it$tag_x, it$tag_y) |> 
      terra::vect(crs = terra::crs(map)) |> 
      terra::buffer(width = 4000)
    # Define occurrence distribution around tag
    r0 <- terra::rast(it$file_occupancy)
    r  <- terra::classify(r0, cbind(0, NA))
    r  <- terra::crop(r, buffer)
    e  <- terra::ext(buffer)
    xlim <- as.numeric(e[1:2])
    ylim <- as.numeric(e[3:4])
    # Make map
    terra::plot(r, 
                xlim = xlim, ylim = ylim,
                axes = FALSE, box = FALSE, legend = FALSE,
                mar = NA, buffer = FALSE, 
                font = 2)
    # patter::map_hr_home(r0, .add = TRUE)
    terra::plot(land, col = "white", add = TRUE, border = NA)
    terra::plot(land, col = scales::alpha("lightgrey", 0.5), add = TRUE, lwd = 0.5)
    points(it$tag_x, it$tag_y, col = "red3", lwd = 1.5)
    points(m$receiver_x, m$receiver_y, col = "black", pch = 4, cex = 0.5, lwd = 1)
    terra::sbar(d = 1000, xy = "bottomleft", labels = "", lonlat = FALSE, halo = FALSE)
    # terra::sbar()
    # mtext(side = 3, 
    #       text = it$individual_id, 
    #       line = -1.75, adj = 0.03, font = 2, cex = 1)
    usr <- par("usr")
    yadj <- 0.2
    rect(
      xleft   = usr[1],
      ybottom = usr[4] - yadj * diff(usr[3:4]),
      xright  = usr[2],
      ytop    = usr[4],
      col     = scales::alpha("dimgrey", 0.7),
      border  = NA
    )
    text(
      x = mean(usr[1:2]),
      y = usr[4] - yadj / 2 * diff(usr[3:4]),
      labels = paste0(
        it$individual_id, " (",
        it$test_duration, ", ",
        it$detection_count, ", ",
        it$detection_gap_max, ")"
      ),
      font = 2,
      cex = 1.2
    )
    box(lwd = 1)
  }) |> invisible()
  par(pp)
  dev.off()
}

#### (optional) Compute the distance between the tag location and the distribution centre
if (TRUE) {
  
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
  iteration[, distance_error := distances]
  
  # Make histogram
  hist(iteration$distance_error)
  
  # Relate distances to explanatory variables e.g., number of detections
  ggplot(iteration, aes(detection_count, distances)) + 
    geom_point() + 
    geom_smooth() + 
    facet_wrap(~dataset, scales = "free")
  ggplot(iteration, aes(detection_gap_max, distances)) + 
    geom_point() + 
    geom_smooth() + 
    facet_wrap(~dataset, scales = "free")
  
}

#### End of code.
###########################
###########################