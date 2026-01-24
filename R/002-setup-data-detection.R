###########################
###########################
#### setup-data-detection.R

#### Aims
# 1) Sets up detection data for analysis with patter

#### Prerequisites
# 1) Raw data provided by M. Futia 


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Set global options
Sys.setenv("JULIA_SESSION" = FALSE)
set.seed(1L)

#### Load essential packages
library(proj.verse)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(leaflet)
library(lubridate)
library(prettyGraphics)
library(tictoc)
files_source_r(here_src())

#### Load data 
# lkt_detections_2013-2017.rds: raw detections (93 fish)
# lkt_detections_2013-2017_filtered.qs: filtered detections (as in Futia et al., 2024)
map                 <- terra::rast(here_input("map.tif"))
epsg_utm            <- qs::qread(here_input("epsg-utm.qs"))
map_bbox            <- qs::qread(here_input("map-bbox.qs"))
detections          <- readRDS(here_data_raw_mf("lkt_detections_2013-2017.rds"))
detections_filtered <- qs::qread(here_data_raw_mf("lkt_detections_2013-2017_filtered.qs"))
moorings            <- readRDS(here_data_raw_mf("OriginalReceiverSummary_2013-2017.rds"))
surgery             <- fread(here_data_raw("mfutia", "model_comparison", "surgery_log.csv"))
champlain_utm       <- qs::qread(here_input("champlain-utm.qs"))


###########################
###########################
#### Data summary

#### Summary
# 93 fish were captured, tagged and detected
# These fish are defined in surgery and detections
# In the filtered detection dataset, there are 79 individuals
# For now, the tagged/detected fish data summaries are based on the raw data
# (optional) TO DO Reorganise this code
# * Process detections first
# * Compute summary statistics only for tagged fish ultimately included in the analysis
# * See also prepare-analysis.R
# This is not currently implemented b/c in the manuscript we provide the full
# dataset size & then explain the filtering steps. 

# In the filtered detection dataset, there are 79 individuals
# For now, the tagged fish summary includes 93 fish 
# (optional) TO DO Move this code around so detections is processed first
# and the supporting tables are only produced for analysed fish 


###########################
###########################
#### Identify tagged fish

#### Examine tagging events
head(surgery)
nrow(surgery)
utils.add::basic_stats(surgery$length / 1000)
range(as.Date(surgery$cap_date, format = "%m/%d/%Y"))
unique(surgery$tag_type)

#### Clean surgery (tagging) events
surgery <- 
  surgery |> 
  mutate(date = as.Date(surgery$cap_date, format = "%m/%d/%Y"),
         lon = detections$deploy_long[match(animal_id, detections$animal_id)], 
         lat = detections$deploy_lat[match(animal_id, detections$animal_id)],
         row = row_number()) |> 
  as.data.table()

#### Write QGIS
surgery |> 
  group_by(lon, lat) |>
  summarise(n = n()) |>
  ungroup() |> 
  as.data.frame() |> 
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) |> 
  sf::st_write(here_fig("local", "qgis", "tagging.shp"), append = FALSE)

#### Write tidy table
surgery |> 
  select(Row = row, 
         ID = animal_id, 
         `Longitude (°)` = lon, 
         `Latitude (°)` = lat, 
         Date = date, 
         Sex = sex, 
         `Total length (mm)` = length) |> 
  tidy_numbers(digits = c(0, 0, 4, 4, 0)) |> 
  tidy_write(here_fig("tables", "fish.txt"))


###########################
###########################
#### Identify detected fish 

#### Define fish (id, size, tagging location)
# Define fish 
fish <- 
  detections |> 
  group_by(animal_id) |> 
  summarise(individual_id = animal_id[1], 
            len = length[1] / 1000, 
            site = cap_site[1],
            lat = deploy_lat[1], 
            lon = deploy_long[1],
            # Checks 
            nlen = n_distinct(length), 
            nlat = n_distinct(lat), 
            nlon = n_distinct(lon)
            ) |> 
  as.data.table()
# Define tagging locations (UTM)
# * This code requires internet
xy <- 
  cbind(fish$lon, fish$lat) |> 
  terra::vect(crs = "EPSG:4326") |> 
  terra::project(epsg_utm) |> 
  terra::crds()
stopifnot(nrow(xy) > 0L)
fish[, x := xy[, 1]]
fish[, y := xy[, 2]]
# Check tagging locations
stopifnot(all(!is.na(terra::extract(map, xy)[, 1])))
terra::plot(map)
points(xy)
# Check tagging sites
table(fish$site)
terra::plot(map)
points(xy[fish$site == "Grand Isle", ])
terra::plot(map)
points(xy[fish$site == "Split Rock", ])
# Check tagging dates
# * Note that fish were tagged at different times
# * If we analyse the data in blocks e.g., months, we need to account for this
# * TO DO Confirm date format %m/%d/%Y
# * We could check that each fish is only associated with detections after tagging
range(detections$detection_timestamp_utc)
range(as.Date(surgery$cap_date, format = "%m/%d/%Y"))

#### Checks
# Each individual is associated with one length (presumably length @ tagging)
# deploy_lat and deploy_long seem to refer to fish tagging locations
stopifnot(all(fish$nlen == 1L))
stopifnot(all(fish$nlat == 1L))
stopifnot(all(fish$nlon == 1L))

#### Clean up
fish <- 
  fish |> 
  select(individual_id, len, site, x, y, lon, lat) |> 
  as.data.table()

#### Comments
# length is total length (mm)
# Hansen et al. 2022 provide an equation to convert TL to FL for lake trout 
# (FL = 0.9143*TL-8.2772)


###########################
###########################
#### Prepare moorings

#### Define receiver coordinates (UTM)
# Receivers were deployed in 31 stations
unique(moorings$StationName)
# With slight changes in position
# > There are 143 unique pairs of lon, lat coordinates
# > Out of 153 receiver deployments
moorings |>
  group_by(deploy_long, deploy_lat) |> 
  slice(1L) |>
  ungroup() |> 
  nrow()
nrow(moorings)
# Define moorings in UTM
# > Note that re-projection requires an internet connection! 
head(moorings)
rxy <- 
  cbind(moorings$deploy_lon, moorings$deploy_lat) |> 
  terra::vect(crs = "EPSG:4326") |> 
  terra::project(epsg_utm) |>
  terra::crds()
stopifnot(nrow(rxy) > 0L)
stopifnot(all(!is.na(terra::extract(map, rxy)$map_value)))

#### Receiver depths 
utils.add::basic_stats(moorings$depth, na.rm = TRUE)
# min  mean median  max   sd   IQR  MAD
# 1 3.4 13.56  11.85 45.7 9.06 11.85 8.82

#### Process moorings for mapping (QGIS)
moorings |>
  select(lat = deploy_lat, lon = deploy_long) |> 
  group_by(lon, lat) |> 
  summarise(n = n()) |> 
  ungroup() |> 
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) |> 
  # Set as polygon with specified radius 
  sf::st_buffer(dist = 1000, nQuadSegs = 1000) |> 
  sf::st_write(here_fig("local", "qgis", "moorings.shp"), append = FALSE)
  
#### Process moorings for tidy table
str(moorings)
moorings |>
  as_tibble() |> 
  janitor::clean_names() |> 
  select(receiver_station = station_name, 
         receiver_sn,
         start = deploy_date_time, end = recover_date_time, 
         lat = deploy_lat, lon = deploy_long, depth) |> 
  arrange(receiver_station, start, receiver_sn) |> 
  mutate(ID = as.character(row_number()), 
         lon = plyr::round_any(lon, 0.0001),
         lon = add_lagging_point_zero(lon, 4), 
         lat = plyr::round_any(lat, 0.0001), 
         lat = add_lagging_point_zero(lat, 4), 
         depth = plyr::round_any(depth, 0.1),
         depth = add_lagging_point_zero(depth, 1),
         ) |> 
  select(ID, 
         Station = receiver_station, 
         Receiver = receiver_sn,
         Start = start,
         End = end, 
         `Longitude (°)` = lon, 
         `Latitude (°)` = lat, 
         `Depth (m)` = depth) |> 
  tidy_write(here_fig("tables", "moorings.txt"))
  
#### Process moorings for modelling
moorings_raw <- 
  moorings |> 
  janitor::clean_names() |> 
  mutate(receiver_station = station_name, 
         receiver_id = row_number(),
         receiver_sn = as.integer(as.character(receiver_sn)),
         receiver_start = as.POSIXct(paste0(deploy_date_time, "00:00:00"), tz = "UTC"), 
         receiver_end = as.POSIXct(paste0(recover_date_time, "00:00:00"), tz = "UTC"), 
         receiver_int = lubridate::interval(receiver_start, receiver_end),
         receiver_x = rxy[, 1],
         receiver_y = rxy[, 2]) |> 
  select(receiver_station, 
         receiver_id, receiver_sn, receiver_start, receiver_end, 
         receiver_int, receiver_x, receiver_y) |>
  as.data.frame()
moorings <- copy(moorings_raw)

#### Check deployment periods
ggplot(moorings) +
  geom_segment(aes(
    x    = receiver_start,
    xend = receiver_end,
    y    = factor(receiver_id),
    yend = factor(receiver_id)
  ), size = 2)

#### Check map by deployment period
# (optional) TO DO - see above.

#### Add detection probability parameters
# This is implemented later


###########################
###########################
#### Prepare detections

#### Examine selected columns
nrow(detections)
head(detections)
table(detections$passFilter)
hist(detections$length)
max(detections$length)

#### Process detections
# (Temporarily retain receiver station labels)
detections <- 
  detections |> 
  mutate(individual_id = as.integer(as.character(animal_id)), 
         timestamp = as.POSIXct(detection_timestamp_utc, tz = "UTC"),
         receiver_id = NA_integer_, 
         receiver_sn = as.integer(as.character(receiver_sn))) |>
  select(individual_id, 
         timestamp,
         receiver_sn,
         receiver_station = StationName) |> 
  as.data.table()

#### Order detections by duration
# This improves speed during algorithm testing 
# (NB: the code below works because duration is unique to each individual)
detections <- 
  detections |>
  group_by(individual_id) |> 
  mutate(duration = as.numeric(difftime(max(timestamp), min(timestamp), units = "days"))) |> 
  arrange(duration, timestamp) |> 
  mutate(-duration) |>
  as.data.table()

#### Define receiver_id (~4 s)
# Match using receiver_sn and the time stamps
tic()
for (i in 1:nrow(moorings)) {
  detections[receiver_sn == moorings$receiver_sn[i] & 
               timestamp %within% moorings$receiver_int[i], receiver_id := moorings$receiver_id[i]]
}
toc()
table(is.na(detections$receiver_id))

#### Drop 'extra' detections
detections <- detections[!is.na(receiver_id), ]


###########################
###########################
#### Clean up

# Define study period
study_start <- min(detections$timestamp)
study_end   <- max(detections$timestamp)
study_int   <- lubridate::interval(study_start, study_end)

#### Clean up moorings 
head(moorings)
nrow(moorings)
moorings <- 
  moorings |> 
  as.data.frame() |>
  mutate(int = lubridate::interval(receiver_start, receiver_end)) |> 
  filter(int_overlaps(int, study_int)) |> 
  select(receiver_station, receiver_id, receiver_start, receiver_end, receiver_x, receiver_y) |> 
  as.data.table()
nrow(moorings)
# Define moorings_stations
moorings_stations <- 
  moorings |> 
  group_by(receiver_station) |> 
  mutate(receiver_x = mean(receiver_x), 
         receiver_y = mean(receiver_y)) |> 
  slice(1L) |> 
  ungroup() |> 
  mutate(receiver_id = row_number()) |> 
  as.data.table()
# Define moorings for simulation analyses
# * We average the locations of the receivers in each Station
# * (Receivers were redeployed in the same area (station) after servicing)
# * The receiver_start and receiver_end columns will be replaced later
#   in line with the simulation timeline (see sim-data.R)
moorings_sim <- 
  moorings_stations |> 
  select(-receiver_station) |> 
  as.data.table()
nrow(moorings_sim)
# Define moorings for real-world analyses
moorings_real <- 
  moorings |> 
  select(-receiver_station) |> 
  as.data.table()
rm(moorings)
# Validation: receiver_ids
stopifnot(length(unique(moorings_sim$receiver_id)) == 31L)
stopifnot(length(unique(moorings_real$receiver_id)) == 153L)
# Validation: positions
pp <- par(mfrow = c(1, 2))
terra::plot(map)
points(moorings_sim$receiver_x, moorings_sim$receiver_y, pch = ".")
terra::plot(map)
points(moorings_real$receiver_x, moorings_real$receiver_y, pch = ".")
par(pp)

#### Clean up detections
# Select columns
head(detections)
detections <-
  detections |> 
  select(individual_id, timestamp, receiver_id, receiver_station) |> 
  as.data.table()
# Record 'raw' detections
detections_raw <- copy(detections)

#### Checks
# number of receivers with detections
(nr <- length(unique(detections$receiver_id))) # 137
stopifnot(nr > 1L)


###########################
###########################
#### Visualise mooring stations

# It is important to visualise mooring stations
# Below we will check for evidence of movement through receiver gates without detection
# This is important to understand for the parameterisation of the acoustic observation model

#### Validate station names
if (FALSE) { 
  
  # Define detection_stations data.table
  detections_stations <- 
    lazy_dt(detections) |> 
    left_join(moorings_raw |> 
                select(receiver_id, receiver_station, receiver_x, receiver_y), 
              by = "receiver_id") |> 
    mutate(receiver_station.x = as.character(receiver_station.x), 
           receiver_station.y = as.character(receiver_station.y)) |> 
    as.data.table()
  
  # Confirm that receiver_stations are valid 
  stopifnot(all(detections_stations$receiver_station.x == detections_stations$receiver_station.y))
  
  # Map receiver stations
  # > This is to check the spatial distribution of receiver stations has been correctly assigned
  if (FALSE) {
    
    # Define base maps
    stations <- sort(unique(moorings_raw$receiver_station))
    champlain_utm_facets <- bind_rows(lapply(stations, function(s) {
      mutate(champlain_utm, receiver_station = s)
    }))
    
    # Map receiver stations in moorings_raw
    tic()
    ggplot() +
      geom_sf(data = champlain_utm_facets) +
      geom_point(data = moorings_raw, aes(receiver_x, receiver_y)) +
      facet_wrap(~receiver_station)
    toc()
    
    # Plot receiver stations in detections (~15 mins!)
    tic()
    png(here_fig("stations-detections.png"), 
        height = 20, width = 20, units = "in", res = 600)
    p <- 
      ggplot() +
      geom_sf(data = champlain_utm_facets) +
      geom_point(data = 
                   detections_stations |> 
                   # Slice by receiver_x and receiver_y for speed
                   group_by(receiver_x, receiver_y) |> 
                   slice(1L) |> 
                   as.data.table(), 
                 aes(receiver_x, receiver_y)) +
      facet_wrap(~receiver_station)
    print(p)
    dev.off()
    toc()
  }
  
}

#### Define map layers (WGS84)
# Define map_ll
map_ll <- terra::project(map, "EPSG:4326")
# Define moorings_real_ll
moorings_real_ll <- 
  moorings_real |> 
  sf::st_as_sf(coords = c("receiver_x", "receiver_y"),
               crs = epsg_utm) |> 
  sf::st_transform(4326)
# Define stations_ll
stations_ll <- moorings_stations |> 
  sf::st_as_sf(
    coords = c("receiver_x", "receiver_y"),
    crs = epsg_utm) |> 
  sf::st_transform(4326)
# Define containers_ll
containers_ll <- 
  moorings_real |> 
  sf::st_as_sf(coords = c("receiver_x", "receiver_y"),
               crs = epsg_utm) |> 
  sf::st_buffer(dist = 1000) |> 
  sf::st_transform(4326)

#### Map receiver stations (static)
terra::plot(map)
points(moorings_real$receiver_x, moorings_real$receiver_y)
basicPlotteR::addTextLabels(moorings_stations$receiver_x, 
                            moorings_stations$receiver_y, 
                            moorings_stations$receiver_station)

#### Map receiver stations (interactive) 
leaflet() |>
  addProviderTiles(providers$Esri.WorldImagery) |>
  addRasterImage(map_ll, opacity = 0.7) |>
  addPolygons(
    data = containers_ll,
    fillOpacity = 0.2,
    weight = 1) |>
  addCircleMarkers(
    data = moorings_real_ll,
    radius = 4,
    stroke = FALSE,
    fillOpacity = 1) |>
  addLabelOnlyMarkers(
    data = stations_ll,
    label = ~receiver_station,
    labelOptions = labelOptions(noHide = TRUE, direction = "top", textOnly = TRUE))

#### Station structure
# Whallon/Split Rock (northern end of southern area of Lake)
# Arnold West/Arnold Central/Arnold East receiver gate (further south)
# Crown Point (southernmost receiver)


###########################
###########################
#### Apply filters 

#### Raw data summary statistics
# 1,735,137 detections
# from 93 individuals
# derived from 153 receiver deployments 
# over a four year period

nrow(detections)
length(unique(detections$individual_id))
range(detections$timestamp)
difftime(max(detections$timestamp), min(detections$timestamp), units = "days")
nrow(moorings_real)

#### False detections 
# Check study duration
c(study_start, study_end)
range(detections_filtered$detection_timestamp_utc)
# Format detections_filtered
detections_filtered |> setDT()
detections_filtered[, individual_id := as.integer(as.character(animal_id))]
# Filter detections using detections_filtered (inner join)
# * This ensures we carry forward false detection filters etc. implemented by Futia et al. (2024)
detections_pre_filter <- copy(detections)
detections <- detections[
  detections_filtered,
  on = .(individual_id = individual_id, timestamp = detection_timestamp_utc),
  nomatch = 0
]
detections <- detections[, .(individual_id, timestamp, receiver_id, receiver_station)]

#### Summarise filtered detection dataset
nrow(detections_pre_filter)
nrow(detections_filtered)
nrow(detections)
length(unique(detections$individual_id))
length(unique(detections_filtered$individual_id))


###########################
###########################
#### Examine detection transitions

#### Potential movements of interest 
# Movements through detection gates can be challenging depending on the formulation of the detection probability model
# Specific movements we should check for, based on examination of the station layout, include: 
# * Movements north from Crown Point to any other receiver without detection at Arnold West/Arnold Central/Arnold East
# * Movements from any receiver to Crown Point without detection at Arnold West/Arnold Central/Arnold East
# * Movements to/from Arnold West/Arnold Central/Arnold East without detection at Crown Point (from the South) or Split Rock (from the North)
# (There are others but based on initial modelling movements in the south appear more problematic) 

# Define dataset
# detections_dataset <- copy(detections_pre_filter)
detections_dataset <- copy(detections)

# Define detection transitions from receiver_station -> receiver_station_next
detections_transitions <- 
  detections_dataset |> 
  mutate(time_id = lubridate::floor_date(timestamp, "months")) |> 
  group_by(individual_id, time_id) |> 
  mutate(receiver_station_next = lead(receiver_station)) |>
  ungroup() |>
  filter(!is.na(receiver_station_next)) |> 
  select(individual_id, timestamp, time_id, receiver_station, receiver_station_next) |> 
  as.data.table()

# There are no detections from Crown Point further north to any receiver without detection at Arnold
# (Only one individual detected at Crown Point)
detections_transitions |>
  filter(receiver_station %in% "Crown Point")

# There are no detections arriving at Crown Point without detection at Arnold
detections_transitions |> 
  filter(receiver_station_next %in% "Crown Point")

# There are some transitions from Crown Point -> Whallon (may be acceptable for model?)
# There are also more problematic transitions to:
# * Saxton
# * Willsboro
# * Burlington
# * Winooski Delta Mid
# * Schuyler

detections_transitions |> 
  filter((receiver_station %in% c("Arnold West", "Arnold Central", "Arnold East") & 
            !(receiver_station_next %in% c("Arnold West", "Arnold Central", "Arnold East", "Crown Point", "Split Rock", "Whallon"))))

# individual_id           timestamp    time_id receiver_station receiver_station_next
# <int>              <POSc>     <POSc>           <fctr>                <fctr>
# 1:         24391 2015-05-11 03:10:15 2015-05-01      Arnold West              Schuyler
# 2:         24327 2015-12-19 12:49:48 2015-12-01      Arnold East    Winooski Delta Mid
# 3:         24385 2016-07-03 18:09:55 2016-07-01      Arnold West             Willsboro
# 4:         24339 2016-12-02 08:53:33 2016-12-01      Arnold East                Saxton
# 5:         24385 2017-05-23 17:36:23 2017-05-01   Arnold Central            Burlington
# 6:         24321 2017-06-07 11:11:40 2017-06-01   Arnold Central              Schuyler

# There are some movements from Whallon -> Arnold without detection on Split Rock which may be problematic
detections_transitions |> 
  filter((receiver_station %in% c("Whallon") & 
            (receiver_station_next %in% c("Arnold West", "Arnold Central", "Arnold East")))) |> 
  arrange(individual_id, time_id)

#### Drop receiver_station
detections[, receiver_station := NULL]


###########################
###########################
#### Write outputs

qs::qsave(fish, here_input("fish.qs"))
qs::qsave(moorings_sim, here_input_sim("main", "moorings-xy.qs"))
qs::qsave(moorings_real, here_input_real("main", "moorings.qs"))
qs::qsave(detections, here_input_real("main", "detections.qs"))
qs::qsave(detections_raw, here_input_real("main", "detections-raw.qs"))


#### End of code. 
###########################
###########################