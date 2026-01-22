###########################
###########################
#### setup-data-map.R

#### Aims
# 1) Sets up the map (used in both simulation & real-world analyses)

#### Prerequisites
# 1) Use ChamplainRegionsGrouped shapefile from M. Futia. 


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
library(dplyr, warn.conflicts = FALSE)
library(leaflet)
library(proj.verse)
library(rnaturalearth)
library(sf)
library(spatial.extensions)
library(tictoc)
files_source_r(here_src())

#### Load data 
champlain  <- terra::vect(here_data_raw_mf("ChamplainRegionsGrouped/ChamplainRegionsGrouped.shp"))
moorings   <- readRDS(here_data_raw_mf("OriginalReceiverSummary_2013-2017.rds"))


###########################
###########################
#### Define study area (~2 s)

#### Visualise study area
# The quality of the champlain shapefile is very high
leaflet() |>
  addProviderTiles(providers$Esri.WorldImagery) |>
  addPolygons(data        = champlain,
              color       = "red",
              weight      = 3,
              fillOpacity = 0)

#### Define UTM SpatVector
# NB: as.numeric(1) is needed for Patter.particle_filter()
epsg_utm         <- "EPSG:3175"
champlain$land   <- as.numeric(1)
champlain_utm    <- champlain |> terra::project(epsg_utm)

#### Build map 
# Use a coarse map for speed sampling initial locations 
tic()
regions <- as_SpatRaster(champlain, .simplify = 0.001, .utm = epsg_utm,
                         .field = "region", .res = 200, .plot = TRUE)
maps    <- as_SpatRaster(champlain, .simplify = 0.001, .utm = epsg_utm,
                         .field = "land", .res = 200, .plot = TRUE)
map     <- maps$SpatRaster
map_len <- terra::ymax(map) - terra::ymin(map)
toc()

#### Examine map properties
# Visualise map simplification
terra::plot(map, col = "blue")
terra::lines(champlain_utm)
# Zoom-in to check resolution
map_zoom <- terra::crop(map, 
                        cbind(1787707, 977397.3) |>
                          terra::vect() |>
                          terra::buffer(width = 10000) |>
                          terra::ext())
terra::plot(map_zoom)
terra::lines(champlain_utm)
# As above, interactively
leaflet() |>
  addProviderTiles(providers$Esri.WorldImagery) |>
  addPolygons(data        = champlain,
              color       = "red",
              weight      = 3,
              fillOpacity = 0) |> 
  addRasterImage(terra::project(map, "WGS84"), opacity = 0.7) 
# Check ncell & compare to dat_gebco() for reference
terra::ncell(map)                  # 107625
terra::ncell(patter::dat_gebco())  # 50160
# Check map size
terra::ncell(map) * 8 / 1e6        # 0.81 MB

#### Map dimensions
# max dimension
terra::ymax(map) - terra::ymin(map) # 175000 m
terra::xmax(map) - terra::xmin(map) # 24600 m
# bbox
map_bbox <- patter:::map_bbox(map)
stopifnot(nrow(map_bbox) == 4L)
terra::plot(map)
points(map_bbox, pch = 21, bg = "red", cex = 5)

#### Write maps
qs::qsave(epsg_utm, here_input("epsg-utm.qs"))
terra::writeRaster(map, here_input("map.tif"), overwrite = TRUE)
terra::writeRaster(regions$SpatRaster, here_input("regions.tif"), overwrite = TRUE)
qs::qsave(map_bbox, here_input("map-bbox.qs"))
file.size(here_input("map.tif")) / 1e6 # MB


###########################
###########################
#### Define QGIS layers

#### Define Lake Champlain (ll)
champlain_ll <- terra::project(champlain, "WGS84")
champlain_ll <- sf::st_as_sf(champlain_ll)

if (curl::has_internet()) {
  
  #### Get a high-res North America polygon
  continent <- ne_download(
    scale = 10,
    type = "admin_0_countries",
    category = "cultural",
    returnclass = "sf"
  )
  continent <- continent[continent$CONTINENT == "North America", ]
  
  #### (optional) Get North American lakes 
  if (FALSE) {
    lakes <- ne_download(scale = 10,
                         type = "lakes",
                         category = "physical",
                         returnclass = "sf"
    ) |> 
      st_make_valid()
    # lakes <- sf::st_crop(lakes, continent)
  }

  
  #### (optional) Get North American rivers 
  if (FALSE) {
    rivers <- ne_download(
      scale = 10,
      type = "rivers_lake_centerlines",
      category = "physical",
      returnclass = "sf"
    )
    # rivers <- sf::st_crop(rivers, continent)
  }

  #### Collect geometries
  names(champlain_ll$region)
  continent <- st_geometry(continent)
  if (FALSE) {
    lakes     <- st_geometry(lakes)
    rivers    <- st_geometry(rivers)
  }

  #### Plot layers 
  plot(continent, col = "gray90")
  if (FALSE) {
    plot(lakes, col = "lightblue", border = NA, add = TRUE)
    plot(rivers, col = "lightblue", lwd = 0.1, add = TRUE)
  }
  # plot(st_geometry(champlain_ll), col = "blue", add = TRUE) # slow
  
  #### Write continent 
  st_write(continent, here_fig("local", "qgis", "continent.shp"), append = FALSE)

  
}

#### Write layers
qsavevect(champlain_utm, here_input("champlain-utm.qs"))
st_write(champlain_ll, here_fig("local", "qgis", "champlain.shp"), append = FALSE)


###########################
###########################
#### Define regions colour scheme

# Define colour scheme 
regions_cs <- tibble::tribble(
  ~region,            ~col,
  "Missisquoi Bay",    "#87e93b", 
  "Northeast Arm",     "#22e45f", 
  "Malletts Bay",      "#df756d", 
  "Main Lake North",   "#9834df", 
  "Main Lake Central", "#6c81de", 
  "Main Lake South",   "#e41ea5", 
  "South Lake",        "#0dcbd9"
) |> 
  mutate(region = factor(region, levels = region)) |>
  as.data.table()

qs::qsave(regions_cs, here_input("regions-colour-scheme.qs"))


###########################
###########################
#### Explore acoustic containers

#### Overview
# To perform inference with patter for trout location, we use acoustic containers
# To optimise the use of containers, we can make them a bit smaller than the default size
# The code below plots the size of the study area to identify a suitable container size

# Define map boundary box
bb <- patter:::map_bbox(map)

# Define map centroid
centroid <- cbind(mean(c(min(bb[, 1]), max(bb[, 1]))), 
                  mean(c(min(bb[, 2]), max(bb[, 2]))))


# Compute max. distance from centroid to edge of map: 118211.8 m
(dist_from_centroid <- max(terra::distance(bb, centroid, lonlat = FALSE)))

# Choose threshold for acoustic containers: 44159.98
(threshold <- dist_from_centroid * 0.5)

# Plot map with centroid & threshold container
terra::plot(map)
centroid |> 
  terra::vect(crs = terra::crs(map)) |> 
  terra::buffer(width = threshold) |> 
  terra::plot(col = scales::alpha("red", 0.25), add = TRUE)
cbind(moorings$deploy_long, moorings$deploy_lat) |> 
  terra::vect(crs = "WGS84") |> 
  terra::project(epsg_utm) |> 
  terra::points()
points(centroid, col = "blue", lwd = 2)


#### End of code. 
###########################
###########################