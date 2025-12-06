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
library(sf)
library(rnaturalearth)
library(proj.verse)
library(tictoc)
files_source_r(here_src())

#### Load data 
champlain  <- terra::vect(here_data_raw_mf("ChamplainRegionsGrouped/ChamplainRegionsGrouped.shp"))


###########################
###########################
#### Define study area (~2 s)

#### Build map 
# Use a coarse map for speed sampling initial locations 
tic()
epsg_utm         <- "EPSG:3175"
champlain$land   <- as.numeric(1)
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
terra::lines(champlain |> terra::project(epsg_utm))
# Zoom-in to check resolution
map_zoom <- terra::crop(map, 
                        cbind(1787707, 977397.3) |>
                          terra::vect() |>
                          terra::buffer(width = 10000) |>
                          terra::ext())
terra::plot(map_zoom)
terra::lines(champlain |> terra::project(epsg_utm))
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

#### Get a high-res North America polgyon
continent <- ne_download(
  scale = 10,
  type = "admin_0_countries",
  category = "cultural",
  returnclass = "sf"
)
continent <- continent[continent$CONTINENT == "North America", ]

#### (optional) Get North American lakes 
lakes <- ne_download(scale = 10,
                     type = "lakes",
                     category = "physical",
                     returnclass = "sf"
) |> 
  st_make_valid()
# lakes <- sf::st_crop(lakes, continent)

#### (optional) Get North American rivers 
rivers <- ne_download(
  scale = 10,
  type = "rivers_lake_centerlines",
  category = "physical",
  returnclass = "sf"
)
# rivers <- sf::st_crop(rivers, continent)

#### Collect geometries
names(champlain_ll$region)
continent <- st_geometry(continent)
lakes     <- st_geometry(lakes)
rivers    <- st_geometry(rivers)

#### Plot layers 
plot(continent, col = "gray90")
plot(lakes, col = "lightblue", border = NA, add = TRUE)
plot(rivers, col = "lightblue", lwd = 0.1, add = TRUE)
# plot(st_geometry(champlain_ll), col = "blue", add = TRUE) # slow

#### Write layers
st_write(continent, here_fig("local", "qgis", "continent.shp"), append = FALSE)
st_write(champlain_ll, here_fig("local", "qgis", "champlain.shp"), append = FALSE)


#### End of code. 
###########################
###########################