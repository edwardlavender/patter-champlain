###########################
###########################
#### develop-model-move-flapper.R

#### Aims
# 1) Develop a movement model for lake trout with flapper analyses

#### Prerequisites
# 1) Process detection data


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
library(tictoc)
files_source_r(here_src())

#### Load data
map         <- terra::rast(here_input("map.tif"))
detections  <- qs::qread(here_input_real("detections.qs"))
moorings    <- qs::qread(here_input_real("moorings.qs"))


###########################
###########################
#### Detection analyses of mobility 

#### Estimation method
# Using detection data, we can estimate mobility
# We compute three speed metrics & for each we compute the min/mean/max value
# Metrics:
# * speed_min: minimum speed, based on movement between nearest container edges 
# * speed_avg: average speed, based on movement between receivers
# * speed_max: maximum speed, based on movement between farthest container edge
# Interpretation:
# * min speeds are 'inappropriately low' b/c indirect routes between receivers are taken
# * maximum values for each metric can inform us about mobility 

#### Estimate mobility (~4 s)
if (requireNamespace("flapper", quietly = TRUE)) {
  
  tic()
  # Define LCP grid (resolution in x and y must be exactly identical)
  bb <- terra::ext(map)[1:4]
  bb <- plyr::round_any(bb, 100)
  grid <- terra::rast(terra::ext(bb), res = 100)
  grid <- terra::resample(map, grid)
  # Prepare observations
  detections[, fct := individual_id]
  msp <- sp::SpatialPoints(moorings[, c("receiver_x", "receiver_y")], sp::CRS(terra::crs(map)))
  msp <- sp::SpatialPointsDataFrame(msp, data.frame(receiver_id = moorings$receiver_id))
  # Compute speeds
  # TO DO Repeat with best-guess detection range
  mvt <- flapper::get_mvt_mobility_from_acoustics(data = detections, 
                                                  fct = "individual_id", 
                                                  moorings = msp, 
                                                  detection_range = 7500, 
                                                  calc_distance = "lcp", 
                                                  bathy = raster::raster(grid),
                                                  step = 120,
                                                  transmission_interval = 160)
  toc()
  
}

#### Results 
# These results are broadly robust to detection_range, step and transmission_interval
# --------------------------------------
#   Estimates (m/s)-----------------------
#   variable min mean  max
# 1 speed_min_ms   0 0.03 0.33
# 2 speed_avg_ms   0 0.10 0.90
# 3 speed_max_ms   0 0.18 1.53
# --------------------------------------
#   Estimates (m/step)--------------------
#   variable  min  mean    max
# 1 speed_min_mstep 0.00  3.15  39.44
# 2 speed_avg_mstep 0.08 12.11 108.29
# 3 speed_max_mstep 0.13 21.07 184.03
# --------------------------------------

#### Conclusions 
# Even if we assume 'minimum distance' movements, maximum speeds of 48 m/2 min are apparent.
# (But travelled  speeds are likely to be higher than this value 
# ... as this assumes the minimum possible travel distance)
# If we assume larger movements, speeds may be up to 109 m/2 min.
# This is a reasonable middle-of-the-road estimate for maximum movement speeds. 
# In the most extreme case (unlikely), movement speeds up to 198 m/2 min are apparent.
# (This is the upper bound for movement speeds suggested by the data.)
# These results are consistent for different transmission intervals
# (But the movement model needs to consider the effect of random transmission)


#### End of code. 
###########################
###########################