###########################
###########################
#### analysis-optimisation-turning-angle.R

#### Aims
# 1) Interactive simulation analysis
# - Easily simulate a trajectory and observations with user-defined parameters
# - Examine estimation of latent locations and static parameters (turning angle)
# - We test specifically whether we can estimate turning angle SD when wide/narrow
#   in both the real-world array and a dense simulated array
# - This code is written for quick & easy application & adjustment

#### Prerequisites
# 1) This code was modified from the following sources:
#    * sim-data.R
#    * ?pf_filter


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Set global options
patter::julia_connect()

#### Load essential packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(JuliaCall)
library(glue)
library(patter)
library(prettyGraphics)
library(proj.verse)
library(tictoc)
files_source_r(here_src())
expect_no_geospatial()

#### Load data
fish              <- qs::qread(here_input("fish.qs"))
moorings          <- qs::qread(here_input_sim("moorings-xy.qs"))
# pars_model_move <- qs::qread(here_input("pars-model-move-best.qs"))
pars_model_obs    <- qs::qread(here_input("pars-model-obs-best.qs"))
map_bbox          <- qs::qread(here_input("map-bbox.qs"))
if (!patter:::os_linux()) {
  map <- terra::rast(here_input("map.tif"))
} else {
  map <- NULL
}


###########################
###########################
#### Simulate data 

#### Setup Julia
seed <- 123L
set_seed(seed)
set_map(here_input("map.tif"))

#### Define n_sim (1L)
n_sim <- 1L

#### Define timeline (one-month)
timeline <- seq(as.POSIXct("2025-01-01 00:00:00", tz = "UTC"), 
                as.POSIXct("2025-01-31 23:58:00", tz = "UTC"), 
                by = "2 mins")

#### Define movement model
state           <- state_trout()
pars_model_move <- list(shape = 3.25, scale = 25, mobility = 216, phi = 1.8)
model_move      <- model_move_trout(pars_model_move)
plot(model_move)

#### Define tagging location
# (Extract map_value via Patter for linux handling)
xinit <- fish[sample.int(n_sim, replace = TRUE), ]
julia_assign("x0", xinit$x)
julia_assign("y0", xinit$y)
xinit[, map_value := julia_eval('[Patter.extract(env, x0[i], y0[i]) for i in eachindex(x0)]')]
xinit <- model_move_xinit(.xinit = xinit, .n_particle = NULL)
stopifnot(all(xinit$map_value == 1L))

#### Simulate movement path
# phi = 0.4 is easy to generate, phi = 0.3 is hard
tic()
count <- counter <- 1
while (counter < 120) {
  path <- tryCatch(sim_path_walk(.map        = map, 
                                 .timeline   = timeline, 
                                 .state      = state, 
                                 .model_move = model_move,
                                 .xinit      = xinit), 
                   error = function(e) e)
  if (inherits(path, "error")) {
    count <- counter <- counter + 1
  } else {
    counter <- Inf
  }
}
count
toc()

#### Define observation model
# analysis_moorings <- "real"
analysis_moorings <- "sim"
if (analysis_moorings == "real") {
  
  # Use real-world array 
  stopifnot(nrow(moorings) == 31L)
  moorings[, receiver_start := min(timeline) -  24 * 60 * 60]
  moorings[, receiver_end := max(timeline) + 24 * 60 * 60]
  model_obs <- model_obs_champlain(moorings, pars_model_obs)
  
} else if (analysis_moorings == "sim") {
  
  #### Use dense simulated array (sense check)
  # This code must be run on MacOS/Windows
  moorings <- sim_array(.map = map, 
                        .timeline = timeline, 
                        .arrangement = "regular", 
                        .n_receiver = 1000L)
  
  #### Crop moorings 
  # For speed, we focus on the receivers in a buffer around the simulated path
  # (We do not need to consider non-detections from receivers very far away)
  region <- 
    terra::ext(c(range(path$x), range(path$y))) |>
    terra::vect(crs = terra::crs(map)) |> 
    terra::buffer(width = 20000)
  # Visually check region size
  terra::plot(map)
  patter:::add_sp_path(path$x, path$y, length = 0)
  terra::lines(region)
  # Filter moorings
  moorings[, map_value := terra::extract(region, cbind(moorings$receiver_x, moorings$receiver_y))[, 2]]
  moorings[!is.na(map_value), ]
  nrow(moorings)
  
  ##### Customise detection probability parameters
  # Check distances between receivers
  dist <- terra::distance(cbind(moorings$receiver_x, moorings$receiver_y), lonlat = FALSE)
  min(dist[dist > 0])
  # Customise parameters
  # plot(1:1000, plogis(4 + -0.01 * 1:1000), type = "l")
  # pars_model_obs <- list(receiver_alpha = 4, receiver_beta = -0.01, receiver_gamma = 800)
  
}
model_obs <- model_obs_champlain(moorings, pars_model_obs)
nrow(model_obs$ModelObsAcousticLogisTrunc)
plot(model_obs)

#### Simulate observations
# Simulate acoustic observations 
# * ~10 s (real-world array)
# * ~30 s (simulated array with 294 receivers)
tic()
obs        <- sim_observations(.timeline = timeline, 
                               .model_obs = model_obs)
acoustics  <- obs$ModelObsAcousticLogisTrunc[[1]]
toc()
# Focus on the period of detections (for speed)
detections <- acoustics[obs == 1L, ]
stopifnot(nrow(detections) > 100L)
timeline <- timeline[timeline >= min(detections$timestamp) & 
                       timeline <= max(detections$timestamp)]
acoustics  <- acoustics[timestamp %in% timeline, ]
# Build containers
containers <- assemble_acoustics_containers(.timeline  = timeline,
                                            .acoustics = acoustics,
                                            .mobility  = pars_model_move$mobility,
                                            .map       = map_bbox)
# Collate observations for forward filter
yobs_fwd <-
  list(ModelObsAcousticLogisTrunc = acoustics,
       ModelObsContainer = containers$forward)

#### Check observations
# High phi leads to good detection time series
range(timeline)
plot(detections$timestamp, detections$sensor_id)


###########################
###########################
#### Prepare optimisation

# Define a function that computes the log-likelihood given input parameters
# * `theta` denotes a parameter/parameter vector
# * (In this case, it is standard deviation of the turning angle distribution)
# * This function is defined locally here for convenience 
# * (We use variables in the global environment)

pf_filter_ll <- function(theta) {
  
  # Safety checks
  if (theta <= 0) {
    return(-Inf)
  }
  
  # Instantiate movement model
  pmm        <- pars_model_move
  pmm$phi    <- theta
  model_move <- model_move_trout(pmm)
  
  # Run filter
  # * (optional) TO DO Define t_resample
  # * Filter settings are defined to minimise computation time over multiple runs
  fwd <- pf_filter(.timeline   = timeline,
                   .state      = state,
                   .xinit      = NULL,
                   .model_move = model_move,
                   .yobs       = yobs_fwd,
                   .n_particle = 2e4L,
                   .n_resample = 500L,
                   .n_move     = 1L,
                   .n_record   = 1L,
                   .direction  = "forward",
                   .progress   = julia_progress(enabled = TRUE),
                   .verbose    = FALSE)
  
  # Return log-lik
  fwd$callstats$loglik
  
}


###########################
###########################
#### Run optimisation

#### Define parameter grid
# Define thetas 
# (We ensure that the grid contains the true parameter value)
if (pars_model_move$phi == 1.8) {
  thetas <- unique(sort(c(seq(0.1, 3.0, by = 0.5), pars_model_move$phi)))
} else if (pars_model_move$mobility == 0.4) {
  thetas <- unique(sort(c(seq(0.1, 1, by = 0.1), seq(1, 3, 0.5), pars_model_move$phi)))
} else {
  stop("A custom theta grid is required for other values of pars_model_move$phi.")
}
# Define iteration grid 
loglik <- 
  CJ(rep_id = 1:3, 
     theta = thetas) |>
  mutate(index = row_number()) |>
  select(index, rep_id, theta) |> 
  as.data.table()
# Check rows (computation time)
nrow(loglik)

#### Run grid search 
# With the real-world array: 
# ~1 min per run (10 threads, SIA-LAVENDED-M)
# ~21-31 mins (phi = c(1.8, 0.4))
# With the simulated array (294 receivers):
# ~4 mins per run 
values <- cl_lapply(split(loglik, loglik$index), function(d) {
  pf_filter_ll(d$theta)
}) |> unlist()
loglik[, value := values]


###########################
###########################
#### Analyse results

# Plot likelihood profile
png(here_fig("model-move", "optim", glue("optim-phi-{analysis_moorings}-{pars_model_move$phi}-{seed}.png")), 
    height = 5, width = 5, units = "in", res = 600)
pp <- par(oma = c(2, 2, 1, 1))
ylim <- range(loglik$value[is.finite(loglik$value)])
pretty_plot(loglik$theta, loglik$value, 
            ylim = ylim,
            xlab = "", ylab = "")
rug(loglik$theta,  ticksize = 0.02, pos = ylim[1], lwd = 1.25)
abline(v = pars_model_move, lty = 3)
mtext(side = 1, "Turning angle standard deviation", line = 2)
mtext(side = 2, "Log likelihood", line = 5)
dev.off()


###########################
###########################
#### Record outputs

# (optional) TO DO Record outputs


#### End of code
###########################
###########################