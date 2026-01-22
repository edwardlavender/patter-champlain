###########################
###########################
#### refine-patter.R

#### Aims
# 1) This script is used to trial/debug the particle filter i.e., convergence failures
#    - For a subset of time series, we run the particle filter via run-filter.jl
#    - We identify cases of convergence failures
#    - In this script, we re-run the filter for those individuals & analyse behaviour of particles
#    - We then revise model parameters/algorithm settings accordingly in an iterative process

#### Prerequisites
# 1) Run run-filter.jl to assess convergence quickly via 001-workflow.sh
# 2) Set JULIA_NUM_THREADS = 12 for this script
# 3) This script is designed to run on SIA-LAVENDED


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
library(JuliaCall)
library(patter)
library(proj.verse)
library(tictoc)
files_source_r(here_src())

#### Load data
map      <- terra::rast(here_input("map.tif"))
regions  <- terra::rast(here_input("regions.tif"))
moorings <- qs::qread(here_input_real("main", "moorings.qs"))


###########################
###########################
#### Select analysis

#### Define analysis 
# analysis <- "sim"
analysis <- "real"
subanalysis <- "main"

#### Define analysis-specific data
here_input_analysis <- switch_here_input_analysis_subanalysis(analysis, subanalysis)
here_fig_analysis   <- switch_here_fig_analysis_subanalysis(analysis, subanalysis)
iteration           <- qs::qread(here_input_analysis("iteration.qs"))
if (analysis == "real") {
  message("Using iteration-1.qs")
  iteration <- qs::qread(here_input_analysis("iteration-1.qs"))
}


###########################
###########################
#### Run filter

#### Connect to Julia
julia_connect()
julia_source(file.path("Julia", "src", "observation-model.jl"))
set_seed()
set_map(map)

#### Examine map
terra::plot(regions)
points(moorings$receiver_x, moorings$receiver_y)
terra::sbar(2000)
# flapper::dist_btw_clicks(lonlat = FALSE)

#### Select individuals (sim)
# For convergence failures, see analysis-patter.R
# (That code must be run on machine where output files live)
# it <- iteration[individual_id == 26 & sensitivity == "best", ] # simulation 

#### Select individuals (real)
it <- iteration[individual_id == 24385 & time_id == as.POSIXct("2017-05-01 00:00:00") & sensitivity == "best", ]; it$index
it$n_particle_filter <- 20000L

#### Read individual-specific data
timeline       <- arrow::read_feather(it$file_timeline)
timeline       <- timeline$timestamp
acoustics      <- arrow::read_feather(it$file_acoustics)
detections     <- acoustics[obs == 1L, ]
containers_fwd <- arrow::read_feather(it$file_containers_fwd)
containers_bwd <- arrow::read_feather(it$file_containers_bwd)

#### Define movement model
state      <- "StateCXY"
model_move <- model_move_trout(it)
model_move

#### Define observation(s) & observation model
# (Containers are added below)
if (TRUE) {
  dist <- 1:7000
  plot(dist, plogis(1.205216 - 0.001672085 * dist), type = "l", col = "green") # coef(models[["F151"]][["GLM"]])
  lines(dist, plogis(0.635159329 - 0.002724118 * dist), col = "darkred")       # coef(models[["F146"]][["GLM"]])
  lines(dist, plogis(0.9039122 - 0.002090107 * dist), col = "red")             # formerly 'restrictive' model between F146 & F151
  lines(dist, plogis(0.9201876645 - 0.0021981015 * dist))                      # midpoint between F141 and F146
  lines(dist, plogis(1.205216 * 0.25 +  0.635159329 * 0.75 
                     - (0.001672085 * 0.25 + 0.002724118 * 0.75) * dist))
  acoustics[, receiver_alpha := 1.205216 * 0.25 +  0.635159329 * 0.75]
  acoustics[, receiver_beta := -(0.001672085 * 0.25 + 0.002724118 * 0.75)]
  # acoustics[, receiver_alpha := 0.9039122]
  # acoustics[, receiver_beta := -0.002090107]
}
yobs <- list(ModelObsAcousticLogisTruncLos = copy(acoustics), 
             ModelObsContainer = NULL)

#### Define starting locations for forward filter)
# ~03:21, 24325, 2016-03-01, best
pinit <- pf_filter_fwd_xinit(iter       = it, 
                             map        = map, 
                             timeline   = timeline, 
                             acoustics  = acoustics, 
                             model_move = model_move)
xinit <- pinit$xinit
terra::plot(map)
points(xinit$x, xinit$y, pch = ".")

#### Run filter 
# Define direction & update yobs 
direction <- "forward"
stopifnot(direction %in% c("forward", "backward"))
if (direction == "forward") {
  yobs$ModelObsContainer <- copy(containers_fwd)
} else {
  yobs$ModelObsContainer <- copy(containers_bwd)
}
# Compute duration before first detection
if (direction == "forward") {
  difftime(min(detections$timestamp), min(timeline))
} else {
  difftime(max(detections$timestamp), max(timeline))
}
# Define arguments
pargs <- list(.timeline   = timeline,
              .xinit      = xinit,
              .state      = state,
              .model_move = model_move,
              .yobs       = yobs,
              .n_move     = it$n_move,
              .n_particle = it$n_particle_filter,
              .n_resample = it$n_resample, 
              .t_resample = sort(unique(which(timeline %in% yobs$ModelObsContainer$timestamp))),
              .n_record   = it$n_particle_smoother,
              .direction  = direction)
# Run filter
pout <- do.call(pf_filter, pargs, quote = TRUE)

#### Timings
# 24325, 2016-03-01, best:
# 0:16:45, default initialisation & LoS
# 0:21:32, new initialisation & LoS @ moment of detection only 

#### Record
# see debug-by-individual.txt


###########################
###########################
#### Animation

#### Segmentation

# NB: Running this for a few steps with .cl = 1L seems to suppress a segmentation
#   fault when it then run in parallel for a larger time series below.
#   If you jump to the parallel version, it can throw a segmentation fault

#### Define steps
# Define focal region
start <- 15000
focal <- 18000:20533 # 11920 
# Define steps, using low resolution before focal region for speed
# (while including all relevant detection container time steps)
steps <- sort(unique(c(seq(start, min(focal) - 1, by = 10), 
                       focal, 
                       # which(timeline %in% yobs$ModelObsContainer$timestamp),
                       which(timeline %in% detections$timestamp)
                     )))
steps <- steps[steps > 0 & steps <= max(focal)]
length(steps)
tnow <- as.numeric(Sys.time())

#### Make animation
# For 12,000 steps: This takes 6 min (12 cl) plus >2 min for 12,000 steps
animate_ac(.iter   = it,
           .map    = map,
           .steps  = steps[1:10L],
           .input  = pargs,
           .output = pout,
           .outdir = here_fig_analysis("debug", it$individual_id, it$time_id, tnow),
           .cl     = 1L)
if (length(steps) > 10L) {
  animate_ac(.iter   = it,
             .map    = map,
             .steps  = steps,
             .input  = pargs,
             .output = pout,
             .outdir = here_fig_analysis("debug", it$individual_id, it$time_id, tnow),
             .cl     = 10L)
}


#### End of code. 
###########################
###########################
