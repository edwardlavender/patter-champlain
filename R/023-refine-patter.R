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
library(patter)
library(proj.verse)
library(tictoc)
files_source_r(here_src())

#### Load data
map <- terra::rast(here_input("map.tif"))


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
set_seed()
set_map(map)

#### Select individuals (sim)
# For convergence failures, see analysis-patter.R
# (That code must be run on machine where output files live)
# it <- iteration[individual_id == 26 & sensitivity == "best", ] # simulation 

#### Select individuals (real)
# it <- iteration[individual_id == 24352 & time_id == as.POSIXct("2017-05-01 00:00:00") & sensitivity == "best", ]
# it <- iteration[individual_id == 24329 & time_id == as.POSIXct("2017-03-01 00:00:00") & sensitivity == "best", ]
# it <- iteration[individual_id == 24331 & time_id == as.POSIXct("2016-03-01 00:00:00") & sensitivity == "best", ] 
# it <- iteration[individual_id == 24333 & time_id == as.POSIXct("2016-02-01 00:00:00") & sensitivity == "best", ] 
# it <- iteration[individual_id == 24352 & time_id == as.POSIXct("2015-03-01 00:00:00") & sensitivity == "best", ] 
it <- iteration[individual_id == 24321 & time_id == as.POSIXct("2016-10-01 00:00:00") & sensitivity == "best", ]; it$index
# it <- iteration[individual_id == 24352 & time_id == as.POSIXct("2016-05-01 00:00:00") & sensitivity == "best", ]; it$index
# it <- iteration[individual_id == 24352 & time_id == as.POSIXct("2016-06-01 00:00:00") & sensitivity == "best", ]; it$index
# it <- iteration[individual_id == 24370 & time_id == as.POSIXct("2015-05-01 00:00:00") & sensitivity == "best", ]; it$index 
# it <- iteration[individual_id == 24387 & time_id == as.POSIXct("2015-04-01 00:00:00") & sensitivity == "best", ]; it$index
# it <- iteration[individual_id == 26805 & time_id == as.POSIXct("2016-04-01 00:00:00") & sensitivity == "best", ]; it$index

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
yobs <- list(ModelObsAcousticLogisTrunc = copy(acoustics), 
             ModelObsContainer = NULL)

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
              .state      = state,
              .model_move = model_move,
              .yobs       = yobs,
              .n_move     = 1000L,
              .n_particle = 50000L,
              .n_resample = 10000,
              .n_record   = it$n_particle_smoother,
              .direction  = direction)
# Run filter
pout <- do.call(pf_filter, pargs, quote = TRUE)

#### Record
# see debug-by-individual.txt


###########################
###########################
#### Animation

# Create animation
# * For 12,000 steps: This takes 6 min (12 cl) plus >2 min for 12,000 steps
# * NB: Running this for a few steps with .cl = 1L seems to suppress a segmentation
#   fault when it then run in parallel for a larger time series below.
#   If you jump to the parallel version, it can throw a segmentation fault
steps <- 1:11054
tnow <- as.numeric(Sys.time())
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