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

#### Select individual
# Read callstats 
# (This code is modified from analysis-patter.R)
if (analysis == "sim") {
  iteration[, file_convergence := file_callstats]
} else if (analysis == "real") {
  iteration[, file_convergence := file_callstats_filter]
}
table(file.exists(iteration$file_convergence))
iteration <- iteration[file.exists(file_convergence), ]
callstats <- lapply(iteration$index, function(i) {
  iteration$file_convergence[iteration$index == i] |> 
    arrow::read_feather() |> 
    mutate(index = i, .before = 1L) |> 
    cbind(iteration[index == i, .(individual_id, time_id, sensitivity, sensitivity_label)]) |> 
    as.data.table()
}) |> 
  rbindlist()
# Identify convergence failures
failures <- 
  callstats |> 
  filter(routine == "filter: forward") |> 
  filter(convergence == FALSE) |> 
  as.data.table()
# Summarise failures
failures |> 
  group_by(sensitivity) |>
  summarise(n())
# Visualise failures
file_convergence

#### Select individual & read data 
# it <- iteration[individual_id == 26 & sensitivity == "best", ]
timeline  <- arrow::read_feather(it$file_timeline)
timeline  <- timeline$timestamp
acoustics <- arrow::read_feather(it$file_acoustics)
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
direction <- "backward"
stopifnot(direction %in% c("forward", "backward"))
if (direction == "forward") {
  yobs$ModelObsContainer <- copy(containers_fwd)
} else {
  yobs$ModelObsContainer <- copy(containers_bwd)
}
# Define arguments
pargs <- list(.timeline   = timeline,
              .state      = state,
              .model_move = model_move,
              .yobs       = yobs,
              .n_move     = 1L,
              .n_particle = it$n_particle_filter,
              .n_record   = it$n_particle_smoother,
              .direction  = direction)
# Run filter
pout <- do.call(pf_filter, pargs, quote = TRUE)

#### Simulation record
## index 176, individual 26, sensitivity == "best"
# - default settings, .n_move = 1L, .n_particle = 50000
# - SIA-LAVENDED: success
# - siam-linux20: failure (assume unlucky)

#### Real-world record (100 test time series)
# TO DO 


###########################
###########################
#### Animation

# Create animation
# * For 12,000 steps: This takes 6 min (12 cl) plus >2 min for 12,000 steps
# * NB: Running this for a few steps with .cl = 1L seems to suppress a segmentation
#   fault when it then run in parallel for a larger time series below.
#   If you jump to the parallel version, it can throw a segmentation fault
steps <- 1:1000L
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