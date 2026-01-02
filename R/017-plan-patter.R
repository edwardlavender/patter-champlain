###########################
###########################
#### plan-patter.R

#### Aims
# 1) This script is used to plan patter runs
#    * Examine outputs for an example individual
#    * Check file sizes
#    * Define implementation strategy 

#### Prerequisites
# 1) Run simulations


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
library(proj.verse)
library(tictoc)
files_source_r(here_src())
expect_no_geospatial()

#### Load data
map <- terra::rast(here_input("map.tif"))


###########################
###########################
#### Select analysis

#### Define analysis 
analysis <- "sim"
# analysis <- "real"
subanalysis <- "main"

#### Define analysis-specific data
here_input_analysis <- switch_here_input_analysis_subanalysis(analysis, subanalysis)
iteration           <- qs::qread(here_input_analysis("iteration.qs"))


###########################
###########################
#### Generate example outputs

#### Define example individual
it <- iteration[1, ]
list.files(it$folder_output)

#### Run algorithms 
# Run 001-run-algorithms.jl for example individual
# Check file sizes (MB)
if (FALSE) {
  dir_size(it$folder_output)
  file.size(it$file_callstats) / 1e6 
  file.size(it$file_diagnostics) / 1e6
  sum(file.size(list.files(it$folder_output[i], pattern = "\\.jld2$", full.names = TRUE))) / 1e6
}

#### Collate states 
# Run 001-run-algorithms.jl for example individual & check file size for it$file_states (MB)
if (FALSE) {
  file.size(it$file_states) / 1e6
}


###########################
###########################
#### Examine outputs

#### Read data 
it <- iteration[1, ]
it_callstats   <- arrow::read_feather(it$file_callstats)
it_diagnostics <- arrow::read_feather(it$file_diagnostics)
it_states      <- arrow::read_feather(it$file_states)
# Check file sizes
# > We should be able to read all diagnostic files into memory 
nrow(it_diagnostics)
nrow(it_states)

#### Examine callstats
# Simulations: 
# * TO DO s for it[1, ] (SIA-LAVENDED, 50000/1500 particles, 1 thread, 10 batches)
# * TO DO s (ETA) for 210 iterations on 20 cores
# Real-world analysis
# * TO DO for it[1, ] 
# * TO DO s (ETA) for 2723 iterations on 20 cores
it_callstats
sum(it_callstats$time) # secs

#### Examine diagnostics
## Summarise ESS
# For sim it[1, ]:
# * TO DO
# * TO DO
it_diagnostics |> 
  group_by(routine) |> 
  reframe(utils.add::basic_stats(ess, na.rm = TRUE))
# Check convergence
it_callstats$convergence[1]
it_callstats$convergence[2]
it_diagnostics |>
  filter(routine == "smoother: two-filter") |> 
  summarise(percent_valid = length(which(!is.na(ess))) / n() * 100, 
            convergence = percent_valid > 75) |> 
  as.data.table()

#### Map states (~7 s)
# Each map requires < 1 MB of storage
tic()
tmp.tif <- tempfile(fileext = ".tif")
tmp     <- map_pou(map, it_states)$ud
terra::writeRaster(tmp, tmp.tif, overwrite = TRUE)
file.size(tmp.tif) / 1e6
# map_dens(map, it_states, .discretise = TRUE)
toc()


###########################
###########################
#### Define implementation strategy 

#### Define the number of particles for algorithms
# We run algorithms on {n_cpu}
# We record {n_record} particles
# We quadrouple the required memory
# - We may start a new run while still finishing another (x2)
# - We need a safety buffer (x3)
# We use a maximum of {mb_available} memory in total
n_record     <- 1500
n_cpu        <- 25
mb_available <- 25e3
mb_required  <- p_mem(n_record, 22320, 4, n_cpu) * 3
p_batch(mb_required, mb_available)

#### Storage requirements for simulations (1500 particles):
# Breakdown:
# * callstats   : 0.00189 
# * diagnostics : 3.550282
# * jld2 files  : 1071.414
# * states      : 1874.882 (with Arrow.write default options)
# * states      : TO DO    (with Arrow.write & compress = ZstdCompressor(level = 9))
# * map         : 0.068324 (see below)
# * Total       : 1878.502 (excluding jld2 files, which we delete)
# To store the outputs of the simulations, we need ~394 GB
1878.502 * 210 / 1e3 

#### Storage requirements for real-world analysis (2000 particles): 
# Breakdown: 
# * callstats   : 
# * diagnostics : 
# * jld2 files  : 
# * states      : 
# * map         : 
# * total       : 
# To store the outputs of the real world analysis, we would need >5000 GB
1878.502 * 2723 / 1e3

#### Define implementation strategy
# Handle iterations in batches e.g., 1:100, 101:200, ..., N
# For each batch:
# * Run particle algorithms in parallel (Julia)
# * Collate states (Julia)
# * Make maps & compute diagnostics (R)
# * Unlink large particle files (R) 
# This minimises disk space requirements
# * Each batch requires only ~200 GB of storage space: 1878.502 * 100 / 1e3 
# * All batches ultimately require only:
# - (0.00189 + 3.550282 + 0.068324) * 210 = 760 MB (simulations)
# - (0.00189 + 3.550282 + 0.068324) * 2723 = 9858.611 MB (real-world analysis)


#### End of code. 
###########################
###########################