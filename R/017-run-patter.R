###########################
###########################
#### analysis-patter.R

#### Aims
# 1) This script provides a generic workflow for analysing simulated/real-world datasets with patter

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
#### Run patter

#### Run algorithms 
# Run 001-run-algorithms.jl via tmux/bash

#### Monitor progress
table(file.exists(iteration$file_callstats))
table(file.exists(iteration$file_diagnostics))

#### Estimate file size (MB)

# File size breakdown for simulations (1500 particles):
# * callstats   : 0.00189 
# * diagnostics : 3.550282
# * jld2 files  : 1071.414
# * states      : 1874.882 
# * map         : 0.068324 (see below)
# * Total       : 1878.502 (excluding jld2 files, which we delete)
# To store the outputs of the simulations, we need ~394 GB
1878.502 * nrow(iteration) / 1e3 

# File size breakdown for real-world analysis (2000 particles): 
# * callstats   : 
# * diagnostics : 
# * jld2 files  : 
# * states      : 
# * map         : 
# * total       : 
# To store the outputs of the real world analysis, we need >5000 GB
1878.502 * 2723 / 1e3

if (FALSE) {
  # Check file sizes for an example row
  i <- 2
  list.files(iteration$folder_output[i])
  dir_size(iteration$folder_output[i])
  file.size(iteration$file_callstats[i]) / 1e6 
  file.size(iteration$file_diagnostics[i]) / 1e6
  sum(file.size(list.files(iteration$folder_output[i], pattern = "\\.jld2$", full.names = TRUE))) / 1e6
  file.size(iteration$file_states[i]) / 1e6
}

#### Collate states 
# Run 001-run-algorithms.jl via tmux/bash
table(file.exists(iteration$file_states))


###########################
###########################
#### Examine example outputs

#### Define example individual
# Read data 
it <- iteration[2, ]
it_callstats   <- arrow::read_feather(it$file_callstats)
it_diagnostics <- arrow::read_feather(it$file_diagnostics)
it_states      <- arrow::read_feather(it$file_states)
# Check file sizes
# > We should be able to read all diagnostic files into memory 
nrow(it_diagnostics)
nrow(it_states)

#### Examine callstats
it_callstats
sum(it_callstats$time) # secs

#### Examine diagnostics
# Summarise ESS
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
#### Identify convergence 

#### Compute total run time
callstats <- lapply(iteration$file_callstats, \(f) arrow::read_feather) |> rbindlist()
sum(callstats$time)

#### Compute total output size (MB, GB)
# Compute folder sizes
iteration[, folder_output_mb := 
            sapply(seq_len(nrow(iteration)),
                   \(i) dir_size(iteration$folder_output[i], recursive = TRUE))]
# Check total size (GB)
sum(iteration$folder_output_mb) / 1e3

#### Identify convergence
# TO DO

#### Filter by convergence (for subsequent steps)
# TO DO


###########################
###########################
#### Summarise callstats

#### Compute computation time (for successful runs)
# TO DO, use ggplot


###########################
###########################
#### Summarise diagnostics

#### Summarise ESS
# (for filters & smoother)
# TO DO

#### Summarise area spanned by 95 % of the distribution
# TO DO, see previous version of this script


###########################
###########################
#### Map states

# TO DO Make maps via cl_lapply()

# TO DO Compute residency in each of the seven regions of interest

# TO DO For simulations
# * check example maps with trajectories
# * Compute map error 
# * Compute residency error 

# TO DO For real-world analysis
# * Create maps for example individuals
# * Collate maps


#### End of code. 
###########################
###########################