###########################
###########################
#### analysis-patter.R

#### Aims
# 1) This script provides generic analysis of patter outputs

#### Prerequisites
# 1) Run previous scripts


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
#### Visualise maps

#### Visualise example maps
# TO DO


#### End of code. 
###########################
###########################