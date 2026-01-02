###########################
###########################
#### analysis-patter-sim.R

#### Aims
# 1) This script provides further specific analysis of simulation outputs 

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
map       <- terra::rast(here_input("map.tif"))
iteration <- qs::qread(here_input_sim("main", "iteration.qs"))


###########################
###########################
#### Analysis 

#### Occupancy error
# Compute mean error between simulated & reconstructed patterns of space use
# TO DO

#### Residency error
# Compute mean error between simulated & reconstructed residency estimates
# TO DO


#### End of code. 
###########################
###########################