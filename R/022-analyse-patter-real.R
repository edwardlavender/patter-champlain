###########################
###########################
#### analysis-patter-real.R

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
iteration <- qs::qread(here_input_real("main", "iteration.qs"))


###########################
###########################
#### Analysis

# TO DO analyse results


#### End of code. 
###########################
###########################