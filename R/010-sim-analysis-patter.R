###########################
###########################
#### sim-analysis-patter.R

#### Aims
# 1) Analyses simulated observations using patter
# 2) This code requires adaptation for Linux. 

#### Prerequisites
# 1) Run simulations


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
library(patter)
library(proj.verse)
library(tictoc)
files_source_r(here_src())

#### Load data
iteration          <- qs::qread(here_input_sim("iteration-patter.qs"))


###########################
###########################
#### Estimate coordinates

#### Select iterations 
iteration <- iteration[sensitivity == "best", ]
iteration[, file_coord := file.path(folder_coord, "coord.qs")]
iteration[, file_output := file_coord]

#### Set vmap
stopifnot(length(unique(iteration$mobility)) == 1L)
set_vmap(.map = map, .mobility = iteration$mobility[1])

#### Estimate coordinates
# TO DO
# * Develop parallelisation (with julia_connect(.socket = TRUE))
# * Add parallelisation check in constructor
# * Develop constructor function for simulations
# * Use simulated timeline 
# * Initialise models correctly from sim with helper functions
coord_list <- 
  cl_lapply_workflow(.iteration   = iteration,
                     .datasets    = datasets,
                     .constructor = constructor_ac_sim, 
                     .algorithm   = estimate_coord_particle)

#### Examine coordinates
list.files(iteration$folder_coord)
lapply_qplot_coord(iteration, 
                   .datasets = list(coordinates = function(x) x$smooth$states),
                   .map = map)


#### End of code. 
###########################
###########################