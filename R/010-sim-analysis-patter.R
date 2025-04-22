###########################
###########################
#### sim-analysis-patter.R

#### Aims
# 1) Analyses simulated observations using patter

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
library(patter.workflows)
library(proj.verse)
library(tictoc)
files_source_r(here_src())
expect_no_geospatial()

#### Load data
iteration <- qs::qread(here_input_sim("iteration-patter.qs"))


###########################
###########################
#### Estimate coordinates

#### (optional) Reset directories
if (FALSE) {
  unlink(iteration$folder_coord, recursive = TRUE)
  dirs.create(iteration$folder_coord)
}

#### Select iterations 
stopifnot(!any(duplicated(iteration$index)))
iteration <- iteration[sensitivity == "best", ]
iteration[, file_coord := file.path(folder_coord, "coord.qs")]
iteration[, file_output := file_coord]

#### Set maps
set_map(here_input("map.tif"))
set_vmap(.vmap = here_input("vmap", iteration$mobility[1], "vmap.tif"))

#### Estimate coordinates
# TO DO
# * Add success methods for cl_lapply_workflow
# * Update constructor function e.g., with xinit 
# * Develop parallelisation (with julia_connect(.socket = TRUE))
coord_list <- 
  cl_lapply_workflow(.iteration   = iteration[1:2L, ],
                     .datasets    = list(),
                     .constructor = constructor_ac_sim, 
                     .algorithm   = patter.workflows::estimate_coord_particle)

#### Examine coordinates
if (FALSE) {
  # This code uses {terra}
  map <- terra::rast(here_input("map.tif"))
  list.files(iteration$folder_coord)
  lapply_qplot_coord(iteration, 
                     .datasets = list(coordinates = function(x) x$smooth$states),
                     .map = map)
}


#### End of code. 
###########################
###########################