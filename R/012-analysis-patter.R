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
patter::julia_connect()

#### Load essential packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(JuliaCall)
library(patter)
library(patter.workflows)
library(proj.verse)
library(tictoc)
files_source_r(here_src())
expect_no_geospatial()


###########################
###########################
#### Select analysis type

analysis <- "sim"
analysis <- "real"

if (analysis == "sim") {

  here_input_analysis     <- here_input_sim
  here_output_analysis    <- here_output_sim_main
  constructor_ac_analysis <- constructor_ac_sim
  
} else if (analysis == "real") {
  
  here_input_analysis     <- here_input_real
  here_output_analysis    <- here_output_real_main
  constructor_ac_analysis <- constructor_ac_real
  
}

iteration <- qs::qread(here_input_analysis("iteration-patter.qs"))
nrow(iteration)


###########################
###########################
#### Estimate coordinates

#### (optional) Set development mode
# * Use one core
# * Set maps 
dev <- TRUE

#### (optional) Reset directories
if (TRUE) {
  # unlink(dirname(iteration$folder_coord), recursive = TRUE)
  unlink(iteration$folder_coord, recursive = TRUE)
  dirs.create(iteration$folder_coord)
}

#### Select iterations 
stopifnot(!any(duplicated(iteration$index)))
table(iteration$mobility)
iteration <- iteration[sensitivity == "best", ]
iteration[, file_diag := file.path(folder_coord, "diagnostics.qs")]
iteration[, file_output := file_diag]

#### Set maps
if (dev) {
  # This is implemented below on each node. 
  set_map(here_input("map.tif"))
  set_vmap(.vmap = here_input("vmap", iteration$mobility[1], "vmap.tif"))
}

#### Set logs
if (dev) {
  log.txt <- TRUE
} else {
  dir.create(here_output_analysis("logs"))
  log.txt <- here_output_analysis("logs", paste0("log-", iteration$mobility[1], ".txt"))
  # unlink(log.txt)
}

#### Setup cluster
if (!dev) {
  # Define number of workers
  ncl <- 2L
  # Define required memory for export
  lobstr::mem_used() * ncl
  # Initialise cluster
  cl  <- parallel::makeCluster(ncl)
  cl_init(iteration = iteration, cl = cl, varlist = ls())
} else {
  cl <- NULL
}

#### Estimate coordinates: time trials
# iteration[1, ], 1 thread, 2.5e4 filter particles, 1e3 smoothing particles, 1e2 smoothing sims
# * 26.8 min XX min: 6.95 min (filter) + 6.82 min (filter) + 11.98 min (smoother)
# * Warning: All smoothing weights (from xbwd[k, t] to xfwd[j, t - 1]) are zero at 8 time step(s) (0.04 %).

#### Estimate coordinates
# TO DO In patter.workflows, update .verbose for parallelisation
# debug(constructor_ac_core)
iteration <- iteration[1:1L, ]
coord_list <- 
  cl_lapply_workflow(.iteration   = iteration,
                     .datasets    = list(),
                     .constructor = constructor_ac_analysis, 
                     .algorithm   = estimate_coord_particle, 
                     .success     = particle_success, 
                     .cleanup     = particle_cleanup,
                     .cl          = cl,
                     .verbose     = log.txt)

#### Collate coordinates across batches
list.files(iteration$folder_coord)
stopifnot(all(file.exists(iteration$file_output)))
iteration[, file_coord := file.path(folder_coord, "coord.qs")]
convergence <- cl_lapply_particle_collate(.iteration = iteration)
# Check the number of algorithm runs for which convergence was achieved 
table(unlist(convergence))

#### Examine coordinates
if (FALSE) {
  
  # Examine example file
  list.files(iteration$folder_coord)
  eg <- qs::qread(iteration$file_coord[1])
  summary(eg)
  # View(eg)
  
  # Visualise coordinates
  map <- terra::rast(here_input("map.tif"))
  lapply_qplot_coord(iteration, 
                     .datasets = list(coordinates = function(smo) smo$states),
                     .map = map)
  
}


#### End of code. 
###########################
###########################