###########################
###########################
#### run-patter.R

#### Aims
# 1) This script provides a generic workflow for running particle algorithms 

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

#### (1) Run Julia workflows 

## (A) Run algorithms 
# Run 001-run-algorithms.jl via tmux/bash

## (B) Collate states 
# Run 001-run-algorithms.jl via tmux/bash

#### (2) Run R workflows 

## (A) Select iteration & load data 
# TO DO Use command_args here
it <- iteration[1, ]

## (B) Load iteration-specific datasets
it_states <- it_diagnostics <- NULL
if (file.exists(it$file_states)) {
  it_states <- arrow::read_feather(it$file_states) |> setDT()
}
if (file.exists(it$file_diagnostics)) {
  it_diagnostics <- arrow::read_feather(it$file_diagnostics) |> setDT()
}

## (B) Make map (~5.06 s)
if (!is.null(it_states)) {
  tic()
  occupancy <- map_pou(.map = map, .coord = it_states, .plot = FALSE)$ud
  terra::writeRaster(occupancy, it$file_occupancy, overwrite = TRUE)
  toc()
}

#### (4) Update diagnostics
# Compute areas spanned by 50 % and 95 % of the distribution 
if (!is.null(it_states) & !is.null(it_diagnostics)) {
  
  #### Compute area (m^2) spanned by 50 % distribution (~7.219 s)
  tic()
  area_core <- particle_hr(.map = map, 
                           .coord = it_states, 
                           .percentage = 0.50,
                           .summarise = FALSE)
  toc()
  
  #### Compute area  (m^2) spanned by 95 % distribution (~6.445 s)
  tic()
  area_home <- particle_hr(.map = map, 
                           .coord = it_states, 
                           .percentage = 0.95,
                           .summarise = FALSE)
  toc()

  #### Update diagnostics
  stopifnot(all(area_core$area <= area_home$area))
  it_diagnostics[, m2_core := area_core$area[match(timestep, area_core$timestep)]]
  it_diagnostics[, m2_home := area_home$area[match(timestep, area_home$timestep)]]
  
  #### Visual check
  if (FALSE) {
    # Check areas (m^2) for the first time step
    it_diagnostics[1, ]
    # Check areas (number of grid cells)
    A <- prod(terra::res(map))
    it_diagnostics$m2_core[1] / A
    it_diagnostics$m2_home[1] / A
    # Check areas (percentage of entire study area, excluding NAs)
    # > This is < 1 % of the size of the study area
    A <- terra::expanse(map)[2]
    it_diagnostics$m2_core[1] / A * 100
    it_diagnostics$m2_home[1] / A * 100
    # Visualise entire area spanned by distribution 
    occupancy_1 <- map_pou(.map = map, .coord = it_states[timestep == 1L, ])$ud
    terra::plot(occupancy_1 > 0)
  }

  #### Overwrite it$file_diagnostics
  write_feather_compressed(it_diagnostics, it$file_diagnostics)
  # arrow::read_feather(it$file_diagnostics)
  
}

#### (5) Clean up 
# We only store callstats, diagnostics and maps
# We delete iteration$file_states (these files are generally >= 678.995 MB)
# file.size(it$file_states) / 1e6 # 678.995 MB for "sim" iteration[1, ]
unlink(it$file_states)



#### End of code. 
###########################
###########################