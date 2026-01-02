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

#### (1) Run algorithms 
# Run 001-run-algorithms.jl via tmux/bash

#### (2) Collate states 
# Run 001-run-algorithms.jl via tmux/bash

#### (3) Make maps
# TO DO

#### (4) Update diagnostics
# Compute areas spanned by 50 % and 95 % of the distribution 
# TO DO 

#### (5) Clean up 
# We only store callstats, diagnostics and maps
# We delete iteration$file_states (these files are large)!


#### End of code. 
###########################
###########################