###########################
###########################
#### develop-model-move-futia.A

#### Aims
# 1) Develop a movement model for lake trout with analyses of VPS data 

#### Prerequisites
# 1) Process VPS data


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Set global options
Sys.setenv("JULIA_SESSION" = FALSE)
set.seed(1L)

#### Load essential packages
library(proj.verse)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
files_source_r(here_src())

#### Load data
# TO DO MF


###########################
###########################
#### Analyse data


# TO DO MF Add vps analysis of (A) step length and (B) turning angle



#### End of code. 
###########################
###########################