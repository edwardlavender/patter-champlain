###########################
###########################
#### setup-data-detection-pr-futia.R

#### Aims
# 1) Setup range testing data from Lake Champlain

#### Prerequisites
# 1) Detection probability data provided by Futia & colleagues


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Set global options
Sys.setenv("JULIA_SESSION" = FALSE)

#### Load essential packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = TRUE)
library(proj.verse)
files_source_r(here_src())

#### Load data
# TO DO MF


###########################
###########################
#### Process data

# TO DO MF Setup detection pr data, 
# following setup-data-detection-pr-klinard.R:
# We need columns select(transmitter_id, timestamp, receiver_id, dB, dist) 
# Output file: data/supp/model-obs/futia-raw.qs


#### End of code. 
###########################
###########################