###########################
###########################
#### analyse-model-comparison.R

#### Aims
# 1) Compare model performance based on accuracy and precision using simulation data
# 2) Evaluate influence on estimates of lake trout distributions using real-world data

#### Prerequisites
# 1) Simulation output from patter project
# 2) Detections from real-world dataset


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Load essential packages
library(proj.verse)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = TRUE)
files_source_r(here_src())

#### Load data
# Simulation datasets
paths <- qs::qread(here_data("export", "sim-paths.qs"))
detections <- qs::qread(here_data("export", "sim-detections.qs"))
residency_moe <- qs::qread(here_data("export", "sim-residency-moe.qs"))

# Real-world datasets
moorings <- qs::qread(here_data("export", "real-moorings.qs"))
detections <- qs::qread(here_data("export", "real-detections.qs"))
residency <- qs::qread(here_data("export", "real-residency.qs"))
