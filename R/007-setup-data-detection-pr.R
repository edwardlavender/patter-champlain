###########################
###########################
#### setup-data-detection-pr.R

#### Aims
# 1) Collates detection probability datasets & prepares data for modelling 

#### Prerequisites
# 1) Process detection probability datasets


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
klinard <- qs::qread(here_data("supp", "model-obs", "klinard-raw.qs"))
futia   <- qs::qread(here_data("supp", "model-obs", "futia-raw.qs"))


###########################
###########################
#### Collate datasets

#### Merge datasets
# (optional) Move computation of detection counts from Klinard & Futia scripts here
# (to avoid code duplication)
klinard[, study := "K"]
futia[, study := "F"]
dcounts <- rbind(klinard, futia)
dcounts[, study := factor(paste0(study, dB))]
levels(dcounts$study)
str(dcounts)

#### Write to file
qs::qsave(dcounts, here_data("supp", "model-obs", "range-testing.qs"))


#### End of code. 
###########################
###########################