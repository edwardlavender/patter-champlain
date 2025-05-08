###########################
###########################
#### setup-data-move-blanchfield.R

#### Aims
# 1) Set up Blanchfield et al. (2023) data for analysis

#### Prerequisites
# 1) Movement data provided by Blanchfield et al. (2023)


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Set global options
Sys.setenv("JULIA_SESSION" = FALSE)

#### Load essential packages
library(proj.verse)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = TRUE)
files_source_r(here_src())

#### Load data
blanchfield <- fread(here_data_raw("model-move", "blanchfield-et-al-2023", 
                                   "Alexie Accelerometer Data 20250307.csv"))


###########################
###########################
#### Setup Blanchfield dataset

# Examine raw data
head(blanchfield)

# Clean data.table
blanchfield <- 
  blanchfield |> 
  select(individual_id = Transmitter, 
         timestamp = DateTime_MST, 
         accel = Accel) |> 
  as.data.table()

# Quality checks
hist(blanchfield$accel, breaks = 100, xlim = range(blanchfield$accel))
nrow(blanchfield)            # 632,503
table(blanchfield$accel > 2) # 6484
table(blanchfield$accel > 3) # 3744
table(blanchfield$accel > 4) # 2567
table(blanchfield$accel > 5)

# Write to file
qs::qsave(blanchfield, here_data("supp", "model-move", "blanchfield.qs"))


###########################
###########################
#### Setup VPS data

# TO DO (MF)


#### End of code. 
###########################
###########################