###########################
###########################
#### prepare-real.R

#### Aims
# 1) Prepare inputs for real-world analyses 
# * Detection data
# * Iteration data.tables
# * Folder structure

#### Prerequisites
# 1) Process detection data


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
library(dplyr, warn.conflicts = FALSE)
library(proj.verse)
files_source_r(here_src())

#### Load data
map             <- terra::rast(here_input("map.tif"))
pars_model_move <- qs::qread(here_input("pars-model-move-best.qs"))
detections      <- qs::qread(here_input_real("detections.qs"))


###########################
###########################
#### Filter dataset

#### Number of individuals in full dataset
length(unique(detections$individual_id))

#### Focus on individuals with sufficient data
# This code is no longer implemented
# We include all individuals
# But focus on individual/month blocks that meet selected criteria (below)
if (FALSE) {
  # Compute metrics of data volume
  durations <- 
    detections |> 
    group_by(individual_id) |> 
    summarise(days = as.numeric(difftime(max(timestamp), min(timestamp), units = "days")), 
              ndays = length(unique(lubridate::floor_date(timestamp, "days")))) |> 
    as.data.table()
  # Summary statistics 
  # * 75 % of individuals were tracked for more than 167 days
  plot(ecdf(durations$days))
  quantile(durations$days, prob = 0.2)
  quantile(durations$days, prob = 0.25)
  # Focus on individuals tracked for longer than one month
  ids        <- durations$individual_id[durations$days > 31]
  detections <- detections[individual_id %in% ids, ]
  length(unique(detections$individual_id))
  # Focus on individuals with more than 31 days with detections
  ids        <- durations$individual_id[durations$ndays > 31]
  detections <- detections[individual_id %in% ids, ]
  length(unique(detections$individual_id))
}

#### Define individual/month units
detections <- 
  detections |> 
  group_by(individual_id) |> 
  mutate(time_id = lubridate::floor_date(timestamp, "months")) |>
  select(individual_id, time_id, timestamp, receiver_id) |>
  as.data.table()

#### Focus on individual/month units with sufficient data

# cf. patter-flapper criteria:
# - individuals must be detected in at least two different weeks 
# - on a total of seven days in a given month
# - (NB: this study included depth observations)

# patter-champlain criteria development:
# - Trout mobility is reasonable relative to the size of the study area
# - It would take a trout one day at top speed to cross the study area
# - Even at slower speeds, locations are uncertain after relatively short times
# - Without regular detections, locations are uncertain
ydim <- terra::ext(map)[4] - terra::ext(map)[3] # 174928.5 m
ydim / (24 * 60/2 * pars_model_move$mobility)   # 1 day to cross area

# Compute the proportion of days per month with detections
durations <- 
  detections |> 
  group_by(individual_id, time_id) |> 
  summarise(
    duration = as.numeric(difftime(max(timestamp), min(timestamp), units = "days")), 
    ndays = length(unique(lubridate::floor_date(timestamp, "days"))),
    pdays = ndays / duration
  ) |> 
  as.data.table()

# Examine the proportion of days with detections
quantile(durations$pdays)
plot(ecdf(durations$pdays), xlim = c(0, 1))
durations[pdays > 0.75, ]

# Select individual/month combinations with detections 75 % of days
detections[, unit_id := .GRP, by = c("individual_id", "time_id")]
durations[, unit_id := .GRP, by = c("individual_id", "time_id")]
length(unique(detections$unit_id))
detections <- detections[unit_id %in% durations$unit_id[durations$pdays >= 0.75], ]

# Redefine unit_ids
detections[, unit_id := .GRP, by = c("individual_id", "time_id")]
length(unique(detections$unit_id))


###########################
###########################
#### Define unitsets

#### Define unitsets
unitsets <- 
  detections |> 
  group_by(unit_id) |> 
  select(unit_id, individual_id, time_id, timestamp) |> 
  slice(1L) |> 
  mutate(
    file_detections = file.path("data", "input", "real", 
                                individual_id, time_id, "detection.qs"),
    folder_home = file.path("data", "output", "real", "main", "runs", 
                            individual_id, time_id), 
    folder_home_patter = file.path(folder_home, "patter")) |>
  as.data.table()

#### Build directories
if (FALSE) {
  unlink(file.path("data", "output", "real", "main", "runs"))
}
dirs.create(dirname(unitsets$file_detections))
dirs.create(unitsets$folder_home)
dirs.create(unitsets$folder_home_patter)

#### Write to file
qs::qsave(unitsets, here_input_real("unitsets.qs"))


###########################
###########################
#### Prepare iteration patter: main analysis

# TO DO Prepare iteration patter
# following prepare-sim.R 


###########################
###########################
#### Prepare iteration heuristics

# TO DO (MF)


#### End of code. 
###########################
###########################