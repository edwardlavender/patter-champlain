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
map        <- terra::rast(here_input("map.tif"))
pars       <- qs::qread(here_input("pars-patter.qs"))
detections <- qs::qread(here_input_real("detections.qs"))


###########################
###########################
#### Define unitsets

# We define unitsets, then the datasets, then iteration
# This matches the simulation workflow
# Not all units (individual/time blocks) pass quality control proceedures
# so not all units appear in the iteration data.table

#### Define units (individual/month combinations)
detections <- 
  detections |> 
  group_by(individual_id) |> 
  mutate(time_id = lubridate::floor_date(timestamp, "months")) |>
  select(individual_id, time_id, timestamp, receiver_id) |>
  as.data.table()

#### Define unitsets
unitsets <- 
  detections |> 
  group_by(individual_id, time_id) |> 
  slice(1L) |> 
  ungroup() |> 
  mutate(unit_id = row_number()) |> 
  select(unit_id, individual_id, time_id, timestamp, receiver_id) |>
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
#### Filter dataset

#### Number of individuals in full dataset
length(unique(detections$individual_id))

#### (optional) Focus on individuals with sufficient data
# This is no longer implemented
# We include all individuals
# But focus on individual/month blocks that meet selected criteria (below)

#### Focus on individual/month units with sufficient data
# Check the number of unit_ids in the raw data:
detections[, unit_id := .GRP, by = c("individual_id", "time_id")]
detections <- filter_detections(detections)

#### Write to file
# unlink(unitsets$file_detections)
detections[, file_detections := unitsets$file_detections[match(unit_id, unitsets$unit_id)]]
cl_lapply(split(detections, detections$unit_id), function(d) {
  qs::qsave(d, d$file_detections[1])
})
detections[, file_detections := NULL]


###########################
###########################
#### Prepare iteration patter: main analysis

#### Define iteration 
# TO DO Limit iteration rows
# * The number of rows in this data.table is too high
# * We should restrict this
#   - A) Restrict individual inclusion criteria 
#   - B) Reconsider simulation priorities 
#   - C) Restrict sensitivity analyses
#        on the basis of simulation results the sensitivity analysis 
iteration <- 
  unitsets |> 
  filter(unit_id %in% detections$unit_id) |> 
  select(unit_id, individual_id, time_id, 
         file_detections,
         folder_home = folder_home_patter) |>
  cross_join(pars) |> 
  mutate(index = row_number(),
         folder_coord = file.path(folder_home, "coord", parameter_id)) |> 
  as.data.table()

#### Build directories 
if (FALSE) {
  unlink(iteration$folder_coord)
}
nrow(iteration)
dirs.create(iteration$folder_coord)

#### Record iteration
qs::qsave(iteration, here_input_real("iteration-patter.qs"))


###########################
###########################
#### Prepare iteration heuristics

# TO DO (MF)


#### End of code. 
###########################
###########################