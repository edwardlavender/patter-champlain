###########################
###########################
#### exports.R

#### Aims
# 1) Export selected datasets to data/exports for convenient use in other projects

#### Prerequisites
# 1) Run previous scripts


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Load essential packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)
files_source_r(here_src())


###########################
###########################
#### Simulation study 

#### Simulation datasets
paths <- 
  here_input_sim("main", "paths.qs") |> 
  qs::qread() |> 
  select(individual_id = path_id, 
         timestep, timestamp, 
         x, y) |> 
  as.data.frame()
detections <- 
  here_input_sim("main", "detections.qs") |> 
  qs::qread() |> 
  select(individual_id, timestamp, receiver_id, receiver_x, receiver_y) |> 
  as.data.frame()
qs::qsave(paths, here_data("export", "sim-paths.qs"))
qs::qsave(detections, here_data("export", "sim-detections.qs"))

#### Real-world datasets
# Read detections
iteration  <- qs::qread(here_input_real("main", "iteration.qs"))
moorings   <- qs::qread(here_input_real("main", "moorings.qs"))
detections <- 
  here_input_real("main", "detections.qs") |> 
  qs::qread() |> 
  filter(individual_id %in% iteration$individual_id) |> 
  filter(timestamp %within% interval(min(iteration$chain_start), 
                                     max(iteration$chain_end))) |> 
  as.data.frame()
range(detections$timestamp)
length(unique(detections$individual_id))
# Read residency
residency  <-
  here_output_real("main", "synthesis", "residency.qs") |> 
  qs::qread() |> 
  select(individual_id, tagging_site = site, 
         year_season = chain_id, season,
         survival_probability, 
         region, region_color = col, 
         residency_percent = perc) |> 
  as.data.frame()
# Write outputs
qs::qsave(moorings, here_data("export", "real-moorings.qs"))
qs::qsave(detections, here_data("export", "real-detections.qs"))
qs::qsave(residency, here_data("export", "real-residency.qs"))


#### End of code. 
###########################
###########################