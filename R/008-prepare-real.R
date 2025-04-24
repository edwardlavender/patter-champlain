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
detections <- qs::qread(here_input_real("detections.qs"))


###########################
###########################
#### Batch datasets

#### Method
# Batch the detection datasets:
# * Mitigate convergence issues
# * Improved memory handling 
# Split at the moment of a detection, roughly into one-month batches 
# Join batches @ detections (last time step, first time step)

#### Define individuals/months
detections <- 
  detections |> 
  group_by(individual_id) |> 
  mutate(unit_id = as.character(cut(timestamp, "months")), 
         unit_id = stringr::str_replace_all(unit_id, "-", "")) |> 
  select(individual_id, unit_id, timestamp, receiver_id) |>
  as.data.table()

#### Update detections 
detections <- 
  lapply(split(detections, detections$individual_id), function(d_id) {
    # d_id <- split(detections, detections$individual_id)[[1]]
    d_batch <- split(d_id, d_id$unit_id)
    n_batch <- length(d_batch)
    # (optional) Merge batches
    # * Merge batches if subsequent batches only contain a few additional observations
    if (n_batch > 1L) {
      for (i in n_batch:2L) {
        # Compute duration of current batch
        duration <- difftime(max(d_batch[[i]]$timestamp), 
                             min(d_batch[[i]]$timestamp), 
                             units = "weeks")
        # If duration is less than threshold (e.g., 1 week), merge with previous batch
        if (duration < 1) {
          # Combine datasets, using unit_id of preceeding batch
          d_batch[[i]][, unit_id := d_batch[[i - 1]]$unit_id[1]]
          d_batch[[i - 1]] <- rbind(d_batch[[i - 1]], d_batch[[i]])
          d_batch[[i]] <- NULL
        }
      }
      d_batch <- plyr::compact(d_batch)
    }
    n_batch <- length(d_batch)
    # Link batches
    # * The next batch should begin with the last observation from the previous batch
    # * This means we always start/stop filter @ a detection
    # * (i.e., in a 'happy' place) 
    if (n_batch > 1L) {
      for (i in 1:(n_batch - 1)) {
        last_row <- d_batch[[i]][.N, ]
        d_batch[[i + 1]] <- rbind(last_row, d_batch[[i + 1]])
      }
    }
    # Update unit_id
    # (use last unit_id as first one for batches 2:N is different)
    for (i in 1:n_batch) {
      d_batch[[i]][, unit_id := paste0(individual_id, "_", unit_id[.N])]
    }
    # Rejoin batches
    rbindlist(d_batch)
  }) |> rbindlist()

#### Check code works
det_1 <- detections[individual_id == 26786, ]
det_1 <- split(det_1, det_1$unit_id)
length(det_1)
(det_1a <- det_1[[1]][.N, ])
(det_1b <- det_1[[2]][1, ])
det_1[[2]][1:3, ]
stopifnot(all.equal(det_1a$individual_id, det_1b$individual_id))
stopifnot(all.equal(det_1a$timestamp, det_1b$timestamp))
stopifnot(all.equal(det_1a$receiver_id, det_1b$receiver_id))

#### Check the number of observations & time range per batch
nobs <- 
  detections |> 
  group_by(unit_id) |> 
  summarise(n = n(), 
            duration = as.numeric(difftime(max(timestamp), min(timestamp)),
                                  units = "days")) |> 
  arrange(n) |>
  as.data.table()
head(sort(nobs$n))
head(sort(nobs$duration)) 


###########################
###########################
#### Define unitsets

#### Record mapping between individual_id and unit_id
# TO DO 
# Clean this code & define unitsets 
detections_units <- 
  detections |> 
  select(individual_id, unit_id) |> 
  group_by(unit_id) |> 
  slice(1L) |> 
  ungroup() |> 
  group_by(individual_id) |>
  mutate(n_batch = n()) |> 
  ungroup() |>
  arrange(individual_id, unit_id) |>
  as.data.table()
# View(detections_units)

#### Build directories
# TO DO


###########################
###########################
#### Write to file

# TO DO


#### End of code. 
###########################
###########################