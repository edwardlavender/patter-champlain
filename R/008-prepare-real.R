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
# Split into approximate individual/month batches
# Split at the moment of a detection
# Join batches @ detections (last time step, first time step)
# (optional) Drop batches with limited data

#### (optional) Focus on select individuals with good data
# Number of individuals 
length(unique(detections$individual_id))
# Frequency distribution for duration of detection time series
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

#### Define individuals/months
detections <- 
  detections |> 
  group_by(individual_id) |> 
  mutate(time_id = lubridate::floor_date(timestamp, "months")) |>
  select(individual_id, time_id, timestamp, receiver_id) |>
  as.data.table()

#### Split detections into individual/month blocks, merging small batches
detections <- 
  lapply(split(detections, detections$individual_id), function(d_id) {
    # d_id <- split(detections, detections$individual_id)[[1]]
    d_batch <- split(d_id, d_id$time_id)
    n_batch <- length(d_batch)
    # (optional) Merge batches
    # * Merge batches if subsequent batches only contain a few additional observations
    if (n_batch > 1L) {
      for (i in n_batch:2L) {
        # Compute duration of current batch
        duration <- difftime(max(d_batch[[i]]$timestamp), 
                             min(d_batch[[i]]$timestamp), 
                             units = "weeks")
        # If duration is less than threshold (e.g., 2 weeks), merge with previous batch
        if (duration < 2) {
          # Combine datasets, using time_id of preceding batch
          d_batch[[i]][, time_id := d_batch[[i - 1]]$time_id[1]]
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
    # Update time_id
    # (use first time_id as first one for batches 2:N is different)
    for (i in 1:n_batch) {
      d_batch[[i]][, time_id := paste0(individual_id, "_", time_id[1])]
    }
    # Rejoin batches
    rbindlist(d_batch)
  }) |> rbindlist()

#### Check code works
det_1 <- detections[individual_id == 26786, ]
det_1 <- split(det_1, det_1$time_id)
length(det_1)
(det_1a <- det_1[[1]][.N, ])
(det_1b <- det_1[[2]][1, ])
det_1[[2]][1:3, ]
stopifnot(all.equal(det_1a$individual_id, det_1b$individual_id))
stopifnot(all.equal(det_1a$timestamp, det_1b$timestamp))
stopifnot(all.equal(det_1a$receiver_id, det_1b$receiver_id))

#### Define unit_id
detections[, unit_id := .GRP, by = c("individual_id", "time_id")]
detections <- detections[, .(unit_id, individual_id, time_id, timestamp, receiver_id)]
# Checks
stopifnot(
  length(unique(detections$unit_id)) ==
  length(unique(paste(detections$individual_id, detections$time_id)))
)

#### Update unit_id
# Compute statistics for each unit_id
# * Number of observations
# * Number of weeks with observations
# * Duration between first and last observation 
unitstats <- 
  detections |> 
  group_by(unit_id) |> 
  summarise(
    n = n(), 
    n_day = length(unique(lubridate::floor_date(timestamp, "days"))),
    n_week = length(unique(lubridate::floor_date(timestamp, "weeks"))),
    duration = as.numeric(difftime(max(timestamp), min(timestamp)),
                          units = "days")) |> 
  arrange(n) |>
  as.data.table()
# Examine summary statistics
head(sort(unitstats$n))
head(sort(unitstats$duration)) 
hist(unitstats$n, breaks = 50)
hist(unitstats$duration, breaks = 50)
# Drop any unit_ids with insufficient data 
head(unitstats)


###########################
###########################
#### Define unitsets

#### Define unitsets
# TO DO 

#### Build directories
# TO DO


###########################
###########################
#### Write to file

# TO DO


#### End of code. 
###########################
###########################