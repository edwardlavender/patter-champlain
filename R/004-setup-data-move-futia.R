###########################
###########################
#### setup-data-move-futia.R

#### Aims
# 1) Set up VPS data for analysis

#### Prerequisites
# 1) Movement data provided by Futia & colleagues


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
# load lake trout detections
di_lkt <- qs::qread(here_data_raw("model-move","drummond-island-2014",
                                  "DI_lkt_2013-2014_May2025.qs"))
tb_lkt <- qs::qread(here_data_raw("model-move","thunder-bay-2024",
                                  "TB_lkt_2021-2024_May2025.qs"))

# load sync tag detections
di_sync <- qs::qread(here_data_raw("model-move","drummond-island-2014",
                                   "DI_sync_2013-2014_May2025.qs"))
tb_sync <- qs::qread(here_data_raw("model-move","thunder-bay-2024",
                                   "TB_sync_2021-2024_May2025.qs"))

###########################
###########################
#### Process data

#### Run HPE filter for Drummond Island data
# Filter out any positions without an associated HPE or HPEm
di_sync_short <- 
  di_sync |> 
  filter(!is.na(HPEm) & !is.na(HPE))

# Remove top 5% HPE estimates
di_sync_filtered <- 
  di_sync_short |>
  filter(HPE <= ceiling(quantile(HPE, 0.95)))

# Create bins for HPE
di_sync_filtered$hpe_interval <- ceiling(di_sync_filtered$HPE)

# Summarize HPEm by HPE bin
di_hpe_filter <- 
  di_sync_filtered |>
  group_by(hpe_interval) |>
  reframe(num_pos = n(), 
          medianHPEm = median(HPEm),
          seventyFiveHPEm = quantile(HPEm, 0.75),
          ninetyHPEm = quantile(HPEm, 0.90), 
          ninetyFiveHPEm = quantile(HPEm, 0.95))

# Calculate the cumsum of total number of positions as HPE increases
di_hpe_filter$cum_sum <- cumsum(di_hpe_filter$num_pos)

# Convert cumsum to percentage of total positions
di_hpe_filter$cum_percent <- di_hpe_filter$cum_sum/nrow(di_sync_short)*100

di_hpe_filter

# Remove lkt detections where HPE associated with > 20 m error for 95% of data
di_hpe_cutoff <- 
  di_hpe_filter |> 
  filter(ninetyFiveHPEm < 20) |> 
  reframe(co = max(hpe_interval))

di_lkt_filtered <- 
  di_lkt |> 
  mutate(hpe_int = ceiling(HPE),
         year_detect = year(PosUTS),
         PosUTS = as.POSIXct(PosUTS),
         # convert length from cm to mm
         FishLen = FishLen/100) |> 
  filter(hpe_int < di_hpe_cutoff[[1]]) |> 
  select(animal_id = TagID,
         FullId = VUEID,
         origin = TreatGroup,
         sex = Sex,
         length_tl = FishLen,
         glatos_release_date_time = RelUTS,
         Time = PosUTS,
         Longitude, Latitude,hpe_int,year_detect)

summary(di_lkt_filtered)

# Write to file
qs::qsave(di_lkt_filtered, file = here_data("supp","model-move",
                                            "DI_lkt_2013-2014_filtered_May2025.qs"))


#### Run HPE filter for Thunder Bay data
# Filter out any positions without an associated HPE or HPEm
tb_sync_short <- 
  tb_sync |> 
  filter(!is.na(HPEm) & !is.na(HPE))

# Remove top 5% HPE estimates
tb_sync_filtered <- 
  tb_sync_short |>
  filter(HPE <= ceiling(quantile(HPE, 0.95)))

# Create bins for HPE
tb_sync_filtered$hpe_interval <- ceiling(tb_sync_filtered$HPE)

# Summarize HPEm by HPE bin
tb_hpe_filter <- 
  tb_sync_filtered |>
  group_by(hpe_interval) |>
  reframe(num_pos = n(), 
          medianHPEm = median(HPEm),
          seventyFiveHPEm = quantile(HPEm, 0.75),
          ninetyHPEm = quantile(HPEm, 0.90), 
          ninetyFiveHPEm = quantile(HPEm, 0.95))

# Calculate the cumsum of total number of positions as HPE increases
tb_hpe_filter$cum_sum <- cumsum(tb_hpe_filter$num_pos)

# Convert cumsum to percentage of total positions
tb_hpe_filter$cum_percent <- tb_hpe_filter$cum_sum/nrow(tb_sync_short)*100

tb_hpe_filter

# Remove lkt detections where HPE associated with > 20 m error for 95% of data
tb_hpe_cutoff <- 
  tb_hpe_filter |> 
  filter(ninetyFiveHPEm < 20) |> 
  reframe(co = max(hpe_interval))

tb_lkt_filtered <- 
  tb_lkt |> 
  mutate(hpe_int = ceiling(HPE),
         Time = as.POSIXct(Time)) |> 
  filter(hpe_int < tb_hpe_cutoff[[1]]) |> 
  select(animal_id,FullId,origin,sex,length_tl,glatos_release_date_time,Time,Longitude,Latitude,hpe_int,year_detect)

summary(tb_lkt_filtered)

# Write to file
qs::qsave(tb_lkt_filtered, file = here_data("supp","model-move",
                                            "TB_lkt_2021-2024_filtered_May2025.qs"))


#### End of code. 
###########################
###########################