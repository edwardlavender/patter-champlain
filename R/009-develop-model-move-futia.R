###########################
###########################
#### develop-model-move-futia.A

#### Aims
# 1) Develop a movement model for lake trout with analyses of VPS data 

#### Prerequisites
# 1) Process VPS data


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Set global options
Sys.setenv("JULIA_SESSION" = FALSE)
set.seed(1L)

#### Load essential packages
library(proj.verse)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
files_source_r(here_src())

#### Load data
# Filtered Drummond Island data
di_lkt_filtered <- qs::qread(file = here_data("supp","model-move",
                                              "DI_lkt_2013-2014_filtered_May2025.qs"))
# Filtered Thunder Bay data
qs::qsave(tb_lkt_filtered, file = here_data("supp","model-move",
                                            "TB_lkt_2021-2024_filtered_May2025.qs"))

###########################
###########################
#### Analyse data
### 1. Lake trout step length
# calculate step length for Drummond Island fish
# calculate distance and time gaps between filtered positions
di_step <- 
  di_lkt_filtered |> 
  group_by(animal_id) |> 
  arrange(Time) |> 
  mutate(site = "Drummond",
         step_gap = as.numeric(Time - lag(Time)),
         step_length = terra::distance(cbind(Longitude, Latitude), cbind(lag(Longitude), lag(Latitude)), 
                                       lonlat = TRUE, pairwise = TRUE)) |> 
  ungroup()

# add season
di_step <- di_step %>% 
  mutate(dayn = yday(Time),
         season = case_when(dayn < 121 | dayn >= 335 ~ "winter", # Dec 1 - Apr 30
                            dayn >= 121 & dayn < 182 ~ "spring", # May 1 - June 30
                            dayn >= 182 & dayn < 274 ~ "summer", # July 1 - Sept 30
                            dayn >= 274 & dayn < 335 ~ "fall"), # Oct 1 - Nov 30
         pos_ssn = factor(season,
                          levels = c("winter", "spring", "summer", "fall")))

# remove positions with transmission gap less than 120 secs & greater than 180 secs 
di_step_cut <- 
  di_step |> 
  filter(step_gap > 120 & step_gap < 180)

# calculate rate of movement (m/s)
di_step_cut <- 
  di_step_cut |> 
  mutate(move_rate_sec = step_length/step_gap,
         move_rate_2min = move_rate_sec*120)

# remove unlikely movements (greatest 5% for each fish)
di_step_95 <- di_step_cut %>% 
  group_by(animal_id) %>% 
  mutate(quant_95 = quantile(move_rate_sec, 0.95)) %>% 
  filter(move_rate_sec < quant_95) %>% 
  ungroup()


# calculate step length for Thunder Bay fish
# calculate distance and time gaps between filtered positions
tb_step <- 
  tb_lkt_filtered |> 
  group_by(animal_id) |> 
  arrange(Time) |> 
  mutate(site = "Thunder",
         step_gap = as.numeric(Time - lag(Time)),
         step_length = terra::distance(cbind(Longitude, Latitude), cbind(lag(Longitude), lag(Latitude)), 
                                       lonlat = TRUE, pairwise = TRUE)) |> 
  ungroup()

# add season
tb_step <- tb_step %>% 
  mutate(dayn = yday(Time),
         season = "fall")

# remove positions with transmission gap less than 120 secs & greater than 240 secs 
tb_step_cut <- 
  tb_step |> 
  filter(step_gap > 120 & step_gap < 240)

# calculate rate of movement (m/s)
tb_step_cut <- 
  tb_step_cut |> 
  mutate(move_rate_sec = step_length/step_gap,
         move_rate_2min = move_rate_sec*120)

# remove unlikely movements (greatest 5% for each fish)
tb_step_95 <- tb_step_cut %>% 
  group_by(animal_id) %>% 
  mutate(quant_95 = quantile(move_rate_sec, 0.95)) %>% 
  filter(move_rate_sec < quant_95) %>% 
  ungroup()

# combine Drummond Island and Thunder Bay data
lkt_step <- 
  di_step_cut |> 
  bind_rows(tb_step_cut)

# plot distribution of rate of movement
lkt_step |> 
  filter(move_rate_2min > 0) |> 
  ggplot2::ggplot() +
  ggplot2::geom_histogram(ggplot2::aes(x = move_rate_2min),
                          color = "black", fill = "gray80") +
  lemon::facet_rep_wrap(~site, scales = "free") +
  ggplot2::labs(x = "Two-minute step length (meters)") +
  ggplot2::theme_classic()


### 2. Lake trout turn angle
# calculate turn angle for Drummond Island fish
di_angle <- 
  di_lkt_filtered |> 
  group_by(animal_id) |>  
  arrange(Time) |> 
  mutate(turn_angle = atan2(Longitude - lag(Longitude,2), Latitude - lag(Latitude,2)) -
           atan2(lag(Longitude,1)- lag(Longitude,2), lag(Latitude,1)- lag(Latitude,2)),
         turn_time = Time - lag(Time, 2)) |>  
  ungroup() |> 
  mutate(turn_angle_deg = if_else(turn_angle < 0, turn_angle*57.2958 + 360, turn_angle*57.2958))


# filter and add season
di_angle_cut <- 
  di_angle |> 
  # remove turns generated with >2 min time steps between detections
  filter(turn_time < 360) |> 
  # assign detections to seasons
  mutate(dayn = yday(Time),
         season = case_when(dayn < 121 | dayn >= 335 ~ "winter", # Dec 1 - Apr 30
                            dayn >= 121 & dayn < 182 ~ "spring", # May 1 - June 30
                            dayn >= 182 & dayn < 274 ~ "summer", # July 1 - Sept 30
                            dayn >= 274 & dayn < 335 ~ "fall"), # Oct 1 - Nov 30
         season = factor(season,
                         levels = c("winter", "spring", "summer", "fall")),
         # add site name
         site = "Drummond")

summary(di_angle_cut)


# calculate turn angle for Thunder Bay fish
tb_angle <- 
  tb_lkt_filtered |> 
  group_by(animal_id) |>  
  arrange(Time) |> 
  mutate(turn_angle = atan2(Longitude - lag(Longitude,2), Latitude - lag(Latitude,2)) -
           atan2(lag(Longitude,1)- lag(Longitude,2), lag(Latitude,1)- lag(Latitude,2)),
         turn_time = Time - lag(Time, 2)) |>  
  ungroup() |> 
  mutate(turn_angle_deg = if_else(turn_angle < 0, turn_angle*57.2958 + 360, turn_angle*57.2958))


# filter to shortest time steps (season not evaluated as all detections occurred during fall)
tb_angle_cut <- 
  tb_angle |> 
  # remove turns generated with >3 min time steps between detections
  filter(turn_time < 540) |> 
  # assign site name
  mutate(site = "Thunder",
         season = "fall")

summary(tb_angle_cut)

# combine data
turn_angles <- 
  di_angle_cut |> 
  bind_rows(tb_angle_cut)

# density plot of turn angle standard deviations by site
turn_angles |> 
  group_by(animal_id, site) |> 
  reframe(turn_sd = sd(turn_angle, na.rm = T)) |>  
  ggplot2::ggplot(ggplot2::aes(x = turn_sd, y = ..density..)) +
  ggplot2::geom_histogram(color = "black", fill = "transparent") +
  ggplot2::geom_density(color = "red", size = 0.75) +
  lemon::facet_rep_wrap(~site) +
  ggplot2::labs(x = "Turn angle standard deviation (radians)",
                y = "Density") +
  ggplot2::theme_classic() 

# Radial histogram of turn angle standard deviation by site
turn_angles |> 
  # focus on Drummond Island Data
  filter(site == "Drummond") |>
  ggplot2::ggplot() +
  ggplot2::geom_histogram(ggplot2::aes(x = turn_angle_deg,
                                       y = ggplot2::after_stat(c(count[group == 1]/sum(count[group == 1]),
                                                                 count[group == 2]/sum(count[group == 2]))),
                                       fill = season),
                          bins = 24, color = "black") +
  ggplot2::coord_polar(theta = "x") +
  lemon::facet_rep_wrap(~season) +
  ggplot2::scale_x_continuous(breaks = seq(0,360, 15)) +
  ggplot2::scale_fill_manual(values = c("#5ab4ac","#8c510a")) +
  ggplot2::labs(x = "Turn angle (degrees)", y = "Frequency", fill = "Season") +
  ggplot2::theme_light()

turn_angles |> 
  # focus on Thunder Bay data with sex identified
  filter(site == "Thunder" &
           sex %in% c("M","F")) |> 
  ggplot2::ggplot() +
  ggplot2::geom_histogram(ggplot2::aes(x = turn_angle_deg, 
                                       y = ggplot2::after_stat(c(count[group == 1]/sum(count[group == 1]),
                                                                 count[group == 2]/sum(count[group == 2]))),
                                       fill = sex),
                          bins = 24, color = "black") +
  ggplot2::coord_polar(theta = "x") +
  lemon::facet_rep_wrap(~ sex) +
  ggplot2::scale_x_continuous(breaks = seq(0,360, 15)) +
  ggplot2::scale_fill_manual(values = c("#5ab4ac","#8c510a")) +
  ggplot2::labs(x = "Turn angle (degrees)", y = "Frequency", fill = "Sex") +
  ggplot2::theme_light()

#### End of code. 
###########################
###########################