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
library(ggplot2)
files_source_r(here_src())

#### Load data
# Filtered Drummond Island data
di_lkt_filtered <- qs::qread(file = here_data("supp","model-move",
                                              "DI_lkt_2013-2014_filtered_May2025.qs"))
# Filtered Thunder Bay data
tb_lkt_filtered <- qs::qread(file = here_data("supp","model-move",
                                              "TB_lkt_2021-2024_filtered_May2025.qs"))


###########################
###########################
#### Analyse data

#### Fit models? (slow)
fit <- FALSE

#### Build Drummond Island step length dataset
# Calculate distance and time gaps between filtered positions
di_step <- 
  di_lkt_filtered |> 
  group_by(animal_id) |> 
  arrange(Time, .by_group = T) |> 
  mutate(site = "Drummond",
         step_gap = as.numeric(difftime(Time, lag(Time), units = "secs")),
         step_length = terra::distance(cbind(Longitude, Latitude), 
                                       cbind(lag(Longitude), lag(Latitude)), 
                                       lonlat = TRUE, pairwise = TRUE)) |> 
  ungroup()
# Add season
di_step <- 
  di_step |> 
  mutate(dayn = yday(Time),
         season = case_when(dayn < 121 | dayn >= 335 ~ "winter", # Dec 1 - Apr 30
                            dayn >= 121 & dayn < 182 ~ "spring", # May 1 - June 30
                            dayn >= 182 & dayn < 274 ~ "summer", # July 1 - Sept 30
                            dayn >= 274 & dayn < 335 ~ "fall"), # Oct 1 - Nov 30
         season = factor(season,
                         levels = c("winter", "spring", "summer", "fall")))
# Remove positions with transmission gap less than 120 secs & greater than 240 secs 
di_step_cut <- 
  di_step |> 
  filter(step_gap > 120 & step_gap < 360)
# Summarise processed position dataset
nrow(di_step)
length(unique(di_step$animal_id))
di_step |> 
  distinct(animal_id, .keep_all = TRUE) |> 
  reframe(utils.add::basic_stats(length_tl))
range(di_step$Time)
sort(unique(lubridate::round_date(di_step$Time, "months")))
# Calculate rate of movement (m/s)
di_step_cut <- 
  di_step_cut |> 
  mutate(move_rate_sec = step_length/step_gap,
         move_rate_2min = move_rate_sec * 120)
# Remove unlikely movements (greatest 5% for each fish)
di_step_95 <- 
  di_step_cut |> 
  group_by(animal_id) |> 
  mutate(quant_95 = quantile(move_rate_sec, 0.95)) |> 
  filter(move_rate_sec < quant_95) |> 
  ungroup()

#### Build Thunder Bay step length dataset
# Calculate distance and time gaps between filtered positions
tb_step <- 
  tb_lkt_filtered |> 
  group_by(animal_id) |> 
  arrange(Time, .by_group = T) |> 
  mutate(site = "Thunder",
         step_gap = as.numeric(difftime(Time, lag(Time), units = "secs")),
         step_length = terra::distance(cbind(Longitude, Latitude), cbind(lag(Longitude), lag(Latitude)), 
                                       lonlat = TRUE, pairwise = TRUE)) |> 
  ungroup()
# Add season
tb_step <- tb_step |> 
  mutate(dayn = yday(Time),
         season = "fall")
# Remove positions with transmission gap less than 120 secs & greater than 240 secs 
tb_step_cut <- 
  tb_step |> 
  filter(step_gap > 120 & step_gap < 360)
# Summarise processed position dataset
nrow(tb_step)
length(unique(tb_step$animal_id))
tb_step |> 
  distinct(animal_id, .keep_all = TRUE) |> 
  reframe(utils.add::basic_stats(length_tl))
range(tb_step$Time)
sort(unique(lubridate::round_date(tb_step$Time, "months")))
# Calculate rate of movement (m/s)
tb_step_cut <- 
  tb_step_cut |> 
  mutate(move_rate_sec = step_length/step_gap,
         move_rate_2min = move_rate_sec*120)
# Remove unlikely movements (greatest 5% for each fish)
tb_step_95 <- 
  tb_step_cut |> 
  group_by(animal_id) |> 
  mutate(quant_95 = quantile(move_rate_sec, 0.95)) |> 
  filter(move_rate_sec < quant_95) |> 
  ungroup()

#### Combine Drummond Island and Thunder Bay data
lkt_step <- 
  di_step_95 |> 
  bind_rows(tb_step_95) |> 
  mutate(move_rate_sec0 = if_else(move_rate_sec == 0, 0.00001, move_rate_sec))

#### Investigate Drummond Island dataset
# Compute summary statistics
di_step_95 |> 
  group_by(sex, season) |> 
  reframe(quart1_step = quantile(move_rate_2min, 0.25),
          med_step = median(move_rate_2min),
          quart3_step = quantile(move_rate_2min, 0.75),
          mean_step = mean(move_rate_2min))
# Fit GLMM
if (fit) {
  di_step0 <- 
    di_step_95 |> 
    mutate(move_rate_sec0 = if_else(move_rate_sec == 0, 0.00001, move_rate_sec))
  di_glmm <- glmmTMB::glmmTMB(move_rate_sec0 ~ sex+season+(1|animal_id),
                              data = di_step0,
                              family = glmmTMB::lognormal(),
                              na.action = na.fail)
  summary(di_glmm) 
  # > significant differences for sex and season 
  # > (greater step length during fall, p < 0.001, and for females, p = 0.003
}

#### Investigate Thunder Bay dataset
# Compute summary statistics
tb_step_95 |> 
  filter(sex %in% c("M", "F")) |> 
  group_by(sex) |> 
  reframe(quart1_step = quantile(move_rate_2min, 0.25),
          med_step = median(move_rate_2min),
          quart3_step = quantile(move_rate_2min, 0.75),
          mean_step = mean(move_rate_2min))
# Fit GLMM 
if (fit) {
  tb_step_sex <- 
    tb_step_95 |> 
    filter(sex %in% c("M", "F"))
  tb_glmm <- glmmTMB::glmmTMB(move_rate_sec ~ sex,
                              family = glmmTMB::lognormal(),
                              data = tb_step_sex,
                              na.action = na.fail)
  summary(tb_glmm) 
  # > Significant difference with greater movement for females (p < 0.001)
}

#### Examine combined dataset
# Examine variation in step length by site during fall
lkt_step |> 
  group_by(site) |> 
  reframe(quart1_step = quantile(move_rate_2min, 0.25),
          med_step = median(move_rate_2min),
          quart3_step = quantile(move_rate_2min, 0.75),
          mean_step = mean(move_rate_2min))
# Fit GLMM
if (fit) {
  fall_step <- 
    lkt_step |> 
    filter(season == "fall" & sex %in% c("M","F"))
  fall_glmm <- glmmTMB::glmmTMB(move_rate_sec0 ~ sex+site+(1|animal_id),
                                family = glmmTMB::lognormal(),
                                data = fall_step,
                                na.action = na.fail)
  summary(fall_glmm) 
  # > significant difference with greater movement for females (p = 0.009) 
  # > and Drummond Island fish (p < 0.001)
}

#### Visualise data using ggplot
# Plot distribution of rate of movement
lkt_step |> 
  filter(move_rate_2min > 0 &
           sex %in% c("M","F")) |> 
  ggplot() +
  geom_histogram(aes(x = move_rate_2min),
                 color = "black", fill = "gray80") +
  lemon::facet_rep_wrap(~site+sex, scales = "free") +
  labs(x = "Two-minute step length (meters)") +
  theme_classic()
# Plot rate of movement by HPE
lkt_step |> 
  ggplot(aes(x = hpe_int, y = move_rate_2min)) +
  geom_point(aes(color = move_rate_2min), alpha = 0.7) +
  scale_color_viridis_c() + 
  lemon::facet_rep_wrap(~site) +
  theme_classic()

#### Overall summaries (step length)
# Overall densities 
plot(density(lkt_step$step_length), ylim = c(0, 0.06))
lines(density(lkt_step$step_length[lkt_step$site == "Drummond"]), col = "blue")
lines(density(lkt_step$step_length[lkt_step$site == "Thunder"]), col = "red")
# Overall summary statistics
# * N obs by site
lkt_step |> 
  group_by(site) |> 
  summarise(n())
# * Movement speeds by site
lkt_step |> 
  # filter(move_rate_sec > 0) |> 
  group_by(site) |> 
  summarise(utils.add::basic_stats(move_rate_sec))
# * The range in median & maximum movement speeds
lkt_step |> 
  group_by(site, animal_id) |> 
  summarise(median = median(move_rate_sec), 
            max = max(move_rate_sec)) |> 
  ungroup() |> 
  group_by(site) |> 
  reframe(median = utils.add::basic_stats(median), 
          max = utils.add::basic_stats(max))
# Record dataset
lkt_step |> 
  select(site, step_length, move_rate_sec) |> 
  as.data.table() |> 
  qs::qsave(here_data("supp", "model-move", "futia-step.qs"))


###########################
###########################
#### Analyse turn angles

#### Calculate turn angle for Drummond Island fish
# Compute angles
di_angle <- 
  di_lkt_filtered |> 
  group_by(animal_id) |>  
  arrange(Time, .by_group = TRUE) |> 
  mutate(turn_angle = calc_angle_vps(pick(Longitude, Latitude)), 
         turn_time = Time - lag(Time, 2)) |>  
  ungroup() |> 
  mutate(turn_angle_deg = turn_angle * pi / 180)
# Filter and add season
di_angle_cut <- 
  di_angle |> 
  # Remove turns generated with >2 min time steps between detections
  filter(turn_time < 360) |> 
  # Assign detections to seasons
  mutate(dayn = yday(Time),
         season = case_when(dayn < 121 | dayn >= 335 ~ "winter", # Dec 1 - Apr 30
                            dayn >= 121 & dayn < 182 ~ "spring", # May 1 - June 30
                            dayn >= 182 & dayn < 274 ~ "summer", # July 1 - Sept 30
                            dayn >= 274 & dayn < 335 ~ "fall"),  # Oct 1 - Nov 30
         season = factor(season,
                         levels = c("winter", "spring", "summer", "fall")),
         # add site name
         site = "Drummond")
summary(di_angle_cut)

#### Calculate turn angle for Thunder Bay fish
# Compute angle 
tb_angle <- 
  tb_lkt_filtered |> 
  group_by(animal_id) |>  
  arrange(Time, .by_group = TRUE) |> 
  mutate(turn_angle = calc_angle_vps(pick(Longitude, Latitude)),
         turn_time = Time - lag(Time, 2)) |>  
  ungroup() |> 
  mutate(turn_angle_deg = turn_angle * pi / 180)
# Filter to shortest time steps 
# (season not evaluated as all detections occurred during fall)
tb_angle_cut <- 
  tb_angle |> 
  # remove turns generated with >3 min time steps between detections
  filter(turn_time < 540) |> 
  # assign site name
  mutate(site = "Thunder",
         season = "fall")
summary(tb_angle_cut)

#### Combine datasets
turn_angles <- 
  di_angle_cut |> 
  bind_rows(tb_angle_cut)

#### Investigate Drummond Island dataset 
# Compute summary statistics
di_angle_cut |> 
  group_by(sex, season) |> 
  reframe(quart1_step = quantile(turn_angle_deg, 0.25),
          med_step = median(turn_angle_deg),
          quart3_step = quantile(turn_angle_deg, 0.75),
          mean_step = mean(turn_angle_deg))
# Fit GLMM
if (fit) {
  di_angle_glmm <- glmmTMB::glmmTMB(turn_angle ~ sex+season+(1|animal_id),
                                    data = di_angle_cut,
                                    family = gaussian(),
                                    na.action = na.fail)
  summary(di_glmm)
}

#### Investigate Thunder Bay dataset
# Compute summary statistics
tb_angle_cut |> 
  filter(sex %in% c("M", "F")) |> 
  group_by(sex) |> 
  reframe(quart1_step = quantile(turn_angle_deg, 0.25),
          med_step = median(turn_angle_deg),
          quart3_step = quantile(turn_angle_deg, 0.75),
          mean_step = mean(turn_angle_deg))
# Fix GLMM
if (fit) {
  tb_angle_sex <- 
    tb_angle_cut |> 
    filter(sex %in% c("M", "F"))
  
  tb_angle_glmm <- glmmTMB::glmmTMB(turn_angle ~ sex,
                                    family = gaussian(),
                                    data = tb_angle_sex,
                                    na.action = na.fail)
  summary(tb_angle_glmm)
}

#### Investigate overall dataset 
# Examine variation in step length angle by site during fall
turn_angles |> 
  group_by(site) |> 
  reframe(quart1_step = quantile(turn_angle_deg, 0.25),
          med_step = median(turn_angle_deg),
          quart3_step = quantile(turn_angle_deg, 0.75),
          mean_step = mean(turn_angle_deg))
# Fit GLMM 
if (fit) {
  fall_turn_angle <- 
    turn_angles |> 
    filter(season == "fall" & 
             sex %in% c("M","F"))
  fall_angle_glmm <- glmmTMB::glmmTMB(turn_angle ~ sex+site+(1|animal_id),
                                      family = gaussian(),
                                      data = fall_turn_angle,
                                      na.action = na.fail)
  summary(fall_angle_glmm) 
}

#### Visualise dataset with ggplot2
# Density plot of turn angle standard deviations by site and sex
turn_angles |> 
  group_by(animal_id, site) |> 
  reframe(turn_sd = sd(turn_angle, na.rm = T)) |>  
  ggplot(aes(x = turn_sd, y = ..density..)) +
  geom_histogram(color = "black", fill = "transparent") +
  geom_density(color = "red", size = 0.75) +
  lemon::facet_rep_wrap(~site) +
  labs(x = "Turn angle standard deviation (radians)",
       y = "Density") +
  theme_classic() 
# Radial histogram of turn angle standard deviation by site (Drummond Island)
turn_angles |> 
  
  filter(site == "Drummond") |>
  ggplot() +
  geom_histogram(aes(x = turn_angle_deg,
                     y = after_stat(c(count[group == 1]/sum(count[group == 1]),
                                      count[group == 2]/sum(count[group == 2]))),
                     fill = season),
                 bins = 24, color = "black") +
  coord_polar(theta = "x") +
  lemon::facet_rep_wrap(~season) +
  scale_x_continuous(breaks = seq(0,360, 15)) +
  scale_fill_manual(values = c("#5ab4ac","#8c510a")) +
  labs(x = "Turn angle (degrees)", y = "Frequency", fill = "Season") +
  theme_light()
# As above for Thunder Island
turn_angles |> 
  filter(site == "Thunder" &
           sex %in% c("M","F")) |> 
  ggplot() +
  geom_histogram(aes(x = turn_angle_deg, 
                     y = after_stat(c(count[group == 1]/sum(count[group == 1]),
                                      count[group == 2]/sum(count[group == 2]))),
                     fill = sex),
                 bins = 24, color = "black") +
  coord_polar(theta = "x") +
  lemon::facet_rep_wrap(~ sex) +
  scale_x_continuous(breaks = seq(0,360, 15)) +
  scale_fill_manual(values = c("#5ab4ac","#8c510a")) +
  labs(x = "Turn angle (degrees)", y = "Frequency", fill = "Sex") +
  theme_light()

#### Overall summaries (turn angle)
# Overall densities
plot(density(turn_angles$turn_angle))
lines(density(turn_angles$turn_angle[turn_angles$site == "Drummond"]), col = "blue")
lines(density(turn_angles$turn_angle[turn_angles$site == "Thunder"]), col = "red")
# Overall summary statistics
turn_angles |> 
  group_by(site) |> 
  summarise(utils.add::basic_stats(turn_angle))
turn_angles |> 
  select(site, turn_angle) |> 
  as.data.table() |> 
  qs::qsave(here_data("supp", "model-move", "futia-angle.qs"))

# #### Overall summarise (turn angle SD)
# # Compute SD in turning angle
# turn_angles_sd <- 
#   turn_angles |> 
#   group_by(site, animal_id) |>
#   summarise(sd = sd(turn_angle)) |>
#   filter(!is.na(sd)) |> 
#   as.data.table()
# plot(density(turn_angles_sd$sd))
# lines(density(turn_angles_sd$sd[turn_angles_sd$site == "Drummond"]), col = "blue")
# lines(density(turn_angles_sd$sd[turn_angles_sd$site == "Thunder"]), col = "red")


#### End of code. 
###########################
###########################