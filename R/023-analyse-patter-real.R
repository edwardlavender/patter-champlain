###########################
###########################
#### analysis-patter-real.R

#### Aims
# 1) This script provides further specific analysis of simulation outputs 

#### Prerequisites
# 1) Run previous scripts


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Set global options
Sys.setenv("JULIA_SESSION" = FALSE)
set.seed(123L)

#### Load essential packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(patter)
library(patter.workflows)
library(prettyGraphics)
library(proj.verse)
library(spatial.extensions)
library(tictoc)
files_source_r(here_src())

#### Load data
map           <- terra::rast(here_input("map.tif"))
champlain_utm <- qreadvect(here_input("champlain-utm.qs"))
regions_cs    <- qs::qread(here_input("regions-colour-scheme.qs"))
fish          <- qs::qread(here_input("fish.qs"))
moorings      <- qs::qread(here_input_real("main", "moorings.qs"))
iteration     <- qs::qread(here_input_real("main", "iteration.qs"))
convergence   <- qs::qread(here_output_real("main", "synthesis", "convergence.qs"))


###########################
###########################
#### Mapping 
# This code maps the overall pattern of space use by population & season

#### Select iterations
# 1. Define success or failure of each block
# * blocks have success if julia = TRUE and the algorithms converged
# * blocks have automatic success if julia = FALSE
# 2. Focus on chains where all blocks succeed
iteration <- 
  iteration |> 
  left_join(convergence |> 
              select("individual_id", "sensitivity", "block_id", success) |> 
              as.data.table(),
            by = c("individual_id", "sensitivity", "block_id")) |> 
  mutate(success = if_else(julia == FALSE, TRUE, success)) |> 
  group_by(individual_id, chain_id, sensitivity) |> 
  filter(all(success)) |> 
  ungroup() |> 
  select("individual_id", "sensitivity", "sensitivity_label", 
         "chain_id", "chain_start", "chain_end", "survival_probability",
         "file_occupancy", "file_residency") |>
  as.data.table()
# Verify that files exist exists for all chains
iteration[!file.exists(file_occupancy), ]
stopifnot(all(file.exists(iteration$file_occupancy)))
stopifnot(all(file.exists(iteration$file_residency)))
# Define sites/seasons
iteration <- 
  iteration |>
  # Define sites (recoded to N, S for brevity on plots) & season
  mutate(site = fish$site[match(individual_id, fish$individual_id)], 
         site = case_match(site, "Grand Isle" ~ "N", "Split Rock" ~ "S"), 
         season = season_factor(chain_start)) |>
  as.data.table()

#### Aggregate maps for each sensitivity/tagging location/season (~11 s)
# Define data.table of sensitivity/tagging location/season combinations
# Rows = sensitivity combinations
# Columns = site/season combinations
map_dt <- 
  iteration |> 
  select(sensitivity, sensitivity_label, site, season) |> 
  group_by(sensitivity, sensitivity_label, site, season) |> 
  slice(1L) |> 
  mutate(
    site_ = tolower(stringr::str_replace_all(site, " ", "-")), 
    season_ = tolower(stringr::str_replace_all(season, " ", "-")),
    file_ud = here_output_real(
      "main", "synthesis", 
      paste0(site_, "-",  season_, "-", iteration$sensitivity[1], ".tif"))
    ) |>
  select(-site_, -season_) |> 
  mutate(
    # Set row and column order on maps
    row = sensitivity_label, 
    column = paste0(season, " (", site, ")"), 
    
    column = factor(column, 
                    levels = c("Winter (N)", "Winter (S)", 
                               "Spring (N)", "Spring (S)", 
                               "Summer (N)", "Summer (S)", 
                               "Fall (N)",   "Fall (S)"))) |> 
  select(sensitivity, sensitivity_label, site, season, row, column, file_ud) |> 
  as.data.table()

#### Aggregate maps for each sensitivity/tagging location/season 
# (~7 s for all 'best' maps)
probs <- seq(0.05, 1, by = 0.05)
cl_lapply(split(map_dt, seq_len(nrow(map_dt))), function(dt) {
  
  # dt <- map_dt[1, ]
  its <- iteration[sensitivity == dt$sensitivity & 
                     site %in% dt$site & 
                     season %in% dt$season,  ]
  stopifnot(nrow(its) > 0L)
  
  # Sum SpatRasters for season over individuals/years, accounting for survival probability
  occupancy <- terra::rast(its$file_occupancy[1]) * its$survival_probability[1]
  for (i in 2:nrow(its)) {
    occupancy <- sum(occupancy, 
                     terra::rast(its$file_occupancy[i]) * its$survival_probability[i], 
                     na.rm = TRUE)
  }
  
  # Renormalise SpatRasters
  occupancy <- spatNormalise(occupancy) # occupancy / nrow(its)
  occupancy <- terra::classify(occupancy, cbind(0, NA))
  stopifnot(isTRUE(all.equal(1, terra::global(occupancy, "sum", na.rm = TRUE)[1, 1])))
  # terra::plot(occupancy)
  
  # (optional) Convert to home raster rasters
  # * For visualisation, it is more effective to use a few colours 
  #   rather than a continuous colour scheme
  # * (Most areas are relatively low probability)
  # * Hence, we use map_hr_prop() with .prop values & then sum up the outputs
  # * For each .prop, map_hr_prop() returns a value between 0 and 1
  # - By summing up, we generate a raster with values between 0 and N levels
  # - I.e., if p is c(0.5, 0.95): 
  #   - a value of 1 means in 95 % contour
  #   - value of 2 means in both (i.e., 0.5 % contour)
  occupancy_contours <- terra::setValues(occupancy, 0)
  for (p in probs) {
    occupancy_contours <- sum(occupancy_contours, 
                              map_hr_prop(occupancy, .prop = p), 
                              na.rm = TRUE)
  }
  occupancy_contours <- terra::classify(occupancy_contours, cbind(0, NA))
  # terra::plot(occupancy_contours)

  # Write to file
  terra::writeRaster(occupancy_contours, dt$file_ud, overwrite = TRUE)
  nothing()
  
})

#### Visualise map for example individual (includes legend)
# lapply_qplot_ud(map_dt, .n_plot = nrow(map_dt))
png(here_fig_real("main", "maps-example-with-legend.png"), 
height = 6, width = 2, units = "in", res = 800)
map_dt[1, ]
occupancy_eg <- terra::rast(map_dt$file_ud[1])
zlim <- c(0, length(probs))
terra::plot(occupancy_eg, range = zlim, 
            col = grDevices::terrain.colors(length(probs), rev = TRUE),
            plg = list(at = pretty(zlim),
                       labels = add_lagging_point_zero(pretty(zlim) / max(zlim))))
terra::lines(champlain_utm)
dev.off()

#### Plot 'best' maps
png(here_fig_real("main", "maps-best.png"), 
    height = 6, width = 10, units = "in", res = 800)
p <- ggmaps(.mapdt = 
              map_dt |> 
              filter(row == "Best") |> 
              mutate(row = factor(row)) |> 
              as.data.table(),
            .zlim = zlim,
            .map = map, 
            .coast = champlain_utm, .geom_coast = list(colour = "black"), 
            .moorings = moorings, .geom_moorings = list(size = 0.5, stroke = 0.5))
print(p)
dev.off()

#### Plot full maps 
# For manuscript, manually rearrange figure.
png(here_fig_real("main", "maps-sensitivity.png"), 
    height = 40, width = 10, units = "in", res = 800)
p <- ggmaps(.mapdt = map_dt,
            .map = map, 
            .coast = champlain_utm, .geom_coast = list(colour = "black"), 
            .moorings = moorings, .geom_moorings = list(size = 0.5, stroke = 0.5))
print(p)
dev.off()


###########################
###########################
#### Residency 

# This code plot the overall pattern of residency by population & season
# See also Fig. 6 in Futia et al. (2024)

#### Define residency
residency <- 
  cl_lapply(iteration$file_residency, qs::qread) |> 
  rbindlist() |> 
  filter(sensitivity == "best") |>
  mutate(site = fish$site[match(individual_id, fish$individual_id)], 
         site = case_match(site, "Grand Isle" ~ "N", "Split Rock" ~ "S"), 
         season = season_factor.ys(chain_id)) |>
  mutate(region = factor(region, levels = levels(regions_cs$region)), 
         col = regions_cs$col[match(region, regions_cs$region)]) |> 
  as.data.table()

#### Precompute mean residencies, accounting for survival probability
# mean/median values are only slightly adjusted by accounting for survivorship
residency_stats <- 
  residency |>
  group_by(site, season, region) |>
  summarise(
    # Compute unweighted/weighted mean
    mean_utd   = mean(estimate),
    mean_wtd   = weighted.mean(estimate, survival_probability),
    # Compute unweighted/weighted median
    # * There are different ways of computed a weighted median
    # * ggplot uses quantreg 
    # * Plotting median_utd and median_wtd on the plot below confirms ggplot2 actions weights appropriately
    median_utd = median(estimate), 
    # median_wtd = matrixStats::weightedMedian(estimate, survival_probability), 
    median_wtd = weighted.median(estimate, survival_probability)) |> 
  ungroup() |> 
  as.data.table()
# Accounting for survivorship makes less than 1 % of difference to mean residency estimates
# (but makes a bit more difference to the median)
hist((residency_stats$mean_utd - residency_stats$mean_wtd) * 100)
hist((residency_stats$median_utd - residency_stats$median_wtd) * 100)

#### Plot the distribution of residencies in each region for the best analysis
png(here_fig_real("main", "residency-best.png"), 
    height = 8, width = 6, units = "in", res = 800)
# Compute the number of observations per panel 
# * We will add this as a label to the plot below
n_per_panel <-
  residency |> 
  group_by(site, season) |> 
  # Count n for region[1] (other regions are the same)
  summarise(n = sum(!is.na(estimate[region == region[1]]))) |> 
  ungroup() |> 
  as.data.table()
# Make plot 
p <- 
  residency |> 
  ggplot(aes(region, estimate)) +
  geom_boxplot(aes(region, estimate, fill = I(col), weight = survival_probability), 
               varwidth = TRUE) +
  # Add jittered points, coloured by survival probability 
  geom_jitter(aes(region, estimate, alpha = survival_probability), 
              size = 0.25, 
              colour = "dimgrey",
              width = 0.1, height = 0) +
  # Add mean values, accounting for survival probability 
  # * The black points line up as expected, demonstrating ggplot2 actions weights properly
  # * Unweighted median values are generally, but not always, similar
  # geom_point(data = residency_stats, aes(region, median_wtd), 
  #            shape = 4, colour = "black", size = 2, stroke = 2) + 
  # geom_point(data = residency_stats, aes(region, median_utd),
  #            shape = 4, colour = "red", size = 1, stroke = 1) +
  # Add label of the number of observations per panel/region
  geom_text(data = n_per_panel,
            aes(x = -Inf, y = Inf, label = n),
            inherit.aes = FALSE,
            hjust = -0.25,   
            vjust = 1.5,
            size = 3) + 
  scale_y_continuous(# limits = c(0, 1),
                     breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = c(0, 1.2), clip = "off") +
  xlab("Region") +
  ylab(expression("Residency")) +
  # labs(alpha = "Survival probability") +
  facet_grid(season ~ site) +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
  )
print(p)
dev.off()

#### As above for the sensitivity analysis
# TO DO. 


#### End of code. 
###########################
###########################