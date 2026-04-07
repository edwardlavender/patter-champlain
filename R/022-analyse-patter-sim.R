###########################
###########################
#### analysis-patter-sim.R

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
library(ggplot2)
library(patter)
library(patter.workflows)
library(proj.verse)
library(spatial.extensions)
library(tictoc)
files_source_r(here_src())

#### Load data
map           <- terra::rast(here_input("map.tif"))
regions_cs    <- qs::qread(here_input("regions-colour-scheme.qs"))
iteration     <- qs::qread(here_input_sim("main", "iteration.qs"))
paths         <- qs::qread(here_input_sim("main", "paths.qs"))
champlain_utm <- qreadvect(here_input("champlain-utm.qs"))
moorings      <- qs::qread(here_input_sim("main", "moorings.qs"))


###########################
###########################
#### Select iterations

# For simulations, for comparison with previous work (Futia et al. 2024) and 
# the real-world analysis, we focus on iterations that generated detections
nrow(iteration)
iteration <- iteration[n_detections > 0L, ]
nrow(iteration)


###########################
###########################
#### Visualise simulated datasets

# This is useful to understand the causes of convergence failures (~2 mins)
# (NB: setting .cl produces an empty image)
if (FALSE) {
  tic()
  png(here_fig_sim("main", "paths.png"), 
      height = 8, width = 12, units = "in", res = 800)
  pp <- par(mfrow = c(5, 20), mar = c(0, 0, 0, 0))
  cl_lapply(split(paths, paths$path_id), function(path) {
    # path <- paths[path_id == 1L]
    terra::plot(map, 
                main = path$path_id[1],
                legend = FALSE, 
                col = "lightgrey", 
                pax = list(labels = FALSE))
    patter:::add_sp_path(path$x, path$y, lwd = 0.25, length = 0.01) |> 
      suppressWarnings()
    nothing()
  })
  dev.off()
  toc()
}


###########################
###########################
#### Compute skill metrics 

#### Compute ncell
nc <- terra::freq(map)[["count"]]

#### Compute occupancy skill
# Compute mean absolute error between simulated & reconstructed patterns of space use
# > This is primarily used to understand sensitivity
# > It is useful for future comparisons against other methods
# > Other metrics e.g., EMD are more interpretable in terms of how 'accurate' maps are
# > But for this work we simply compare simulated tracks & associated maps
overwrite <- FALSE
file_occupancy_skill <- here_output_sim_main("synthesis", "occupancy-skill.qs")
if (!file.exists(file_occupancy_skill) | overwrite) {
  
  occupancy_skill <- iteration[file.exists(file_occupancy), ]
  stopifnot(nrow(occupancy_skill) > 0L)
  
  skills <- 
    split(occupancy_skill, seq_len(nrow(occupancy_skill))) |> 
    cl_lapply(function(d) {
      
      #### Read rasters
      # d <- occupancy_skill[1, ]
      obs <- terra::rast(d$file_occupancy_sim)
      mod <- terra::rast(d$file_occupancy)
      
      #### Visual check on similarity 
      if (FALSE) {
        pp <- par(mfrow = c(1, 2))
        terra::plot(obs)
        terra::plot(mod)
        par(pp)
      }
      
      #### Compute MAE
      # This is the average absolute difference in prob density per pixel
      MAE <- skill_me(.obs = obs, .mod = mod) 
      
      #### Compute TVD
      # This is the fraction of total probability mass in the wrong place (0 --> 1)
      # But this penalises even small differences
      # TVD <- 0.5 * nc * MAE
      
      #### Compute Jensen–Shannon divergence
      # This ranges between 0 and 1 (maximally different)
      # JSD <- spatJS(obs, mod)
      
      #### Compute EMD (slow)
      # For obs, set cells outside 95 % contour to optimise the calculation 
      if (FALSE) {
        obs_home <- map_hr_home(obs)
        obs_home <- terra::classify(obs, cbind(0, NA))
        obs      <- terra::mask(obs, obs_home)
        obs      <- terra::classify(obs, cbind(NA, 0))
        obs      <- terra::aggregate(obs, fact = 5, "sum", na.rm = TRUE)
        obs      <- obs / terra::global(obs, "sum", na.rm = TRUE)[1, 1]
        # For sim, as above 
        mod_home <- map_hr_home(mod)
        mod_home <- terra::classify(mod, cbind(0, NA))
        mod      <- terra::mask(mod, mod_home)
        mod      <- terra::classify(mod, cbind(NA, 0))
        mod      <- terra::aggregate(mod, fact = 5, "sum", na.rm = TRUE)
        mod      <- mod / terra::global(mod, "sum", na.rm = TRUE)[1, 1]
        terra::res(mod); terra::ncell(mod)
        # Compute EMD (~20 s)
        # * Set threshold to 5 km 
        # * This speeds up calculation (20 s -> 2.1 s)
        # * We limit small discrepancies at large distances affecting the result
        tic()
        EMD <- move::emd(raster::raster(obs), raster::raster(mod), threshold = 5000)
        toc()
      }
 
      #### Return selected metric
      MAE
      
    }) |> 
    unlist()
  
  occupancy_skill[, skill := skills]
  qs::qsave(occupancy_skill, file_occupancy_skill)
  occupancy_skill
  
} else {
  occupancy_skill <- qs::qread(file_occupancy_skill)
}

#### Compute residency skill
# Compute error between simulated & reconstructed residency estimates _by region_
overwrite <- FALSE
file_residency_skill <- here_output_sim_main("synthesis", "residency-skill.qs")
if (!file.exists(file_residency_skill) | overwrite) {
  
  residency_skill <- iteration[file.exists(file_residency), ]
  stopifnot(nrow(residency_skill) > 0L)
  
  residency_skill <- 
    split(residency_skill, seq_len(nrow(residency_skill))) |> 
    cl_lapply(function(d) {
      # d <- iteration_skill[1, ]
      # Load residency for simulated path 
      sim <- qs::qread(d$file_residency_sim)
      # Read output (estimated) residency
      out <- qs::qread(d$file_residency)
      # Merge datasets with "simulation" and "estimated" columns for residency statistics
      out <- 
        out |> 
        left_join(sim |> 
                    select("region", simulation = "estimate") |> 
                    as.data.table(),
                  by = "region")
      # Compute residency skill/error (estimated time - truth)
      # * <0: underestimation of residency
      # * >0: overestimation of residency  
      out[, skill := estimate - simulation]
      out
      
    }) |> 
    rbindlist() |> 
    # Add colours by region
    mutate(region = factor(region, levels = levels(regions_cs$region)), 
           col = regions_cs$col[match(region, regions_cs$region)]) |> 
    as.data.table()
  
  qs::qsave(residency_skill, file_residency_skill)
  residency_skill
  
} else {
  residency_skill <- qs::qread(file_residency_skill)
}
residency_skill[, perc := skill * 100]

#### Average residency skill over regions using MOE (Futia et al., 2024)
# MOE = mean(abs(truth_r - estimate_r) * truth_r), averaged over all regions r
# MOE varies between 0 (no error) and 100 % (completely wrong)
residency_skill_moe <- 
  residency_skill |> 
  group_by(individual_id, sensitivity) |> 
  mutate(algorithm = "patter") |> 
  # Compute MOE as the mean absolute error over all regions, weighted by prop. time per region
  mutate(moe = weighted.mean(abs(simulation - estimate), simulation) * 100) |> 
  # This code is equivalent: 
  # mutate(moe = sum(abs(simulation - estimate) * simulation) / sum(simulation) * 100) |> 
  slice(1L) |> 
  ungroup() |>
  select("individual_id", "algorithm", "sensitivity", "sensitivity_label", "moe") |> 
  as.data.table()


###########################
###########################
#### Example plots

#### Select individual 
it        <- iteration[1, ]
path      <- qs::qread(it$file_path)
occupancy <- terra::rast(it$file_occupancy)
residency <- residency_skill[individual_id == it$individual_id & sensitivity == "best", ]
stopifnot(it$sensitivity == "best")

#### Map simulated path for example individual
png(here_fig_sim("main", "example-path.png"), 
    height = 5, width = 5, units = "in", res = 800)
# Create base map with correct legend for path 
terra::plot(map,
            col = viridis::inferno(nrow(path)),
            range = range(path$timestep), 
            type = "continuous",
            legend = TRUE, 
            pax = list(labels = FALSE, lwd.ticks = 0))
# Cover base map colouration
terra::plot(map, col = "white", legend = FALSE, add = TRUE)
# Add path & coastline
patter:::add_sp_path(path$x, path$y, lwd = 0.75, length = 0, 
                     col = viridis::inferno(nrow(path))) |> 
  suppressWarnings()
# Add receivers
points(moorings$receiver_x, moorings$receiver_y, pch = 4, cex = 0.35)
terra::lines(champlain_utm, lwd = 0.5)
dev.off()

#### Map occupancy distribution for example individual
png(here_fig_sim("main", "example-occupancy.png"), 
    height = 5, width = 5, units = "in", res = 800)
terra::plot(occupancy, pax = list(labels = FALSE, lwd.ticks = 0))
points(moorings$receiver_x, moorings$receiver_y, pch = 4, cex = 0.35)
terra::lines(champlain_utm, lwd = 0.5)
dev.off()

#### Map occupancy quantiles for example individual
# Compute occupancy quantiles (as in analyse-patter-real.R)
probs              <- seq(0.05, 1, by = 0.05)
occupancy_contours <- terra::setValues(occupancy, 0)
for (p in probs) {
  occupancy_contours <- sum(occupancy_contours, 
                            map_hr_prop(occupancy, .prop = p), 
                            na.rm = TRUE)
}
occupancy_contours <- terra::classify(occupancy_contours, cbind(0, NA))
# Make map
png(here_fig_sim("main", "example-occupancy-quantiles.png"), 
    height = 5, width = 5, units = "in", res = 800)
terra::plot(occupancy_contours, pax = list(labels = FALSE, lwd.ticks = 0))
zlim <- c(0, length(probs))
terra::plot(occupancy_contours, 
            range = zlim, 
            col = grDevices::terrain.colors(length(probs), rev = TRUE),
            pax = list(labels = FALSE, lwd.ticks = 0),
            plg = list(at = pretty(zlim),
                       labels = prettyGraphics::add_lagging_point_zero(pretty(zlim) / max(zlim))))
points(moorings$receiver_x, moorings$receiver_y, pch = 4, cex = 0.35)
terra::lines(champlain_utm, lwd = 0.5)
dev.off()

#### Map residency error for example individual
# Update spatial layer with skill 
champlain_utm$skill <- residency$skill[match(champlain_utm$region, residency$region)]
champlain_utm$skill_perc <- champlain_utm$skill * 100
# Choose y limits to be symmetrical to force diverging colour scale centred at zero
mx <- max(abs(champlain_utm$skill_perc))
# Make map
png(here_fig_sim("main", "example-residency.png"), 
    height = 5, width = 5, units = "in", res = 800)
terra::plot(champlain_utm, y = "skill_perc",
            range = c(-mx, mx), 
            col = terra::map.pal("differences", 100), type = "continuous", 
            pax = list(labels = FALSE, lwd.ticks = 0), lwd = 0.5)
points(moorings$receiver_x, moorings$receiver_y, pch = 4, cex = 0.35)
dev.off()


###########################
###########################
#### Visualise skill

#### Choose individuals
# We plot maps for a few example individuals
# We choose maps without convergence failures
length(unique(iteration$sensitivity))
iteration[!file.exists(file_occupancy), .(index, individual_id, sensitivity)]
iteration |> 
  filter(file.exists(file_occupancy)) |> 
  group_by(individual_id) |> 
  filter(n() == 7L) |> 
  ungroup() |> 
  slice_sample(n = 3L) |> 
  pull(individual_id)
ids <- c(9, 10, 25)

#### Visualise occupancy skill (best maps)
# (rows: individuals; columns: trajectory, reconstructed map)
# (A) Define map_dt
columns <- c("A", "B") # simulation, output
map_dt <- lapply(columns, function(.column) {
  iteration |> 
    filter(individual_id %in% ids & sensitivity == "best") |> 
    mutate(file_ud = file_occupancy, 
           row = individual_id,
           column =  .column) |> 
    select("file_ud", "row", "column") |> 
    as.data.table()
}) |> rbindlist()
# (B) Define simulated paths
map_paths <- 
  paths |> 
  filter(path_id %in% map_dt$row) |> 
  mutate(row = path_id, 
         column = columns[1]) |> 
  select("row", "column", 
         "timestep", "x", "y") |> 
  as.data.table()
stopifnot(nrow(map_paths) > 1L)
# (C) Make maps
png(here_fig_sim("main", "maps-best.png"), 
    height = 5, width = 10, units = "in", res = 800)
p <- ggmaps(.mapdt = map_dt, .map = map, 
            .coast = champlain_utm, .geom_coast = list(colour = "black"),
            .path = map_paths, 
            .geom_path = list(linewidth = 0.2, arrow = grid::arrow(length = unit(0.001, "cm"))))
print(p)
dev.off()

#### Visualise occupancy skill (sensitivity maps)
# (Rows: individuals; columns: sensitivities)
# (A) Define map_dt
map_dt <- 
  iteration |> 
  filter(individual_id %in% ids) |> 
  mutate(file_ud = file_occupancy, 
         row = individual_id, 
         column = sensitivity_label) |> 
  select("file_ud", "row", "column") |> 
  as.data.table()
# (B) Make maps
png(here_fig_sim("main", "maps-sensitivity.png"), 
    height = 12, width = 12, units = "in", res = 800)
p <- ggmaps(.mapdt = map_dt, .map = map, 
            .coast = champlain_utm, .geom_coast = list(colour = "black"))
print(p)
dev.off()

#### As above but with trajectories & sensitivities on one plot
# First column: trajectories; then sensitivities
# (A) Define map_dt
columns <- c("Trajectory", levels(iteration$sensitivity_label))
columns <- factor(columns, levels = columns)
map_dt <- 
  rbind(
    # Trajectory column (1)
    iteration |> 
      filter(individual_id %in% ids & sensitivity == "best") |> 
      mutate(file_ud = file_occupancy, # here_input("map.tif")
             row = individual_id, 
             column = factor("Trajectory", levels = columns)) |> 
      select("file_ud", "row", "column"),
    # Sensitivity columns 
    iteration |> 
      filter(individual_id %in% ids) |> 
      mutate(file_ud = file_occupancy, 
             row = individual_id, 
             column = factor(sensitivity_label, levels = columns)) |> 
      select("file_ud", "row", "column") |> 
      as.data.table()
  ) |> 
  arrange(row, column) |> 
  mutate(file_ud = if_else(column == "Trajectory", NA_character_, file_ud)) |>
  as.data.table()
# (B) Define simulated paths
map_paths <- 
  paths |> 
  filter(path_id %in% map_dt$row) |> 
  mutate(row = path_id, 
         column = columns[1]) |> 
  select("row", "column", 
         "timestep", "x", "y") |> 
  as.data.table()
stopifnot(nrow(map_paths) > 1L)
# (C) Make maps
png(here_fig_sim("main", "maps-full.png"), 
    height = 14, width = 6, units = "in", res = 800)
p <- ggmaps(.mapdt = map_dt, .map = map, 
            .coast = champlain_utm, .geom_coast = list(colour = "black"),
            .path = map_paths, 
            .geom_path = list(linewidth = 0.2, arrow = grid::arrow(length = unit(0.001, "cm"))), 
            .scale_path = scale_colour_gradientn(colours = viridis::inferno(nrow(path)))) + 
  # Use bold
  theme(strip.text = element_text(face = "bold"))
print(p)
dev.off()

#### As above but using probability mass quantiles 
# Define probs and zlim (as above)
probs <- seq(0.05, 1, by = 0.05)
zlim  <- c(0, length(probs))
# Update map_dt with occupancy-quantiles.tif rasters
map_dt[, file_ud_quantile := file.path(dirname(file_ud), "occupancy-quantiles.tif")]
cl_lapply(split(map_dt, seq_len(nrow(map_dt))), function(d) {
  if (!is.na(d$file_ud)) {
    occupancy          <- terra::rast(d$file_ud)
    occupancy_contours <- terra::setValues(occupancy, 0)
    for (p in probs) {
      occupancy_contours <- sum(occupancy_contours, 
                                map_hr_prop(occupancy, .prop = p), 
                                na.rm = TRUE)
    }
    occupancy_contours <- terra::classify(occupancy_contours, cbind(0, NA))
    terra::writeRaster(occupancy_contours, d$file_ud_quantile, overwrite = TRUE)
  }
  nothing()
})
map_dt[, file_ud := file_ud_quantile]
map_dt[, file_ud_quantile := NULL]
# Make figure
op <- options(terra.pal = grDevices::terrain.colors(length(probs), rev = TRUE))
png(here_fig_sim("main", "maps-full-quantiles.png"), 
    height = 14, width = 6, units = "in", res = 800)
p <- ggmaps(.mapdt = map_dt, .map = map, .zlim = zlim,
            .coast = champlain_utm, .geom_coast = list(colour = "black"),
            .path = map_paths, 
            .geom_path = list(linewidth = 0.2, arrow = grid::arrow(length = unit(0.001, "cm"))), 
            .scale_path = scale_colour_gradientn(colours = viridis::inferno(nrow(path)))) + 
  # Use bold
  theme(strip.text = element_text(face = "bold"))
print(p)
dev.off()
options(op)

#### Visualise occupancy skill (MAE)
# This is principally useful as a measure of sensitivity
# It is hard to understand 'accuracy' from MAE
# But this provides a baseline for other methods 
png(here_fig_sim("main", "maps-skill.png"), 
    height = 5, width = 8, units = "in", res = 800)
p <- 
  occupancy_skill |> 
  ggplot() + 
  geom_violin(aes(sensitivity_label, skill, fill = sensitivity_label), 
              linewidth = 0.25, 
              scale = "count") +
  scale_y_continuous(labels = prettyGraphics::sci_notation) + 
  xlab("Sensitivity") + 
  ylab(expression("Mean absolute error (" * italic(ME) * ")")) + 
  labs(fill = "Analysis") +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) 
print(p)
dev.off()

#### Visualise residency skill by region, for 'best' analyses
png(here_fig_sim("main", "residency-skill-best.png"), 
    height = 5, width = 8, units = "in", res = 800)
p <- 
  residency_skill |>
  filter(sensitivity == "best") |> 
  as_tibble() |> 
  ggplot() + 
  geom_boxplot(aes(region, perc, fill = I(col)), 
               linewidth = 0.5, size = 1, varwidth = TRUE) + 
  geom_hline(yintercept = 0, linetype = 3) + 
  # scale_y_continuous(expand = c(0, 0), limits = c(-1, 1)) + 
  xlab("Region") + 
  ylab("Residency error (%)") + 
  labs(fill = "Region") +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)), 
        axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()

#### Visualise residency skill, by region, including sensitivity
# This plot is by region & sensitivity
png(here_fig_sim("main", "residency-skill-sensitivity.png"), 
    height = 6, width = 12, units = "in", res = 800)
p <- 
  residency_skill |>
  ggplot() + 
  geom_boxplot(aes(region, perc, fill = sensitivity_label), 
               linewidth = 0.25, size = 0.5, varwidth = TRUE) + 
  geom_hline(yintercept = 0, linetype = 3) + 
  # scale_y_continuous(expand = c(0, 0), limits = c(-1, 1)) + 
  xlab("Region") + 
  ylab("Residency error (%)") + 
  labs(fill = "Analysis") +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)), 
        axis.text.x = element_text(angle = 45, hjust = 1)) 
print(p)
dev.off()

#### Summarise residency skill overs region simply (%)
# Range in absolute mean error for 'best' analyses
residency_skill |> 
  filter(sensitivity == "best") |> 
  summarise(utils.add::basic_stats(abs(simulation - estimate) * 100))
# Overall residency skill for 'best' analyses
residency_skill |> 
  filter(sensitivity == "best") |> 
  # filter(!(estimate == 0 & simulation == 0)) |> 
  reframe(utils.add::basic_stats(perc))
# Residency skill for 'best' analyses split by region
residency_skill |> 
  filter(sensitivity == "best") |> 
  group_by(region) |> 
  # filter(!(estimate == 0 & simulation == 0)) |> 
  reframe(utils.add::basic_stats(perc))

#### Visualise residency skill averaged over tracks using MOE, including sensitivity
# Plot distribution of MOE % across all tracks by sensitivity, as in Futia et al. 2024
# * Each point is the MOE for a single track
# * We expect similar median MOE compared to best model (4 %)
# * We expect reduced maximum MOE (from 40.4% -> 10 %)
# * I.e., we expect reduced variation (increased precision)
# * We show the variation with a boxplot over all tracks
png(here_fig_sim("main", "residency-skill-moe.png"), 
    height = 4, width = 8, units = "in", res = 800)
p <- 
  residency_skill_moe |>
  ggplot(aes(sensitivity_label, moe, fill = sensitivity_label)) + 
  geom_boxplot(linewidth = 0.25, size = 0.5, varwidth = TRUE, outliers = FALSE) + 
  geom_jitter(size = 0.25, colour = "dimgrey", width = 0.1, height = 0) +
  # scale_y_continuous(expand = c(0, 0), limits = c(-1, 1)) + 
  xlab("Sensitivity") + 
  ylab("MOE (%)") + 
  labs(fill = "Analysis") +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)), 
        axis.text.x = element_text(angle = 45, hjust = 1)) 
print(p)
dev.off()

#### (optional) Compare MOE to heuristic methods (Futia et al., 2024)

## Read MOE scores for best model in Futia et al. (2024)
# * AInt_model_performance.qs includes error by region for each track 
# * Int_model_moe.qs has the mean weighted error (sum across regions by individuals) that was used for Fig 3.
residency_skill_moe_int <- 
  here_data_raw_mf("Int_model_moe.qs") |> 
  qs::qread() |> 
  mutate(
    individual_id = animal_id, 
    algorithm = as.character(model), 
    sensitivity = "Int",
    sensitivity_label = factor("Int", levels = c("Int", levels(residency_skill_moe$sensitivity_label))), 
    moe = mean_wt_abs_err
  ) |> 
  select("individual_id", "algorithm", "sensitivity", "sensitivity_label", "moe") |> 
  as.data.table()

## Compute summary statistics
# For patter:
# * TO DO
# For Int algorithm: 
# * Average MOE ± SD = 5.1 ± 8.1 %
# * Max MOE =40.4 % 
residency_skill_moe_full <- rbind(residency_skill_moe, residency_skill_moe_int)
residency_skill_moe_full |> 
  group_by(algorithm, sensitivity) |> 
  reframe(utils.add::basic_stats(moe))

#### (optional) Update ggplot of MOE including Int model from Futia et al. (2024)
png(here_fig_sim("main", "residency-skill-moe-full.png"),
    height = 6, width = 12, units = "in", res = 800)
p <-
  residency_skill_moe_full |>
  ggplot(aes(sensitivity_label, moe, fill = sensitivity_label)) +
  geom_boxplot(linewidth = 0.25, size = 0.5, varwidth = TRUE, outliers = FALSE) +
  geom_jitter(size = 0.25, colour = "dimgrey", width = 0.1, height = 0) +
  # scale_y_continuous(expand = c(0, 0), limits = c(-1, 1)) +
  xlab("Sensitivity") +
  ylab("MOE (%)") +
  labs(fill = "Analysis") +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()


#### End of code. 
###########################
###########################