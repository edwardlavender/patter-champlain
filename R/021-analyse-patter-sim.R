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


###########################
###########################
#### Compute skill metrics 

#### Compute ncell
nc <- terra::freq(map)[["count"]]

#### Compute occupancy skill
# Compute mean absolute error between simulated & reconstructed patterns of space use
# > This is primarily used to understand sensitivity
# > It is useful for future comparisons against other methods
# > Other metrics e.g., EMD are more intepretable in terms of how 'accurate' maps are
# > But for this work we simply compare simulated tracks & associated maps
overwrite <- TRUE
file_occupancy_skill <- here_output_sim_main("synthesis", "occupancy-skill.qs")
if (!file.exists(file_occupancy_skill) | overwrite) {
  
  occupancy_skill <- iteration[file.exists(file_occupancy), ]
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
  
} else {
  occupancy_skill <- qs::qread(file_occupancy_skill)
}

#### Compute residency skill
# Compute error between simulated & reconstructed residency estimates _by region_
overwrite <- TRUE
file_residency_skill <- here_output_sim_main("synthesis", "residency-skill.qs")
if (!file.exists(file_residency_skill) | overwrite) {
  
  residency_skill <- iteration[file.exists(file_residency), ]
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
  
} else {
  residency_skill <- qs::qread(file_residency_skill)
}


###########################
###########################
#### Visualise skill

#### Visualise occupancy skill (best maps)
# (rows: individuals; columns: trajectory, reconstructed map)
# (A) Define map_dt
columns <- c("A", "B") # simulation, output
map_dt <- lapply(columns, function(.column) {
  iteration |> 
    filter(individual_id %in% 1:3L & sensitivity == "best") |> 
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
  filter(individual_id %in% 1:3L) |> 
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

#### Visualise occupancy skill (full)
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
  labs(fill = "Sensitivity") +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) 
print(p)
dev.off()

#### Visualise residency skill (error by region, for 'best' analyses)
png(here_fig_sim("main", "residency-skill-best.png"), 
    height = 5, width = 8, units = "in", res = 800)
p <- 
  residency_skill |>
  filter(sensitivity == "best") |> 
  as_tibble() |> 
  ggplot() + 
  geom_boxplot(aes(region, skill, fill = I(col)), 
               linewidth = 0.5, size = 1) + 
  geom_hline(yintercept = 0, linetype = 3) + 
  # scale_y_continuous(expand = c(0, 0), limits = c(-1, 1)) + 
  xlab("Region") + 
  ylab(expression("Residency error (" * italic(RE) * ")")) + 
  labs(fill = "Region") +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)), 
        axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()

#### Visualise residency skill (full)
# This plot is by region & sensitivity
png(here_fig_sim("main", "residency-skill-sensitivity.png"), 
    height = 6, width = 12, units = "in", res = 800)
p <- 
  residency_skill |>
  ggplot() + 
  geom_boxplot(aes(region, skill, fill = sensitivity_label), 
               linewidth = 0.25, size = 0.5) + 
  geom_hline(yintercept = 0, linetype = 3) + 
  # scale_y_continuous(expand = c(0, 0), limits = c(-1, 1)) + 
  xlab("Region") + 
  ylab(expression("Residency error (" * italic(RE) * ")")) + 
  labs(fill = "Sensitivity") +
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