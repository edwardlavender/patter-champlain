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
library(proj.verse)
library(tictoc)
files_source_r(here_src())

#### Load data
map       <- terra::rast(here_input("map.tif"))
iteration <- qs::qread(here_input_sim("main", "iteration.qs"))


###########################
###########################
#### Compute skill metrics 

#### Compute occupancy skill
# Compute mean absolute error between simulated & reconstructed patterns of space use
overwrite <- TRUE
file_occupancy_skill <- here_output_sim_main("synthesis", "occupancy-skill.qs")
if (!file.exists(file_occupancy_skill) | overwrite) {
  
  occupancy_skill <- iteration[file.exists(file_occupancy), ]
  skills <- 
    split(occupancy_skill, seq_len(nrow(occupancy_skill))) |> 
    cl_lapply(function(d) {
      skill_me(.obs = terra::rast(d$file_occupancy_sim), 
               .mod = terra::rast(d$file_occupancy)) 
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
    rbindlist()
  
  qs::qsave(residency_skill, file_residency_skill)
  
} else {
  residency_skill <- qs::qread(file_residency_skill)
}


###########################
###########################
#### Visualise skill

#### Visualise occupancy skill (full)
# TO DO CLEAN UP DRAFT CODE
occupancy_skill |> 
  ggplot() + 
  geom_violin(aes(sensitivity_label, skill, fill = sensitivity_label), 
              linewidth = 0.25, 
              scale = "count") +
  # scale_y_continuous(
  #   limits = c(0, 1.7e-5),
  #   expand = c(0, 0),
  #   labels = prettyGraphics::sci_notation) + 
  xlab("Sensitivity") + 
  ylab(expression("Mean absolute error (" * italic(ME) * ")")) + 
  labs(fill = "Sensitivity") +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) 

#### Visualise residency skill (full)
# This plot is by region & sensitivity
# TO DO CLEAN UP DRAFT CODE
residency_skill |>
  ggplot() + 
  geom_violin(aes(sensitivity_label, skill, fill = region), 
              linewidth = 0.25,
              scale = "count") + 
  geom_hline(yintercept = 0, linetype = 3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-1, 1)) + 
  xlab("Sensitivity") + 
  ylab(expression("Residency error (" * italic(RE) * ")")) + 
  labs(fill = "Region") +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) 


#### End of code. 
###########################
###########################