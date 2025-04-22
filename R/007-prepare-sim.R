###########################
###########################
#### prepare-sim.R

#### Aims
# 1) Prepares an iteration data.table and directories for simulations

#### Prerequisites
# 1) Previous scripts
# 2) Following ?patter.workflows::`patter.workflows-package`


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
pars_model_move <- qs::qread(here_input("pars-model-move.qs"))
pars_model_obs  <- qs::qread(here_input("pars-model-obs.qs"))


###########################
###########################
#### Define unitsets 

#### Number of simulations (individuals)
n_sim <- 30L

#### Define unitsets
# For consistency, we include the same columns as in the real-world analysis
unitsets <-
  data.table(unit_id = 1:30L, 
             individual_id = 1:n_sim, 
             time_id = 1:n_sim) |> 
  mutate(
    file_detection = file.path("data", "input", "sim", individual_id, time_id, "detection.qs"),
    folder_home = file.path("data", "output", "sim", "runs", individual_id, time_id), 
    folder_home_patter = file.path(folder_home, "patter")) |>
  as.data.table()

#### Build directories
dirs.create(dirname(unitsets$file_detection))
dirs.create(unitsets$folder_home)
dirs.create(unitsets$folder_home_patter)

#### Write to file
qs::qsave(unitsets, here_input_sim("unitsets.qs"))


###########################
###########################
#### Prepare iteration patter

#### Define parameter combinations

## (A) Movement parameters
# Define parameter uncertainty 
adj <- 0.5
inflate <- 1 + adj
deflate <- 1 - adj
# Collect 'best-guess' parameters
mobility <- pars_model_move$mobility
shape    <- pars_model_move$shape
scale    <- pars_model_move$scale
phi      <- pars_model_move$phi
# Define more restrictive/flexible models
restrictive <- gamma_rescale(shape, scale, fact = deflate)
flexible    <- gamma_rescale(shape, scale, fact = inflate)
# Collect movement parameters (best, restrictive, flexible)
pars_movement <- data.table(mobility = c(mobility, mobility * deflate, mobility * inflate),
                              shape = c(shape, restrictive[1], flexible[1]),
                              scale = c(scale,  restrictive[2],flexible[2]),
                              phi = c(phi, phi * deflate, phi * inflate))
# Visualise models
# Best model (step-length)
plot_dbn("gamma", 
         xlim = c(0, 300), 
         pars = list(shape = pars_movement$shape[1], scale = pars_movement$scale[1]))
# Restrictive model (step-length)
plot_dbn("gamma", 
         xlim = c(0, 300), 
         pars = list(shape = pars_movement$shape[2], scale = pars_movement$scale[2]), 
         add = TRUE, col = "red")
# Flexible model (step-length)
plot_dbn("gamma", 
         xlim = c(0, 300), 
         pars = list(shape = pars_movement$shape[3], scale = pars_movement$scale[3]), 
         add = TRUE, col = "darkgreen")

## (B) Detection parameters
# We assume these are known
# To minimise computation time, we only explore the effects of uncertainty in movement
# We find this more interesting
pars_detection <- as.data.table(pars_model_obs)

## (C) Collect parameter combinations
pars <- 
  rbind(
    cbind(sensitivity = "best", pars_movement[1, ], pars_detection[1, ]),
    cbind(sensitivity = "move(-)", pars_movement[2, ], pars_detection[1, ]),
    cbind(sensitivity = "move(+)", pars_movement[3, ], pars_detection[1, ])
  ) |> 
  mutate(parameter_id = row_number()) |> 
  as.data.table()

#### Define iteration 
iteration <- 
  unitsets |> 
  select(unit_id, individual_id, time_id, folder_home = folder_home_patter) |>
  cross_join(pars) |> 
  mutate(index = row_number(),
         folder_coord = file.path(folder_home, "coord", parameter_id)) |> 
  as.data.table()

#### Build directories 
nrow(iteration)
dirs.create(iteration$folder_coord)

#### Record iteration
qs::qsave(iteration, here_input_sim("iteration-patter.qs"))


###########################
###########################
#### Prepare iteration heuristics

# TO DO


#### End of code. 
###########################
###########################