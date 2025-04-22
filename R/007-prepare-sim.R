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
pars_model_move <- qs::qread(here_input("pars-model-move-full.qs"))
pars_model_obs  <- qs::qread(here_input("pars-model-obs-full.qs"))


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
# Check formatting
stopifnot(nrow(pars_model_move) == 3L)
stopifnot(nrow(pars_model_obs) == 1L)
# Define parameters 
pars <- 
  rbind(
    cbind(sensitivity = "best", pars_model_move[1, ], pars_model_obs[1, ]),
    cbind(sensitivity = "move(-)", pars_model_move[2, ], pars_model_obs[1, ]),
    cbind(sensitivity = "move(+)", pars_model_move[3, ], pars_model_obs[1, ])
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