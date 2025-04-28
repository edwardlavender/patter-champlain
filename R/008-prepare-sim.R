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
set.seed(123L)

#### Load essential packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(proj.verse)
files_source_r(here_src())

#### Load data
pars <- qs::qread(here_input("pars-patter.qs"))


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
    file_detections = file.path("data", "input", "sim", 
                                individual_id, time_id, "detection.qs"),
    folder_home = file.path("data", "output", "sim", "main", "runs", 
                            individual_id, time_id), 
    folder_home_patter = file.path(folder_home, "patter")) |>
  as.data.table()

#### Build directories
dirs.create(dirname(unitsets$file_detections))
dirs.create(unitsets$folder_home)
dirs.create(unitsets$folder_home_patter)

#### Write to file
qs::qsave(unitsets, here_input_sim("unitsets.qs"))


###########################
###########################
#### Prepare iteration patter: main analysis

#### Define iteration 
iteration <- 
  unitsets |> 
  select(unit_id, individual_id, time_id, 
         file_detections,
         folder_home = folder_home_patter) |>
  cross_join(pars) |> 
  mutate(index = row_number(),
         folder_coord = file.path(folder_home, "coord", parameter_id)) |> 
  as.data.table()

#### Build directories 
if (FALSE) {
  unlink(iteration$folder_coord)
}
nrow(iteration)
dirs.create(iteration$folder_coord)

#### Record iteration
qs::qsave(iteration, here_input_sim("iteration-patter.qs"))


###########################
###########################
#### Prepare iteration patter: optimisation analysis

#### Define iteration
# We consider a subset of N individuals
# For each individual, we run the filter/optimisation 3 times 
# We save optim() outputs in sim/optimisation/individual_id/rep_id/
n_id  <- 3L
n_rep <- 3L
iteration <- 
  iteration |> 
  filter(sensitivity == "best") |> 
  slice(1:n_id) |> 
  cross_join(data.table(rep_id = 1:n_rep)) |> 
  mutate(index = row_number()) |> 
  mutate(file_output = file.path("data", "output", "sim", "optim", "runs", 
                                 individual_id, rep_id, "optim.qs")) |> 
  # Select columns, including parameters required by constructor_ac_sim()
  select(index, unit_id, individual_id, rep_id, 
         phi, mobility, receiver_alpha, receiver_beta, receiver_gamma, , 
         file_detections, file_output) |> 
  as.data.table()

#### Update iteration with initial parameter values
# Explore possible step length distributions 
CJ(shape = seq(1, 10, by = 2),
   scale = seq(20, 30, by = 1)) |>
  mutate(row = paste(shape, scale, sep = ", ")) |> 
  tidyr::expand_grid(x = seq(0, pars$mobility[1], length.out = 100)) |> 
  ggplot(aes(x, dgamma(x, shape = shape, scale = scale))) +
  geom_line() +
  facet_wrap(~row, scales = "free_y")
# Select suitable parameters for sampling distribution 
hist(rnorm(100, mean = pars$shape[1], sd = 2))
hist(rnorm(100, mean = pars$scale[1], sd = 2))
# Simulate init parameters 
iteration[, shape := rnorm(.N, mean = pars$shape[1], sd = 2)]
iteration[, scale := rnorm(.N, mean = pars$scale[1], sd = 2)]
stopifnot(all(iteration$shape > 0) & all(iteration$scale > 0))
# Examine distributions
iteration |>
  select(shape, scale) |> 
  mutate(row = paste0(shape, scale, ", ")) |> 
  tidyr::expand_grid(x = seq(0, pars$mobility[1], length.out = 200)) |> 
  ggplot(aes(x, dgamma(x, shape = shape, scale = scale))) +
  geom_line() +
  facet_wrap(~row, scales = "free_y")
# Record iteration
qs::qsave(iteration, here_input_sim("iteration-patter-optim.qs"))

#### Build directories
dirs.create(dirname(iteration$file_output))


###########################
###########################
#### Prepare iteration heuristics

# TO DO (MF)


#### End of code. 
###########################
###########################