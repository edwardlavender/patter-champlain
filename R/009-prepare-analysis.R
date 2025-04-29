###########################
###########################
#### prepare-analysis.R

#### Aims
# 1) Prepares detection datasets, iteration data.tables and directories for analyses
#    (This includes analyses of both simulation and real-world datasets)

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
#### Select analysis

#### Define analysis 
# analysis <- "sim"
# analysis <- "real"
stopifnot(analysis %in% c("sim", "real"))

#### Define analysis-specific routines
here_input_analysis <- switch_here_input_analysis(analysis)

#### Define analysis-specific data
detections <- qs::qread(here_input_analysis("detections.qs"))


###########################
###########################
#### Define unitsets

#### Define units (individual/month combinations)
detections <- 
  detections |> 
  group_by(individual_id) |> 
  mutate(time_id = lubridate::floor_date(timestamp, "months")) |>
  select(individual_id, time_id, timestamp, receiver_id) |>
  as.data.table()

#### Define unitsets
unitsets <- 
  detections |> 
  group_by(individual_id, time_id) |> 
  slice(1L) |> 
  ungroup() |> 
  mutate(unit_id = row_number()) |> 
  select(unit_id, individual_id, time_id, timestamp, receiver_id) |>
  mutate(
    file_detections = file.path("data", "input", analysis, 
                                individual_id, time_id, "detection.qs"),
    folder_home = file.path("data", "output", analysis, "main", "runs", 
                            individual_id, time_id), 
    folder_home_patter = file.path(folder_home, "patter")) |>
  as.data.table()

#### Build directories
if (FALSE) {
  unlink(file.path("data", "output", analysis, "main", "runs"))
}
dirs.create(dirname(unitsets$file_detections))
dirs.create(unitsets$folder_home)
dirs.create(unitsets$folder_home_patter)


###########################
###########################
#### Filter dataset

#### Number of individuals in full dataset
length(unique(detections$individual_id))

#### (optional) Focus on individuals with sufficient data
# This is no longer implemented
# We include all individuals
# But focus on individual/month blocks that meet selected criteria (below)

#### Focus on individual/month units with sufficient data
# Check the number of unit_ids in the raw data:
detections[, unit_id := .GRP, by = c("individual_id", "time_id")]
detections <- filter_detections(detections)
unitsets   <- unitsets[unit_id %in% detections$unit_id, ]

#### Write unitsets/detections
qs::qsave(unitsets, here_input_analysis("unitsets.qs"))
detections[, file_detections := unitsets$file_detections[match(unit_id, unitsets$unit_id)]]
cl_lapply(split(detections, detections$unit_id), function(d) {
  qs::qsave(d, d$file_detections[1])
})


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

#### Check nrows
# TO DO Limit iteration rows
# * For the real analysis, the number of rows in this data.table is too high
# * We should restrict this
#   - A) Restrict individual inclusion criteria 
#   - B) Reconsider simulation priorities 
#   - C) Restrict sensitivity analyses
#        on the basis of simulation results the sensitivity analysis 
nrow(iteration)

#### Build directories 
if (FALSE) {
  unlink(iteration$folder_coord)
}
nrow(iteration)
dirs.create(iteration$folder_coord)

#### Write 
qs::qsave(iteration, here_input_analysis("iteration-patter.qs"))

###########################
###########################
#### Prepare iteration patter: optimisation analysis

if (analysis == "sim") {
  
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
    mutate(file_output = file.path("data", "output", analysis, "optim", "runs", 
                                   individual_id, rep_id, "optim.qs")) |> 
    # Select columns, including parameters required by constructor_ac_core()
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
  qs::qsave(iteration, here_input_analysis("iteration-patter-optim.qs"))
  
  #### Build directories
  dirs.create(dirname(iteration$file_output))
  
}


###########################
###########################
#### Prepare iteration heuristics

# TO DO (MF)


###########################
###########################
#### Write outputs




#### End of code. 
###########################
###########################