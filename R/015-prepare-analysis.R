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
library(truncdist)
files_source_r(here_src())

#### Load data
pars <- qs::qread(here_input("pars-patter.qs"))


###########################
###########################
#### Select analysis

#### Define analysis 
# analysis <- "sim"
analysis <- "real"
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
  select(unit_id, individual_id, time_id, timestamp) |>
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
# A) Assign unit_id from unitsets
detections <- 
  detections |> 
  left_join(unitsets |> 
              select(unit_id, individual_id, time_id) |>
              as.data.table(), 
            by = c("individual_id", "time_id")) |> 
  as.data.table()
# B) Filter detections 
detections <- filter_detections(detections)
# C) Update unitsets
unitsets <- unitsets[unit_id %in% detections$unit_id, ]

#### Checks
# Visually validate matching between unitsets & detections
unitsets[unit_id == 14, ]
detections[unit_id == 14, ]
# Validate all unit_ids present in each dataset
stopifnot(all(unitsets$unit_id %in% detections$unit_id) & 
            all(detections$unit_id %in% unitsets$unit_id))
# Validate matching between all unitsets and detections
cl_lapply(split(unitsets, seq_len(nrow(unitsets))), function(sim) {
  vdetections <- detections[unit_id == sim$unit_id, ]
  stopifnot(all(sim$unit_id == vdetections$unit_id))
  stopifnot(all(sim$individual_id == vdetections$individual_id))
  stopifnot(all(sim$time_id == vdetections$time_id))
})
# Validate time_id assignment in detections
stopifnot(all(detections$time_id == 
                lubridate::floor_date(detections$timestamp, "months")))
# Check the number of days with detections meets minimum criteria
# * This is based on the threshold specified in filter_detections.R
ck <- 
  detections |> 
  group_by(unit_id) |> 
  summarise(ck = length(unique(lubridate::yday(timestamp)))) |> 
  pull(ck) |> 
  sort()
stopifnot(all(ck >= 14))
# Check the number of detections
# * We want to catch time series with 'too few' observations
# * What is 'too few' here is somewhat arbitrary
# * The goal is to catch potential mistakes in data processing
# * E.g., that would otherwise allow 1 or 2 row datasets forward for analysis
ck <- 
  detections |> 
  group_by(unit_id) |> 
  summarise(ck = n()) |> 
  pull(ck) |> 
  sort()
stopifnot(all(ck >= 50))

#### Summarise detection dataset
# cf. raw data summary statistics (setup-data-detection.R)
nrow(detections)
length(unique(detections$individual_id))
range(detections$timestamp)
int <- lubridate::interval(min(detections$timestamp),max(detections$timestamp))
lubridate::time_length(int, "months")
lubridate::time_length(int, "years")
length(unique(detections$receiver_id))
length(unique(paste(detections$individual_id, detections$time_id)))

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

#### Validation
# Validate match between iteration and file_detections
cl_lapply(split(iteration, seq_len(nrow(iteration))), function(sim) {
  # sim <- iteration[1, ]
  # print(sim$index)
  vdetections <- qs::qread(sim$file_detections)
  stopifnot(all(sim$unit_id == sim$unit_id))
  stopifnot(all(sim$individual_id == vdetections$individual_id))
  stopifnot(all(sim$time_id == vdetections$time_id))
  stopifnot(all(vdetections$time_id == lubridate::floor_date(vdetections$timestamp, "months")))
})


###########################
###########################
#### Prepare iteration patter: optimisation analysis

# We will explore estimation of latent locations & static parameters
# (focusing on shape/scale parameters of gamma distribution)

# We explore two optimisation routines:
# A) optim()
# - programmatically quick & easy to extend for multiple parameters
# - but did not work not work well in initial tests 
# - see analysis-optimisation.R
# * grid-search
# - scales poorly with increasing numbers of parameters
# - but parallelisable

# We'll run optimisation for the following settings:
# * We consider a subset of N individuals
# * For each individual, we run the filter/optimisation 3 times 
# * We'll record outputs in:
# * sim/optim/individual_id/rep_id/parameter_id (parameter_id = 1)
# * sim/grid/individual_id/rep_id/parameter_id  (multiple parameters)

if (analysis == "sim") {
  
  #### Copy iterations
  iteration_main <- copy(iteration)
  
  #### Define parameters
  n_id  <- 3L
  n_rep <- 3L
  
  #### Explore possible step length distributions 
  CJ(shape = seq(1, 10, by = 2),
     scale = seq(20, 30, by = 1)) |>
    mutate(row = paste(shape, scale, sep = ", ")) |> 
    tidyr::expand_grid(x = seq(0, pars$mobility[1], length.out = 100)) |> 
    ggplot(aes(x, dgamma(x, shape = shape, scale = scale))) +
    geom_line() +
    facet_wrap(~row, scales = "free_y") + 
    theme(axis.text.y = element_blank())
  # Examine sampling distributions
  hist(rtrunc(100, "norm", lower = 0, mean = pars$shape[1], sd = 2))
  hist(rtrunc(100, "norm", lower = 0, mean = pars$scale[1], sd = 2))
  
  
  ###########################
  #### optim analysis
  
  #### Define iteration
  iteration <- 
    iteration_main |> 
    filter(sensitivity == "best") |> 
    slice(1:n_id) |> 
    cross_join(data.table(rep_id = 1:n_rep)) |> 
    mutate(index = row_number(), 
           parameter_id = 1L, 
           # Simulate starting values for shape/scale for optimisation
           shape = rtrunc(n(), "norm", lower = 0, mean = pars$shape[1], sd = 2),
           scale = rtrunc(n(), "norm", lower = 0, mean = pars$scale[1], sd = 2),
           file_output = file.path("data", "output", analysis, "optim", "runs", 
                                   individual_id, rep_id, parameter_id, "optim.qs")) |> 
    # Select columns, including parameters required by constructor_ac_core()
    select(index, unit_id, individual_id, rep_id, 
           shape, scale, phi, mobility, 
           receiver_alpha, receiver_beta, receiver_gamma, 
           file_detections, file_output) |> 
    as.data.table()
  
  # Examine starting distributions for optimisation
  iteration |>
    select(shape, scale) |> 
    mutate(row = paste0(shape, scale, ", ")) |> 
    tidyr::expand_grid(x = seq(0, pars$mobility[1], length.out = 200)) |> 
    ggplot(aes(x, dgamma(x, shape = shape, scale = scale), 
               colour = row, group = row)) +
    geom_line() +
    theme(axis.text.y = element_blank(), 
          legend.position = "none")
  
  # Record iteration
  qs::qsave(iteration, here_input_analysis("iteration-patter-optim.qs"))
  
  # Build directories
  dirs.create(dirname(iteration$file_output))
  
  
  ###########################
  #### Grid-search

  # Define parameter grid
  # * We know the true parameter values
  # * We could consider the same bounds of uncertainty as in real-world analyses
  # * This is relatively well defined 
  shapes <- sort(c(pars$shape[1], seq(min(pars$shape), max(pars$shape), length.out = 10)))
  scales <- sort(c(pars$scale[1], seq(min(pars$scale), max(pars$scale), length.out = 10)))
  grid   <- CJ(shape = shapes, scale = scales) |> 
    mutate(parameter_id = row_number(), 
           row =  paste(shape, scale, sep = ", ")) |>
    select(parameter_id, row, shape, scale) |>
    as.data.table()
  # Visualise parameter grid
  grid |>
    tidyr::expand_grid(x = seq(0, max(pars$mobility[1]), length.out = 100)) |> 
    ggplot(aes(x, dgamma(x, shape = shape, scale = scale), 
               colour = row, group = row)) +
    geom_line() +
    theme(axis.text.y = element_blank(), 
          legend.position = "none")
  
  # Define iteration data.table
  iteration <- 
    iteration_main |> 
    filter(sensitivity == "best") |> 
    select(-parameter_id, -shape, -scale) |> 
    slice(1:n_id) |> 
    cross_join(data.table(rep_id = 1:n_rep)) |> 
    cross_join(grid) |> 
    mutate(index = row_number()) |> 
    mutate(file_output = file.path("data", "output", analysis, "grid", "runs", 
                                   individual_id, rep_id, parameter_id, "grid.qs")) |> 
    # Select columns, including parameters required by constructor_ac_core()
    select(index, unit_id, individual_id, rep_id, parameter_id, 
           shape, scale, mobility, phi, 
           receiver_alpha, receiver_beta, receiver_gamma,
           file_detections, file_output) |> 
    as.data.table()
  
  # Record iteration
  qs::qsave(iteration, here_input_analysis("iteration-patter-grid.qs"))
  
  # Build directories
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