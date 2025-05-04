###########################
###########################
#### analysis-optimisation.R

#### Aims
# 1) Analyses simulated observations using patter
#    specifically testing how well we can estimate latent locations & parameters

#### Prerequisites
# 1) Run simulations


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Set global options
patter::julia_connect()

#### Load essential packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(patter)
library(patter.workflows)
library(proj.verse)
files_source_r(here_src())

#### Load data
pars_model_move <- qs::qread(here_input("pars-model-move-best.qs"))
iteration_optim <- qs::qread(here_input_sim("iteration-patter-optim.qs"))
iteration_grid  <- qs::qread(here_input_sim("iteration-patter-grid.qs"))


###########################
###########################
#### Set up analysis 

#### Define iteration 
# iteration <- iteration_optim
iteration <- iteration_grid

#### (optional) Clean up
if (FALSE) {
  unlink(iteration$file_output)
}

#### (optional) Set up cluster
# TO DO Set up cluster
map <- here_input("map.tif")
set_map(map)


###########################
###########################
#### Run optimisation via optim()

if (FALSE) {
  
  #### Run filter & optimisation
  # Iterate over iteration and run optimisation
  # We assume the function can be initialised at the initial parameters in this data.table!
  # For optimisation, we try optim()
  # This is more scalable for multiple parameters
  # With just two parameters, a grid search approach may be quicker
  cl_lapply_workflow(.iteration = iteration[1, ], 
                     .datasets = NULL, 
                     .constructor = constructor_pf_filter_loglik_optim, 
                     .algorithm = pf_filter_loglik_optim)
  
  #### Results for iteration[1, ]
  
  # Time:
  # qs::qread(file.path(dirname(iteration$file_output[1]), "callstats.qs"))
  # (10 julia threads, 2.5e4 particles, 1e3 smoothing particles)
  
  # Estimates:
  # qs::qread(iteration$file_output[1])
  # $par
  # [1]  7.297578 25.300327
  # $value
  # [1] -8455.062
  # $counts
  # function gradient 
  # 305       NA 
  # $convergence
  # [1] 10 # degeneracy of the Nelder–Mead simplex
  # $message
  # NULL
  # $hessian
  # [,1]     [,2]
  # [1,] -213315.39 22294.09
  # [2,]   22294.09 37835.44
  
  # Visually compare 'true' versus estimated parameters in this example:
  plot_dbn("gamma", lower = 0, upper = 216, xlim = c(0, 216), 
           pars = list(shape = 3.25, scale = 25))
  plot_dbn("gamma", lower = 0, upper = 216, xlim = c(0, 216), 
           pars = list(shape = 7.297578, scale = 25.300327), 
           add = TRUE, col = "red")
  
  #### Preliminary conclusions
  # optim doesn't work well if we just use default settings 
  # (For optim, careful customisation may be required)
  # We'll try a grid search instead
  # We'll leverage biological knowledge to search within a sensible parameter space
  
}


###########################
###########################
#### Grid-based search

if (TRUE) {
  
  # Prepare iteration
  iteration <- iteration[individual_id == 1L, ]
  iteration[, smooth := FALSE]
  nrow(iteration)
  
  # Run grid search
  # * Use standard constructor_ac_sim constructor with iteration$smooth = FALSE
  cl_lapply_workflow(.iteration   = iteration[1:2, ], 
                     .datasets    = NULL, 
                     .constructor = constructor_ac_sim, 
                     .algorithm   = estimate_coord_particle, 
                     .success     = function(x) x$forward$callstats$convergence, 
                     .cleanup     = particle_cleanup)
  
  # Collate log-likelihood
  iteration_ll <- iteration[file.exists(file_output), ]
  logliks      <- cl_lapply_iteration_file(
    iteration_ll, 
    .file = "file_output", 
    .fun = function(.sim, .input) {
      .input$forward$callstats$loglik
    }) |> unlist()
  iteration_ll[, loglik := logliks]
  
  # Plot marginal log likelihood profile(s)
  # A) For each individual/rep, for shape plot max(loglik over all scales) ~ shape value
  pshape <-
    iteration_ll |> 
    group_by(individual_id, rep_id, shape) |> 
    summarise(loglik = max(loglik)) |>
    ungroup() |>
    as.data.frame() |>
    ggplot(aes(x = shape, y = loglik)) + 
    geom_point() + 
    facet_wrap(~individual_id)
  # B) As above for scale
  pscale <- 
    iteration_ll |> 
    group_by(individual_id, rep_id, scale) |> 
    summarise(loglik = max(loglik)) |>
    ungroup() |>
    as.data.frame() |>
    ggplot(aes(x = scale, y = loglik)) + 
    geom_point() + 
    facet_wrap(~individual_id)
  # C) Plots 
  gridExtra::grid.arrange(pshape, pscale)

  # Plot joint log likelihood profile
  iteration_ll |> 
    ggplot(aes(x = shape, y = scale, z = loglik)) + 
    geom_contour_filled() +
    geom_contour(color = "black", bins = 10) +
    theme_minimal()
}


#### End of code. 
###########################
###########################