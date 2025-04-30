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
library(patter)
library(proj.verse)
files_source_r(here_src())

#### Load data
pars_model_move <- qs::qread(here_input("pars-model-move-best.qs"))
iteration       <- qs::qread(here_input_sim("iteration-patter-optim.qs"))


###########################
###########################
#### Set up analysis 

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


###########################
###########################
#### Grid-based search

# TO DO
# * Define new iteration data.table in prepare-analysis.R 
# * Options: marginal parameter estimates, 2D grid (expensive)
# * Add output/sim/grid/ folder
# * Run optimisation (a couple of times for each point)
# * (optional) Get SE from span around optimum after smoothing


#### End of code. 
###########################
###########################