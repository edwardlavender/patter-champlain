###########################
###########################
#### sim-analysis-patter-01.R

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
#### Run analysis 

#### (optional) Set up cluster
# TO DO Set up cluster
map <- here_input("map.tif")
set_map(map)

#### Run filter & optimisation
cl_lapply_workflow(.iteration = iteration[1, ], 
                   .datasets = NULL, 
                   .constructor = constructor_pf_filter_loglik_optim, 
                   .algorithm = pf_filter_loglik_optim)


#### End of code. 
###########################
###########################