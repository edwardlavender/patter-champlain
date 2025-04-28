###########################
###########################
#### develop-model.R

#### Aims
# 1) Collect state-space model parameters 

#### Prerequisites
# 1) Define model parameters


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
#### Collate pars

pm <- pars_model_move
po <- pars_model_obs

pars <- rbind(
  
  data.table(sensitivity = "best", 
             mobility = pm$mobility[1], 
             shape = pm$shape[1], 
             scale = pm$scale[1], 
             phi = pm$phi[1], 
             pars_model_obs),
  
  data.table(sensitivity = "step(-)", 
             mobility = pm$mobility[2], 
             shape = pm$shape[2], 
             scale = pm$scale[2], 
             phi = pm$phi[1], 
             pars_model_obs),
  
  data.table(sensitivity = "step(+)", 
             mobility = pm$mobility[3], 
             shape = pm$shape[3], 
             scale = pm$scale[3], 
             phi = pm$phi[1], 
             pars_model_obs),
  
  data.table(sensitivity = "angle(-)", 
             mobility = pm$mobility[1], 
             shape = pm$shape[1], 
             scale = pm$scale[1], 
             phi = pm$phi[2], 
             pars_model_obs),
  
  data.table(sensitivity = "angle(+)", 
             mobility = pm$mobility[1], 
             shape = pm$shape[1], 
             scale = pm$scale[1], 
             phi = pm$phi[3], 
             pars_model_obs)
) |> 
  mutate(parameter_id = row_number(), .before = "sensitivity") |> 
  as.data.table()

qs::qsave(pars, here_input("pars-patter.qs"))

#### End of code.
###########################
###########################
