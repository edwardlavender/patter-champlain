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
library(prettyGraphics)
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
  
  data.table(sensitivity    = "best", 
             mobility       = pm$mobility[1], 
             shape          = pm$shape[1], 
             scale          = pm$scale[1], 
             phi            = pm$phi[1], 
             receiver_alpha = po$receiver_alpha[1], 
             receiver_beta  = po$receiver_beta[1], 
             receiver_gamma = po$receiver_gamma[1]),

  data.table(sensitivity    = "step(-)", 
             mobility       = pm$mobility[2], 
             shape          = pm$shape[2], 
             scale          = pm$scale[2], 
             phi            = pm$phi[1], 
             receiver_alpha = po$receiver_alpha[1], 
             receiver_beta  = po$receiver_beta[1], 
             receiver_gamma = po$receiver_gamma[1]),
  
  data.table(sensitivity    = "step(+)", 
             mobility       = pm$mobility[3], 
             shape          = pm$shape[3], 
             scale          = pm$scale[3], 
             phi            = pm$phi[1], 
             receiver_alpha = po$receiver_alpha[1], 
             receiver_beta  = po$receiver_beta[1], 
             receiver_gamma = po$receiver_gamma[1]),
  
  data.table(sensitivity    = "angle(-)", 
             mobility       = pm$mobility[1], 
             shape          = pm$shape[1], 
             scale          = pm$scale[1], 
             phi            = pm$phi[2], 
             receiver_alpha = po$receiver_alpha[1], 
             receiver_beta  = po$receiver_beta[1], 
             receiver_gamma = po$receiver_gamma[1]),
  
  data.table(sensitivity    = "angle(+)", 
             mobility       = pm$mobility[1], 
             shape          = pm$shape[1], 
             scale          = pm$scale[1], 
             phi            = pm$phi[3], 
             receiver_alpha = po$receiver_alpha[1], 
             receiver_beta  = po$receiver_beta[1], 
             receiver_gamma = po$receiver_gamma[1]),
  
  data.table(sensitivity    = "ac(-)", 
             mobility       = pm$mobility[1], 
             shape          = pm$shape[1], 
             scale          = pm$scale[1], 
             phi            = pm$phi[1], 
             receiver_alpha = po$receiver_alpha[2], 
             receiver_beta  = po$receiver_beta[2], 
             receiver_gamma = po$receiver_gamma[2]),
  
  data.table(sensitivity    = "ac(+)", 
             mobility       = pm$mobility[1], 
             shape          = pm$shape[1], 
             scale          = pm$scale[1], 
             phi            = pm$phi[1], 
             receiver_alpha = po$receiver_alpha[3], 
             receiver_beta  = po$receiver_beta[3], 
             receiver_gamma = po$receiver_gamma[3])
  
) |> 
  mutate(parameter_id = row_number(), .before = "sensitivity") |> 
  as.data.table()

qs::qsave(pars, here_input("pars-patter.qs"))

pars |> 
  mutate(parameter_id = as.character(parameter_id), 
         mobility = as.character(mobility)) |> 
  tidy_numbers(digits = c(4, 4, 4, 4, 4, 0)) |> 
  tidy_write(here_fig("tables", "pars.txt"))


#### End of code.
###########################
###########################