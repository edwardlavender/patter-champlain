###########################
###########################
#### 03-clean.R

#### Aims
# 1) Clean the project 
# * This script unlinks the data/ folder and rebuilds the essential directories
# * I.e., directories that are, generally, not build by the scripts
# * This code is used during project development to maintain the directory structure
#   and ensure all processed data/ files are up to date.
# > Run this code
# > Rerun all scripts, restarting R before each script

#### Prerequisites
# 1) NA


###########################
###########################
#### Implement cleanup

#### Set up
rm(list = ls())
library(proj.verse)
files_source_r()

#### Unlink data/ folder
# Only unlink the data/ folder for complete project rebuilds
if (FALSE) {
  unlink(here_data(), recursive = TRUE)
}

#### Rebuild directories 

dirs.create("data")

dirs.create(here_input())
dirs.create(here_input_real())
dirs.create(here_input_sim())
dirs.create(here_input("vmap"))

dirs.create(here_debug())

dirs.create(here_data("inst"))

dirs.create(here_output_sim())

dirs.create(here_output_sim_main())
dirs.create(here_output_sim_main("logs", "R"))
dirs.create(here_output_sim_main("logs", "R-CMD-BATCH"))
dirs.create(here_output_sim_main("runs"))
dirs.create(here_output_sim_main("synthesis"))

dirs.create(here_output_sim_optim())
dirs.create(here_output_sim_optim("logs", "R"))
dirs.create(here_output_sim_optim("logs", "R-CMD-BATCH"))
dirs.create(here_output_sim_optim("runs"))
dirs.create(here_output_sim_optim("synthesis"))

dirs.create(here_output_sim_grid())
dirs.create(here_output_sim_grid("logs", "R"))
dirs.create(here_output_sim_grid("logs", "R-CMD-BATCH"))
dirs.create(here_output_sim_grid("runs"))
dirs.create(here_output_sim_grid("synthesis"))

dirs.create(here_output_real_main())
dirs.create(here_output_real_main("logs", "R"))
dirs.create(here_output_real_main("logs", "R-CMD-BATCH"))
dirs.create(here_output_real_main("runs"))
dirs.create(here_output_real_main("synthesis"))

dirs.create(here_data("supp"))
dirs.create(here_data("supp", "model-move"))
dirs.create(here_data("supp", "model-obs"))

#### Record tree
# See 02-clone.R


#### End of code. 
###########################
###########################
