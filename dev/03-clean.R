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

# data/
dirs.create("data")

# data/input/
dirs.create(here_input())
dirs.create(here_input_real())
dirs.create(here_input_sim())
dirs.create(here_input("vmap"))

# data/debug/ 
dirs.create(here_debug())

# data/inst/
dirs.create(here_data("inst"))

# data/output/sim/main/
dirs.create(here_output_sim())
dirs.create(here_output_sim_main())
dirs.create(here_output_sim_main("logs"))
dirs.create(here_output_sim_main("runs"))
dirs.create(here_output_sim_main("synthesis"))

# data/output/real/main/
dirs.create(here_output_real_main())
dirs.create(here_output_real_main("logs"))
dirs.create(here_output_real_main("runs"))
dirs.create(here_output_real_main("synthesis"))

# data/output/sim/main/
dirs.create(here_data("supp"))
dirs.create(here_data("supp", "model-move"))
dirs.create(here_data("supp", "model-obs"))

# data/fig/
dirs.create(here_fig())
dirs.create(here_fig("local"))
dirs.create(here_fig("local", "qgis"))
dirs.create(here_fig("model", "model-move"))
dirs.create(here_fig("model", "model-obs"))
dirs.create(here_fig("tables"))

#### Record tree
# See 02-clone.R


#### End of code. 
###########################
###########################
