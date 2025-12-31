###########################
###########################
#### rebuild.R

#### Aims
# 1) Rerun selected scripts & rebuild data outputs, after cleaning 

#### Prerequisites
# 1) Clean project 


###########################
###########################
#### Rerun selected scripts

library(proj.verse)

# List scripts
scripts <- list.files(here_r(), full.names = TRUE)

# Check scripts
scripts

# (1) Set up project 
# > This code can be run on Windows/MacOS/Linux
# > Run scripts 1:13 sequentially, in isolation
# > These scripts can be run independently locally/on a server to set up files etc.
cl_lapply(scripts[1:13L], function(script) {
  print(script)
  callr::rscript(script)
})

# (2) Simulate data
# > Data simulation is implemented from R via patter
# > This code should be run locally (on SIA-LAVENDED)
# > Then copy data/input/sim onto linux server 
# callr::rscript(script[14])

# (3) Prepare analysis
# > This script can be run locally and/or on the server

# (4) Trial patter
# > This script runs patter via R and should be run locally

# (5) Run patter
# > Run Julia scripts at scale via bash on server 
# > Monitor workflow via run-patter.R on server

# (6) Synthesise results
# > These scripts need to be run on the server, where patter outputs live


#### End of code. 
###########################
###########################