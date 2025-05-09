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

# Run selected scripts, sequentially, in isolation (~x mins)
cl_lapply(scripts[1:14], function(script) {
  print(script)
  callr::rscript(script)
})


#### End of code. 
###########################
###########################