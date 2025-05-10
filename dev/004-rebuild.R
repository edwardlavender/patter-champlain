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

# Rebuild on Windows/MacOS:
# > Run scripts 1:14 sequentially, in isolation (~6.5 mins)
# > Run subsequent scripts manually, e.g., for sim and real analyses
cl_lapply(scripts[1:14], function(script) {
  print(script)
  callr::rscript(script)
})

# Rebuild on Linux: 
# > Find in files julia_connect()
# > Reun scripts 1:13 (JULIA_SESSION = FALSE) (~5 mins)
# > Restart RStudio Project
# > Run other scripts manually 
cl_lapply(scripts[1:13], function(script) {
  print(script)
  callr::rscript(script)
})


#### End of code. 
###########################
###########################