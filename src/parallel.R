cl_init <- function(iteration, cl, varlist) {
  parallel::clusterExport(cl = cl, varlist = varlist)
  parallel::clusterEvalQ(cl = cl, {
    
    # Load packages 
    library(data.table)
    library(dtplyr)
    library(dplyr, warn.conflicts = FALSE)
    library(JuliaCall)
    library(patter)
    library(patter.workflows)
    library(proj.verse)
    library(tictoc)
    files_source_r(here_src())
    
    # Initialise Julia 
    particle_startup(.sim = NULL, .cl = ncl)
    expect_no_geospatial()
    
    # Set maps
    set_map(here_input("map.tif"))
    set_vmap(.vmap = here_input("vmap", iteration$mobility[1], "vmap.tif"))
    invisible(NULL)
  })
  
}