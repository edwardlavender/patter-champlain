###########################
###########################
#### run-patter.R

#### Aims
# 1) This script provides a generic workflow for running particle algorithms 

#### Prerequisites
# 1) Generate particle algorithm outputs via run-algorithms.jl
# 2) This code should be run on the same machine (where the particle outputs live)


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Set global options
Sys.setenv("JULIA_SESSION" = FALSE)
set.seed(123L)

#### Load essential packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(proj.verse)
library(tictoc)
files_source_r(here_src())

#### Load data
# NA


###########################
###########################
#### Select analysis

#### Define analysis 
analysis <- "sim"
# analysis <- "real"
subanalysis <- "main"

#### Define analysis-specific data
here_input_analysis <- switch_here_input_analysis_subanalysis(analysis, subanalysis)
iteration           <- qs::qread(here_input_analysis("iteration.qs"))


###########################
###########################
#### Run patter

#### (1) Run Julia workflow
# Run 001-run-algorithms.jl via tmux/bash

#### (2) Run R workflows (mapping)

# For "sim" iteration[1, ]:
# * Memory required per iteration : 251.46 MB 
# * Time required per iteration   : 1.012 s
# * ETA for 210 iteration         : 3.5 mins on 1 cl (1.012 * 210 / 60)

# For "real" iteration[1, ]:
# * Memory required per iteration : TO DO
# * Time required per iteration   : TO DO
# * ETA for 2723 iteration        : TO DO
# > We can safely run ≈ TO DO N CPUs

# Define cluster
tic()
cl <- parallel::makeCluster(2L)
parallel::clusterEvalQ(cl, {
  library(data.table)
  library(dtplyr)
  library(dplyr, warn.conflicts = FALSE)
})
toc()

# Make maps
tic()
cl_lapply(split(iteration, seq_len(nrow(iteration)))[1], 
          .cl = cl,
          .fun = function(it) {
  
  # Read data
  # tic()
  # it <- iteration[1, ]
  map         <- terra::rast("./data/input/map.tif")
  it_timeline <- arrow::read_feather(it$file_timeline)$timestamp
  
  # Compute a data.table of POU (x, y, mark = probability mass)
  coord <- 
    list.files(it$folder_output, full.names = TRUE, pattern = "pou-") |> 
    lapply(arrow::read_feather) |>
    rbindlist() |> 
    arrange(timestep, x, y) |> 
    group_by(x, y) |> 
    summarise(mark = sum(mark)) |> 
    ungroup() |> 
    mutate(mark = mark / length(it_timeline)) |> 
    as.data.table()
  
  # Verify that weights sum to one 
  stopifnot(all.equal(1, sum(coord$mark)))
  
  # Map occupancy 
  occupancy <- terra::rasterize(coord, map, values = coord$mark)
  occupancy <- terra::classify(occupancy, cbind(NA, 0))
  occupancy <- terra::mask(occupancy, map)
  names(occupancy) <- "map_value"
  # terra::plot(occupancy)

  # Write to file
  terra::writeRaster(occupancy, it$file_occupancy, overwrite = TRUE)
  
  # (optional) Cleanup pou-{i}.feather files
  # unlink(list.files(it$folder_output, full.names = TRUE, pattern = "pou-"))
  
  # lobstr::mem_used()
  # toc()
  invisible(NULL)
  
})
toc()

# terra::rast(iteration$file_occupancy[1])


#### End of code. 
###########################
###########################