###########################
###########################
#### run-patter.R

#### Aims
# 1) This script provides a generic workflow for running particle algorithms 

#### Prerequisites
# 1) Generate particle algorithm outputs via run-algorithms.jl
# 2) This code should be run on the same machine (where the particle outputs live)

# TO DO REVISE THIS CODE
# SPLIT by chain_id & for each chain
# - loop over blocks
# - read in all pou files 
# - filter by chain time series ENSURING WE ONLY KEEP ONE FILE!
# - define output files


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

## For "sim" iteration[1, ]:
# * Memory required per iteration : 251.46 MB 
# * Time required per iteration   : 1.012 s
# * ETA for 210 iteration         : 3.5 mins on 1 cl (1.012 * 210 / 60)

## For "real" iteration[1, ]:
# * Memory required per iteration : TO DO
# * Time required per iteration   : TO DO
# * ETA for 2723 iteration        : TO DO
# > We can safely run ≈ TO DO N CPUs

## Subset iterations
# We will produce maps for all iterations with:
# a) julia = FALSE
# b) julia = TRUE with pou-{i}.feather files
#    (these were produced when both forward and backward filters were run successfully)
iteration[, julia_success := sapply(iteration$folder_output_block, function(folder) {
  length(list.files(folder, "pou-")) > 0L
})]
iteration[julia == FALSE, julia_success := as.logical(NA)]
table(iteration$julia_success[iteration$julia])
iteration <- iteration[!julia | (julia & julia_success), ]
stopifnot(nrow(iteration) > 0L)

## Define cluster
tic()
cl <- parallel::makeCluster(5L)
parallel::clusterEvalQ(cl, {
  library(data.table)
  library(dtplyr)
  library(dplyr, warn.conflicts = FALSE)
})
toc()

## Make maps
tic()
cl_lapply(split(iteration, seq_len(nrow(iteration))), 
          .cl = cl,
          .fun = function(it) {
  
  # Read data
  # tic()
  # it <- iteration[1, ]
  map         <- terra::rast("./data/input/map.tif")
  regions     <- terra::rast("./data/input/regions.tif")
  it_block    <- seq(it$time_start, it$time_end, by = "2 mins")
  it_timeline <- arrow::read_feather(it$file_timeline)
  it_timeline <- it_timeline[!is.na(block), ]
  
  if (it$julia) {
    
    # Compute a data.table of POU (x, y, mark = probability mass)
    coord <- 
      # List files 
      # (convergence failures handled above)
      list.files(it$folder_output_block, full.names = TRUE, pattern = "pou-") |> 
      lapply(arrow::read_feather) |>
      rbindlist() |> 
      arrange(timestep, x, y) |> 
      # Filter time series by relevant period
      # TO DO 
      filter(timestep %in% it_timeline$timestep)
      group_by(x, y) |> 
      summarise(mark = sum(mark)) |> 
      ungroup() |> 
      mutate(mark = mark / nrow(it_timeline)) |> 
      as.data.table()
    
    # Verify that weights sum to one 
    stopifnot(isTRUE(all.equal(1, sum(coord$mark))))
    
    # Map occupancy 
    occupancy <- terra::rasterize(coord, map, values = coord$mark)
    
  } else {
    
    # For time blocks with julia = FALSE, we use a uniform map
    occupancy <- terra::setValues(map, 1)
    occupancy <- terra::mask(occupancy, map)
    occupancy <- spatNormalise(occupancy)
    
  }
  
  occupancy <- terra::classify(occupancy, cbind(NA, 0))
  occupancy <- terra::mask(occupancy, map)
  names(occupancy) <- "map_value"
  # terra::plot(occupancy)
  
  # Estimate residency
  residency <- 
    occupancy |> 
    terra::zonal(regions, fun = "sum", na.rm = TRUE) |> 
    lazy_dt() |> 
    mutate(unit_id           = it$unit_id,
           individual_id     = it$individual_id,
           time_id           = it$time_id,
           sensitivity       = it$sensitivity,
           sensitivity_label = it$sensitivity_label) |>
    select("unit_id", "individual_id", "time_id", "sensitivity", "sensitivity_label",
           region = "map_value", estimate = "map_value.1") |>
    as.data.table()

  # Write to file
  qs::qsave(residency, it$file_residency)
  terra::writeRaster(occupancy, it$file_occupancy, overwrite = TRUE)
  
  # (optional) Cleanup pou-{i}.feather files
  # unlink(list.files(it$folder_output_block, full.names = TRUE, pattern = "pou-"))
  
  # lobstr::mem_used()
  # toc()
  invisible(NULL)
  
})
toc()

# terra::rast(iteration$file_occupancy[1])


#### End of code. 
###########################
###########################