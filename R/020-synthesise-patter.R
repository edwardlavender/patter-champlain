###########################
###########################
#### synthesise-patter.R

#### Aims
# 1) This script collates particle outputs into occupancy maps and residency estimates
# - Particles are produced in batches by block
# - We are interested in chain-level (seasonal) patterns (each chain comprises multiple blocks)
# - Here, we collate outputs for each chain into a map and residency estimates

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
library(spatial.extensions)
library(tictoc)
files_source_r(here_src())

#### Load data
map <- terra::rast("./data/input/map.tif")


###########################
###########################
#### Select analysis

#### Define analysis 
analysis <- "sim"
analysis <- "real"
subanalysis <- "main"

#### Define analysis-specific data
here_input_analysis <- switch_here_input_analysis_subanalysis(analysis, subanalysis)
iteration           <- qs::qread(here_input_analysis("iteration.qs"))


###########################
###########################
#### Synthesise outputs

#### Examine timings & folder sizes
#
## For "sim" iteration:
# * 12 min on 1 cl
#
## For "real" iteration[1, ]:
# * 15 min on 5 cl

#### Subset iterations
# We will produce maps for all chains where _all blocks_ have:
# a) julia = FALSE
# b) julia = TRUE with {n_batch} pou-{i}.feather files
#    (these were produced when both forward and backward filters were run successfully)
# i.e., We drop any chains in which one or more blocks failed
iteration[, chain_group := paste(individual_id, chain_id, sensitivity)]
iteration[, julia_success := sapply(iteration$folder_output_block, function(folder) {
  length(list.files(folder, "pou-")) == iteration$n_batch[1]
})]
iteration[julia == FALSE, julia_success := as.logical(NA)]
table(iteration$julia_success[iteration$julia])
iteration <-
  iteration |>
  group_by(chain_group) |> 
  filter(all(!julia | (julia & julia_success))) |> 
  as.data.table()
stopifnot(nrow(iteration) > 0L)

#### Pre-compute blank (starting) occupancy map
# For each chain, this is a starting point that is updated block-by-block
occupancy_zero <- terra::setValues(map, 0)
occupancy_zero <- terra::mask(occupancy_zero, map)
terra::writeRaster(occupancy_zero, 
                   here_input("occupancy-zero.tif"), 
                   overwrite = TRUE)
# terra::plot(occupancy_zero)

#### Pre-compute uniform map
# We assume uniform distribution over the study area when julia = FALSE
# (I.e., a long time after the last detection)
occupancy_uniform <- terra::setValues(map, 1)
occupancy_uniform <- terra::mask(occupancy_uniform, map)
occupancy_uniform <- spatNormalise(occupancy_uniform)
terra::writeRaster(occupancy_uniform, 
                   here_input("occupancy-uniform.tif"), 
                   overwrite = TRUE)
# terra::plot(occupancy_uniform)

#### Define cluster
tic()
cl <- parallel::makeCluster(5L)
parallel::clusterEvalQ(cl, {
  library(data.table)
  library(dtplyr)
  library(dplyr, warn.conflicts = FALSE)
})
toc()

#### Synthesise outputs
# For each chain, compute file_occupancy and file_residency
toc()
cl_lapply(
  split(iteration, iteration$chain_group), 
  .cl = cl,
  .fun = function(chain) {
  
  #### Read datasets
  # chain           <- split(iteration, iteration$chain_group)[[1]]
  map               <- terra::rast("./data/input/map.tif")
  regions           <- terra::rast("./data/input/regions.tif")
  occupancy         <- terra::rast("./data/input/occupancy-zero.tif")
  occupancy_uniform <- terra::rast("./data/input/occupancy-uniform.tif")
  
  #### Update map for each block in the chain 
  # This is implemented iteratively to keep memory use down
  for (i in seq_len(nrow(chain))) {
    
    # Define iteration 
    it <- chain[i, ]
    
    # Define block timeline 
    # * block_timeline is timeline of interest 
    # * Each chain is made up of several non-overlapping block_timeline(s) 
    block_timeline <- seq(it$block_start, it$block_end, by = "2 mins")
    
    # Update occupancy map for block 
    if (it$julia) {
      
      # Define modelled timeline
      # * it_timeline is the modelled timeline for a block (may be longer than the block)
      it_timeline <- arrow::read_feather(it$file_timeline)
      
      # Compute the total weight per cell (cell, mark) over all time steps
      marks <- 
        # Load all POU files
        list.files(it$folder_output_block, full.names = TRUE, pattern = "pou-") |> 
        lapply(arrow::read_feather) |>
        rbindlist() |> 
        lazy_dt(immutable = FALSE) |>
        # Assign timestamps and filter by block_timeline
        mutate(timestamp = it_timeline$timestamp[match(timestep, it_timeline$timestep)]) |> 
        filter(timestamp %in% block_timeline) |> 
        # For each cell, compute the overall weight 
        mutate(cell = terra::cellFromXY(map, cbind(x, y))) |> 
        group_by(cell) |> 
        summarise(mark = sum(mark)) |> 
        as.data.table()
      
      # Update map
      occupancy[marks$cell] <- occupancy[marks$cell] + marks$mark
        
    } else {
      
      # For un-modelled blocks, assume uniform weights & update map for all time steps
      occupancy <- occupancy + (occupancy_uniform * length(block_timeline))
    
    }
    
  }
  
  #### Process map
  occupancy <- terra::classify(occupancy, cbind(NA, 0))
  occupancy <- terra::mask(occupancy, map)
  occupancy <- spatial.extensions::spatNormalise(occupancy)
  names(occupancy) <- "map_value"
  # terra::plot(occupancy)
  
  #### Estimate residency by region
  residency <- 
    occupancy |> 
    terra::zonal(regions, fun = "sum", na.rm = TRUE) |> 
    lazy_dt() |> 
    mutate(unit_id           = it$unit_id,
           individual_id     = it$individual_id,
           chain_id          = chain$chain_id[1],
           sensitivity       = it$sensitivity,
           sensitivity_label = it$sensitivity_label) |>
    select("unit_id", "individual_id", "chain_id", "sensitivity", "sensitivity_label",
           region = "map_value", estimate = "map_value.1") |>
    as.data.table()
  
  #### Write to file
  terra::writeRaster(occupancy, chain$file_occupancy[1], overwrite = TRUE)
  qs::qsave(residency, chain$file_residency[1])

  #### (optional) Cleanup pou-{i}.feather files
  # unlink(list.files(it$folder_output_block, full.names = TRUE, pattern = "pou-"))
  
  invisible(NULL)
  
})
toc()

# terra::rast(iteration$file_occupancy[1])


#### End of code. 
###########################
###########################