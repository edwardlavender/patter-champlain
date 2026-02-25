###########################
###########################
#### prepare-analysis.R

#### Aims
# 1) Prepares detection datasets, iteration data.tables and directories for analyses
#    (This includes analyses of both simulation and real-world datasets)

#### Prerequisites
# 1) Previous scripts
# 2) Following ?patter.workflows::`patter.workflows-package`


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
library(ggplot2)
library(lubridate)
library(patter)
library(proj.verse)
library(spatial.extensions)
library(tictoc)
library(truncdist)
files_source_r(here_src())

#### Load data
pars       <- qs::qread(here_input("pars-patter.qs"))
fish       <- qs::qread(here_input("fish.qs"))


###########################
###########################
#### Select analysis

#### Define analysis 
analysis <- "sim"
# analysis <- "real"
subanalysis <- "main"

#### Define analysis-specific data
here_input_analysis <- switch_here_input_analysis_subanalysis(analysis, subanalysis)
detections          <- qs::qread(here_input_analysis("detections.qs"))
moorings            <- qs::qread(here_input_analysis("moorings.qs"))
if (analysis == "sim") {
  paths <- qs::qread(here_input_sim("main", "paths.qs"))
}


###########################
###########################
#### Process detections

#### Process timestamps (~5 s)
# Round observations to the nearest two minutes & drop duplicate observations
detections <- 
  detections |> 
  lazy_dt() |> 
  mutate(timestamp = round_date(timestamp, "2 mins")) |> 
  group_by(individual_id, timestamp, receiver_id) |>
  slice(1L) |>
  ungroup() |> 
  arrange(individual_id, timestamp, receiver_id) |> 
  as.data.table()

#### Filter individuals
if (analysis == "real") {
  # We focus on the study period of 2014 -> 2017 
  # We focus on individuals detected in the study period 
  # (For those individuals, we retain data for 2013
  #  as it may support initialisation of filter)
  # A) Identify individuals to drop: 
  sort(unique(fish$date))
  individuals_to_exclude <- 
    detections |> 
    group_by(individual_id) |> 
    summarise(drop = all(timestamp < as.POSIXct("2014-12-01 00:00:00"))) |>
    filter(drop) |> 
    pull(individual_id)
  # B) Filter detections accordingly 
  detections <- detections[!(individual_id %in% individuals_to_exclude), ]
}

nrow(detections)


###########################
###########################
#### Define unitsets

#### Overview 
# unitsets defines the individual/time blocks for analysis (defined by unit_id)
# For simulations, we model each individual and the one-month simulated time series
# For real-world analysis, we model each individual and every seasonal block from
# the start to the end of the study period

if (analysis == "sim") {
  
  #### Define individuals
  # We model 1:100 individuals
  individuals <- sort(unique(paths$path_id))
  
  #### Define study timeframe
  # We simulated data over one month
  # interval(min(paths$timestamp), max(paths$timestamp))
  # 2025-01-01 UTC--2025-01-31 23:58:00 UTC
  timeframe <-
    tribble(
      ~chain_id,  ~chain_interval,
      "2014-Jan", interval(min(paths$timestamp), max(paths$timestamp), tzone = "UTC")
    ) |> 
    mutate(chain_start = int_start(chain_interval), 
           chain_end = int_end(chain_interval)) |> 
    select(-chain_interval) |> 
    as.data.table()
  
  #### Define time blocks
  # All data were simulated over the same one-month period
  blocks <- 
    min(paths$timestamp) |>
    format("%y-%b") |>
    tolower()

} else if (analysis == "real") {

  #### Define individuals
  # We focus on individuals detected in at least two seasons 
  # (as defined in the detections data.table)
  individuals <- sort(unique(detections$individual_id))
  
  #### Define study timeframe
  # We are interested in seasonal patterns, defined as: 
  # * Winter: December 1–March 31
  # * Spring: April 1–May 31
  # * Summer: June 1–September 30
  # * Fall: October 1–November 30
  # For the start of the study period, we use Winter 2014
  # * This is the first full season of data for most fish (tagged during fall 2013)
  sort(unique(fish$date[lubridate::year(fish$date) > 2013]))
  # For the end of the study period, we use Spring 2015
  # * When was the array dismantled? 
  max(detections$timestamp)
  # Define study period
  timeframe <-
    tribble(
      ~chain_id,     ~chain_interval,
      "2014-winter", interval("2014-12-01 00:00:00", "2015-03-31 23:58:00", tzone = "UTC"),
      "2015-spring", interval("2015-04-01 00:00:00", "2015-05-31 23:58:00", tzone = "UTC"),
      "2015-summer", interval("2015-06-01 00:00:00", "2015-09-30 23:58:00", tzone = "UTC"),
      "2015-fall",   interval("2015-10-01 00:00:00", "2015-11-30 23:58:00", tzone = "UTC"),
      "2016-winter", interval("2015-12-01 00:00:00", "2016-03-31 23:58:00", tzone = "UTC"),
      "2016-spring", interval("2016-04-01 00:00:00", "2016-05-31 23:58:00", tzone = "UTC"),
      "2016-summer", interval("2016-06-01 00:00:00", "2016-09-30 23:58:00", tzone = "UTC"),
      "2016-fall",   interval("2016-10-01 00:00:00", "2016-11-30 23:58:00", tzone = "UTC"),
      "2017-winter", interval("2016-12-01 00:00:00", "2017-03-31 23:58:00", tzone = "UTC"),
      "2017-spring", interval("2017-04-01 00:00:00", "2017-05-31 23:58:00", tzone = "UTC")
    ) |> 
    mutate(chain_start = int_start(chain_interval), 
           chain_end = int_end(chain_interval)) |> 
    select(-chain_interval) |> 
    as.data.table()
  
  #### Define time blocks
  # We use monthly blocks over the study period
  blocks <- 
    seq(as.POSIXct("2014-12-01 00:00:00"), 
                 as.POSIXct("2017-05-31 23:58:00"), 
                 by = "months") |> 
    format("%y-%b") |>
    tolower()
  
}

#### Define unitsets
# For the real-world analysis, there are 2070 individual-month blocks 
unitsets <- 
  CJ(individual_id = individuals, block_id = blocks) |>
  mutate(block_start := parse_date_time(block_id, "y-b", tz = "UTC"),
         block_end := block_start + months(1) - 60 * 2,
         julia = TRUE) |> 
  arrange(individual_id, block_start) |>
  mutate(unit_id = 1:n()) |>
  select("unit_id", "individual_id", 
         "block_id", "block_start", "block_end", 
         "julia") |>
  as.data.table()

#### Assign chains
# We implement the algorithm in blocks (e.g., individuals/months):
# * For each block, Julia produces callstats, diagnostics and pou-{i} files
# * The pou-{i} file contains :timestep, :id, :x, :y, :mark
# We aggregate patterns over multiple blocks (chains)
# * For the simulation analysis, blocks = chains
# * For the real-world analysis, chains are individuals/seasons
# Hence, we assign chain_id to the iteration data.table
# * After analysis, we aggregate maps/residency over all time stamps in each chain
unitsets <- 
  unitsets |>
  left_join(
    timeframe,
    join_by(between(block_start, chain_start, chain_end))) |> 
  select("unit_id", "individual_id", 
         "chain_id", "chain_start", "chain_end", 
         "block_id", "block_start", "block_end",
         "julia") |> 
  as.data.table()

#### For the real-analysis, identify the subset of blocks that require modelling
if (analysis == "real") {
  
  # Identify the last detection for each individual
  detections_end <- 
    detections |> 
    group_by(individual_id) |> 
    summarise(timestamp_last_detection = max(timestamp)) |> 
    as.data.table()
  
  # Identify which blocks require modelling
  # * We model all blocks for _at least_ month month after the last detection
  # * (I.e., for the rest of the month and the next whole month)
  unitsets <- 
    unitsets |> 
    left_join(detections_end, by = "individual_id") |> 
    mutate(timestamp_one_block_after_last_detection = timestamp_last_detection + months(1), 
           julia = block_start <= timestamp_one_block_after_last_detection) |> 
    as.data.table()

  # Visual checks
  unitsets[individual_id == individual_id[1], ]
  
  # Cleanup
  unitsets <- 
    unitsets |> 
    select(-c(timestamp_last_detection, timestamp_one_block_after_last_detection)) |> 
    as.data.table()
  
}

# For the real-world analysis, there are 2070 blocks, of which 1510 require modelling (73 %)
nrow(unitsets)
table(unitsets$julia)
# View(unitsets)


###########################
###########################
#### Prepare iterations

#### Review iteration parameters
## (1) Guess the number of batches
# We can compute the number of batches for:
#   n_particles * n timesteps * 4 states * 100 CPUs if 50e3 MB memory available
# We have used one-month blocks (22320 time steps)
#   The actual number of time series varies because process time series
#   to start/end with detections where possible. 
#   So these choices are reviewed below. 
p_batch(p_mem(1500, 22320, 4, 100), 50e3)
p_batch(p_mem(2000, 22320, 4, 100), 50e3)

#### Define iteration 
iteration <- 
  unitsets |> 
  cross_join(pars) |> 
  mutate(index = row_number(), .before = 1L) |> 
  mutate(
    
    # Define input files 
    # * Some files depend on both unit_id & sensitivity parameters
    # * For convenience, we store all files in an {individual_id}/{unit_id}/{parameter_id} directory 
    folder_input          = file.path("data", "input", analysis, subanalysis, "blocks", 
                                      individual_id, block_id, parameter_id), 
    file_timeline         = file.path(folder_input, "timeline-algorithm.feather"),
    file_acoustics        = file.path(folder_input, "acoustics.feather"),
    file_containers_fwd   = file.path(folder_input, "containers-fwd.feather"),
    file_containers_bwd   = file.path(folder_input, "containers-bwd.feather"),
    file_t_resample_fwd   = file.path(folder_input, "t-resample-fwd.feather"),
    file_t_resample_bwd   = file.path(folder_input, "t-resample-bwd.feather"),
    file_vmap             = file.path("data", "input", "vmap", as.integer(mobility), "vmap.tif"),
    
    # Define block output files (unit-specific & sensitivity specific)
    # * We use .feather to record outputs
    # * We can write these from Julia & read them into R correctly
    folder_output_block   = file.path("data", "output", analysis, subanalysis, "blocks", 
                                      individual_id, block_id, parameter_id),
    file_callstats_filter = file.path(folder_output_block, "callstats-forward-filter.feather"),
    file_callstats        = file.path(folder_output_block, "callstats.feather"),
    file_diagnostics      = file.path(folder_output_block, "diagnostics.feather"),
    file_states           = file.path(folder_output_block, "states.feather"),
    
    # Define chain output files 
    folder_output_chain   = file.path("data", "output", analysis, subanalysis, "chains", 
                                      individual_id, chain_id, parameter_id),
    file_occupancy        = file.path(folder_output_chain, "occupancy.tif"),
    file_residency        = file.path(folder_output_chain, "residency.qs"),
    file_path_sim         = if_else(rep(analysis == "sim", n()),
                                    file.path(folder_output_chain, "path-sim.qs"),
                                    NA_character_), 
    file_occupancy_sim    = if_else(rep(analysis == "sim", n()),
                                    file.path(folder_output_chain, "occupancy-sim.tif"),
                                    NA_character_), 
    file_residency_sim    = if_else(rep(analysis == "sim", n()),
                                    file.path(folder_output_chain, "residency-sim.qs"),
                                    NA_character_), 
    
    # Add modelling columns
    n_move              = 1000L,
    n_particle_filter   = ifelse(analysis == "sim", 10000L, 20000L), 
    n_particle_smoother = ifelse(analysis == "sim", 1500L, 2000L),
    n_resample          = as.numeric(1000.0),
    n_batch             = ifelse(analysis == "sim", 9L, 25L)
  ) |> 
  as.data.table()

#### Record the number of detections
if (analysis == "sim") {
  # Add iteration$n_detections column 
  # (Note that for simulations we only need to match by individual_id, not block_start)
  iteration <- 
    iteration |> 
    left_join(
      detections |> 
        group_by(individual_id) |> 
        summarise(n_detections = n()) |> 
        as.data.table(), 
      by = "individual_id") |> 
    mutate(n_detections = if_else(is.na(n_detections), 0, n_detections), 
           n_detections = as.integer(n_detections)) |> 
    as.data.table()
  # For simulations, 96/100 runs produced detections
  iteration |> 
    filter(sensitivity == "best") |> 
    summarise(length(which(n_detections > 0L)))
  # The number of detections was ~1000
  iteration |> 
    filter(sensitivity == "best") |> 
    summarise(utils.add::basic_stats(n_detections))
}

#### Filter iterations
if (analysis == "real") {
  # To keep computations manageable, we do not run a real-world sensitivity analysis
  iteration <- iteration[sensitivity == "best", ]
}

#### Check iterations 
# Check nrow is feasible! 
nrow(iteration)
# Check vmap files exist
stopifnot(all(file.exists(iteration$file_vmap)))

#### Build directories 
if (FALSE) {
  unlink(iteration$folder_input, recursive = TRUE)
  unlink(iteration$folder_output_block, recursive = TRUE)
  unlink(iteration$folder_output_chain, recursive = TRUE)
}
dirs.create(iteration$folder_input)
dirs.create(iteration$folder_output_block)
dirs.create(iteration$folder_output_chain)


###########################
###########################
#### Create Julia inputs

#### Duration
# "sim": 2 m 05 s s on SIA-LAVENDED or 3 m 16 s on siam-linux20 (10 cl, 2 chunks per core) 
# "real": 5 m 38 s or 16 m 19 s on siam-linux20 (10 cl, 2 chunks per core)
# (There is some speed benefit of chunking)

#### Write options (derived for analysis = "sim")
# write.csv
# * Simple, avoids issues in Julia e.g., with time stamps, but:
# * ~2 min to write all files for simulations (below)
# * 71.3612 MB per iteration
# * 15,230.96 MB for all simulations
# fwrite
# * faster but causes issues with timestamps
# arrow::write_feather()
# * works with Julia, if we set object types in Julia
# * 46 s for simulations 
# * 0.69476 MB per iteration
# * 184.3025 MB for simulations
# write_feather_compressed() and compression_level = 9
# --> 46 s
# --> 0.270096 MB
# --> 81.95327 MB
# write_feather_compressed() and compression_level = 22
# * 2 min 32 s
# * 0.268848 MB per iteration
# * 81.21565 MB for all iterations 
# > We use write_feather_compressed() and compression_level = 9
# > The marginal gains of max compression are v. limited
#   compared to the speed cost of writing files (important for real-world)

#### Define iterations for Julia
iteration_julia <- iteration[julia == TRUE, ]
nrow(iteration_julia)

#### Write files 
overwrite <- FALSE
if (!all(file.exists(iteration_julia$file_timeline)) | overwrite) {
  
  pbo <- pbapply::pboptions(nout = 2L)
  cl_lapply(
    split(iteration_julia, seq_len(nrow(iteration_julia))), 
    .cl = 10L,
    .chunk = TRUE,
    .fun = function(d) {
      
      #### (optional) selected rows 
      # d <- iteration_julia[1, ]
      # d <- iteration_julia[index == 8, ]
      # if (file.exists(d$file_residency_sim)) {
      #   return(NULL)
      # }
      print(d$index)
      
      #### Define detections
      # We identify the detection time series for the block
      # To ensure blocks are independent, we try to start/end blocks with detection:
      # - This prevents loss of information when blocks are treated independently 
      # - Hence, we do not filter the real detection time series 
      #   (i.e., by excluding pre-2014 data which could inform initialisation locations)
      # - (However, note this is not always possible e.g., at the start of the time series)
      # To achieve blocks, we select all detections in the range: 
      #  (block_id, block_id - one block, block_id + one block)
      #  (These criteria stop us jumping too far backward/forward in time)
      
      dets   <- detections[individual_id == d$individual_id, ]
      if (nrow(dets) > 0L) {
        
        # (A) Extract detections during the relevant time block (may be zero)
        during <- before <- after <- NULL
        during <- dets[timestamp >= d$block_start & timestamp <= d$block_end, ]
        
        # (B) Ensure the detection time series covers the full block:
        # (i) Add detection(s) immediately before the block start, if needed
        if (nrow(during) == 0L || min(during$timestamp) > d$block_start) {
          before <- dets[timestamp >= d$block_start - months(1) & timestamp < d$block_start, ]
          if (nrow(before) > 0L) {
            before <- before[timestamp == max(timestamp), ]
          }
        }
        # (ii) Add the detection(s) immediately after the block end, if needed
        if (nrow(during) == 0L || max(during$timestamp) < d$block_end) {
          after <- dets[timestamp > d$block_end & timestamp <= d$block_end + months(1), ]
          if (nrow(after) > 0L) {
            after <- after[timestamp == min(timestamp), ]
          }
        }
        dets <- rbind(during, before, after)
      }
      dets <- 
        dets |> 
        arrange(timestamp, receiver_id) |> 
        as.data.table()
      
      #### Define file_timeline(s)
      # Define timeline for algorithm run
      # (This may be longer than the timeline for the block)
      timeline <- seq(min(c(d$block_start, dets$timestamp)), 
                      max(c(d$block_end, dets$timestamp)), 
                      by = "2 mins")
      timeline <- data.table(timestep = 1:length(timeline), 
                             timestamp = timeline)
      write_feather_compressed(timeline, d$file_timeline)
      
      #### Define moorings, with detection probability parameters
      moors <- 
        moorings |> 
        lazy_dt(immutable = TRUE) |> 
        mutate(receiver_alpha = d$receiver_alpha, 
               receiver_beta  = d$receiver_beta, 
               receiver_gamma = d$receiver_gamma) |> 
        filter(
          int_overlaps(
            interval(receiver_start, receiver_end),
            interval(min(timeline$timestamp), max(timeline$timestamp))
          )
        ) |>
        as.data.table()
      
      #### Define file_acoustics
      # This uses the timeline for the algorithm run
      accs <- assemble_acoustics(.timeline   = timeline$timestamp, 
                                 .detections = dets, 
                                 .moorings   = moors)
      write_feather_compressed(accs, d$file_acoustics)
      
      #### Define acoustic containers (file_containers_fwd, file_containers_bwd)
      # To optimise the use of acoustic containers:
      # - We use a threshold that is slightly below the default
      # - See setup-data-map.R
      containers_fwd <- containers_bwd <- NULL
      if (nrow(dets) > 1L) {
        threshold  <- 44159.98
        containers <- assemble_acoustics_containers(.timeline  = timeline$timestamp, 
                                                    .acoustics = accs,
                                                    .mobility  = d$mobility, 
                                                    .map       = NULL, 
                                                    .threshold = threshold)
        containers_fwd <- containers$forward
        containers_bwd <- containers$backward
        stopifnot(nrow(containers_fwd) > 0L & nrow(containers_bwd) > 0L)
        stopifnot(max(c(containers_fwd$radius, containers_bwd$radius)) <= threshold)
        write_feather_compressed(containers_fwd, d$file_containers_fwd)
        write_feather_compressed(containers_bwd, d$file_containers_bwd)
      }
      
      #### Define t_resample
      if (isTRUE(nrow(containers_fwd) > 0L)) {
        t_resample_fwd <- 
          data.table(timestep = sort(unique(which(timeline$timestamp %in% 
                                                    containers_fwd$timestamp))))
        stopifnot(nrow(t_resample_fwd) > 0L)
        write_feather_compressed(t_resample_fwd, d$file_t_resample_fwd)
      }
      if (isTRUE(nrow(containers_bwd) > 0L)) {
        t_resample_bwd <- 
          data.table(timestep = sort(unique(which(timeline$timestamp %in% 
                                                    containers_bwd$timestamp))))
        stopifnot(nrow(t_resample_bwd) > 0L)
        write_feather_compressed(t_resample_bwd, d$file_t_resample_bwd)
      }
    
      nothing()
    })
  pbapply::pboptions(pbo)
  
}

#### Review the number of time steps/batches
# For the real-world time series, each block contains up to ~22320 time steps
# But modelled time blocks may be longer because we start with detections before/after blocks
# The worst-case scenario is we jump back almost a whole month and forwards a whole month
# (i.e., 22320 * 3 timesteps)
# Here, we review the number of time steps and check the number of batches is sufficient
nt <- pbapply::pbsapply(iteration_julia$file_timeline, \(f) nrow(arrow::read_feather(f)))
range(nt)
stopifnot(min(nt) > 20000 & max(nt) < 22320 * 3)
nb <- p_batch(p_mem(2000, max(nt), 4, 100), 50e3)
range(nb)
stopifnot(nb <= iteration_julia$n_batch[1])

#### Check total size of input directories
# For sim, with write_feather_compressed():
# * 362.6076 MB for all iterations
# For real, with write_feather_compressed():
# * 510.7291 MB for all iterations
dir_size(iteration_julia$folder_input[1])
dir_size(file.path("data", "input", analysis, subanalysis), recursive = TRUE)

#### Write iteration
# Write full dataset 
qs::qsave(iteration, here_input_analysis("iteration.qs"))
# Write subsetted dataset (julia == TRUE) for modelling 
write_feather_compressed(iteration_julia, 
                         here_input_analysis("iteration.feather"))
# Write subsampled datasets
if (analysis == "real") {
  
  # Sample 100 rows for the initial validation analysis 
  n    <- nrow(iteration_julia[sensitivity == "best", ])
  size <- 100L
  pos  <- sample.int(n, size)
  stopifnot(length(pos) == 100L)
  
  # Define an initial validation dataset
  iteration_1 <- lapply(split(iteration_julia, iteration_julia$sensitivity), function(d) {
    d[pos, ]
  }) |> rbindlist()
  
  # Define remaining dataset
  iteration_2 <- lapply(split(iteration_julia, iteration_julia$sensitivity), function(d) {
    d[-pos, ]
  }) |> rbindlist()

  # Check nrow
  stopifnot(sum(c(nrow(iteration_1), nrow(iteration_2))) == nrow(iteration_julia))
  
  # Write to file 
  qs::qsave(iteration_1, here_input_analysis("iteration-1.qs"))
  qs::qsave(iteration_2, here_input_analysis("iteration-2.qs"))
  write_feather_compressed(iteration_1, here_input_analysis("iteration-1.feather"))
  write_feather_compressed(iteration_2, here_input_analysis("iteration-2.feather"))
  
  # Checks
  nrow(arrow::read_feather(here_input_analysis("iteration-1.feather")))
  nrow(arrow::read_feather(here_input_analysis("iteration-2.feather")))

}


###########################
###########################
#### Create chain files

# (optional) TO DO Move this code to run-patter.R 
# where other chain files are created

# ~ 1.9 mins (1 cl), 03m 40s
# ~ 0.5 mins (10 cl, chunk = TRUE)

if (analysis == "sim") {
  
  chains    <- copy(iteration)
  overwrite <- FALSE
  if (!all(file.exists(iteration_julia$file_residency_sim)) | overwrite) {
    
    pbo <- pbapply::pboptions(nout = 2L)
    cl_lapply(
      split(chains, seq_len(nrow(chains))), 
      .cl = 10L,
      .chunk = TRUE,
      .fun = function(d) {
        
        # Read maps
        # d <- chains[1, ]
        print(d$index)
        .map     <- terra::rast(here_input("map.tif"))
        .regions <- terra::rast(here_input("regions.tif"))
        
        # Define file_path_sim
        path <- paths[path_id == d$individual_id, ]
        qs::qsave(path, d$file_path_sim)
        
        # Define file_occupancy_sim
        path[, cell := terra::cellFromXY(.map, cbind(x, y))]
        if (any(is.na(path$cell))) {
          # Check for or filter NA cells
          # This is presumably a floating point issue on siam-linux20
          msg <- glue::glue("For path$path_id {path$path_id[1]}, there are {length(which(is.na(path$cell)))} NA cell(s)!")
          stop(msg)
          # warning(msg)
          path <- path[!is.na(cell)]
        }
        occupancy_sim <- map_pou(.map, path, .plot = FALSE)$ud
        terra::writeRaster(occupancy_sim, d$file_occupancy_sim, overwrite = TRUE)
        
        # Define file_residency_sim 
        residency_sim <- 
          occupancy_sim |> 
          terra::zonal(.regions, fun = "sum", na.rm = TRUE) |> 
          lazy_dt() |> 
          mutate(path_id = path$path_id[1]) |>
          select("path_id", region = "map_value", estimate = "map_value.1") |>
          as.data.table()
        
        # Check residency
        # * If needed, use reduced tolerance to handle floating point issue
        # * This is necessary if simulations are run on SIA-LAVENDED and copied onto server
        # * This does not seem to be necessary if simulations & this code are run on same machine
        stopifnot(isTRUE(all.equal(sum(residency_sim$estimate), 1))) # tolerance = 0.01
        qs::qsave(residency_sim, d$file_residency_sim)
        
      })
    pbapply::pboptions(pbo)
    
    #### Spot checks
    if (analysis == "sim") {
      stopifnot(all(file.exists(iteration_julia$file_path_sim)))
      stopifnot(all(file.exists(iteration_julia$file_occupancy_sim)))
      stopifnot(all(file.exists(iteration_julia$file_residency_sim)))
      # lapply_qplot_sim(iteration_julia, .n_plot = 4L)
    }
    
  }
  
}


#### End of code. 
###########################
###########################