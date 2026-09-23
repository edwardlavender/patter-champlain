###########################
###########################
#### plan-patter.R

#### Aims
# 1) This script is used to plan patter runs
#    * Examine outputs for an example individual
#    * Check file sizes
#    * Define implementation strategy 

#### Prerequisites
# 1) Run simulations
# 2) This code is designed for SIA-LAVENDED 


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Set global options
Sys.setenv(JULIA_SESSION = "FALSE")
set.seed(123L)

#### Load essential packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(glue)
library(JuliaCall)
library(proj.verse)
library(tictoc)
files_source_r(here_src())

#### Load data
map <- terra::rast(here_input("map.tif"))


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
#### Generate example outputs

#### Define example individual
i  <- 1
it <- iteration[i, ]
list.files(it$folder_output_block)

#### Run algorithms 
# Run 001-run-algorithms.jl for example individual
# Check file sizes (MB)
if (FALSE) {
  dir_size(it$folder_output_block)
  file.size(it$file_callstats) / 1e6 
  file.size(it$file_diagnostics) / 1e6
  sum(file.size(list.files(it$folder_output_block, pattern = "\\.jld2$", full.names = TRUE))) / 1e6
  sum(file.size(list.files(it$folder_output_block, pattern = "pou", full.names = TRUE)) / 1e6)
}

#### Collate states 
# Run 001-run-algorithms.jl for example individual & check file size for it$file_states (MB)
if (FALSE) {
  
  # Connect to Julia 
  # * Connecting to Julia can cause issues with geospatial operations below (even on SIA-LAVENDED)
  # * If needed:
  #   - Implement this code
  #   - Restart R
  #   - Continue analysis below using produced file (it$file_states)
  patter::julia_connect()
  
  # Make it$file_states & compute MB (~11 s) 
  tic()
  julia_source("./Julia/src/utils.jl")
  julia_assign("iteration", iteration)
  julia_command(glue('iter = iteration[{i}, :]'))
  julia_command('write_file_states(iter)')
  file.size(it$file_states) / 1e6 
  toc()
  
  # For sim iteration[1, ], there is little benefit in rewriting as .qs:
  it_states <- arrow::read_feather(it$file_states)
  tmp.qs <- tempfile(fileext = ".qs")
  qs::qsave(it_states, tmp.qs, preset = "custom", algorithm = "zstd", compress_level = 9L) #  628.6245 MB
  qs::qsave(it_states, tmp.qs, preset = "archive") # 625.6087
  file.size(tmp.qs) / 1e6
  #  qs::qread(tmp.qs)
}


###########################
###########################
#### Read outputs 

#### Read outputs 
it <- iteration[1, ]
it_timeline    <- arrow::read_feather(it$file_timeline)$timestamp
it_callstats   <- arrow::read_feather(it$file_callstats) |> setDT()
it_diagnostics <- arrow::read_feather(it$file_diagnostics) |> setDT()
it_states      <- NULL
if (file.exists(it$file_states)) {
  it_states    <- arrow::read_feather(it$file_states) |> setDT()
}

# Check object sizes
# > We should be able to read all diagnostic files into memory 
nrow(it_diagnostics)             # 66960 rows
lobstr::obj_size(it_diagnostics) # 1.34-3.5 MB
# nrow(it_states)                # 44640000 rows


###########################
###########################
#### Examine computation time

#### Check callstats
# sim: SIA-LAVENDED, 50000/1500 particles, 1 thread, 10 batches
# real: SIA-LAVENDED, 75000/2000 particles, 1 thread, 10 batches
it_callstats
sum(it_callstats$time)

#### Simulations, iteration[1, ]:

# timestamp           routine              n_particle n_iter   loglik convergence  time
# <dttm>              <chr>                     <int>  <dbl>    <dbl> <lgl>       <dbl>
# 1 2026-01-02 10:46:43 filter: forward           50000      1 -108537. TRUE         910.
# 2 2026-01-02 11:01:53 filter: backward          50000      1 -106402. TRUE         642.
# 3 2026-01-02 11:12:36 smoother: two-filter       1500    NaN     NaN  TRUE        3433.

# > Total : 4984.926 s (1.38 hours) for iteration[1, ]
# > ETA   : 0.60 days for 210 iterations on 20 cores
(4984.926 * 210)  / 60 / 60 / 24 / 20

#### Real-world analysis, iteration[1, ]

# timestamp           routine              n_particle n_iter   loglik convergence  time
# <dttm>              <chr>                     <int>  <dbl>    <dbl> <lgl>       <dbl>
# 1 2026-01-02 12:33:19 filter: forward           75000      1 -195703. TRUE        1069.
# 2 2026-01-02 12:51:08 filter: backward          75000      1 -195507. TRUE        1079.
# 3 2026-01-02 13:09:07 smoother: two-filter       2000    NaN     NaN  TRUE        6631.

# > Total : 8779.24 s (2.438678 hours) for iteration[1, ]
# > ETA   : 13.83442 days for 2723 iterations on 20 cores
(8779.24 * 2723)  / 60 / 60 / 24 / 20


###########################
###########################
#### Examine convergence 

#### For sim iteration[1, ]:
# percent_valid convergence
# <num>      <lgcl>
#  99.98656        TRUE

#### For real iteration[1, ]:
# percent_valid convergence
# <num>      <lgcl>
#  100        TRUE

it_callstats$convergence[1]
it_callstats$convergence[2]
it_diagnostics |>
  filter(routine == "smoother: two-filter") |> 
  summarise(percent_valid = length(which(!is.na(ess))) / n() * 100, 
            convergence = percent_valid > 75) |> 
  as.data.table()


###########################
###########################
#### Examine ESS

#### For sim iteration[1, ]:
# routine                min   mean median    max     sd    IQR    MAD
# <chr>                <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
# 1 filter: backward      152. 14137.  8793. 49933. 13181. 17798. 9236. 
# 2 filter: forward       277. 14095.  8833. 49935. 13123. 17606. 9285. 
# 3 smoother: two-filter    1    214.   182.  1500    159.   111.   82.8

#### For real iteration[1, ]:
# routine                 min   mean median    max     sd    IQR    MAD
# <chr>                 <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
# 1 filter: backward     380.   20271. 11796. 74866. 19741. 26389. 12700.
# 2 filter: forward      401.   20498. 12186. 74769. 19864. 26697. 13239.
# 3 smoother: two-filter   1.63   694.   703.  2000    376.   643.   469.

it_diagnostics |> 
  group_by(routine) |> 
  reframe(utils.add::basic_stats(ess, na.rm = TRUE))


###########################
###########################
#### Mapping 

#### Map states 
if (Sys.getenv("JULIA_SESSION") == "FALSE") {
  
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
  stopifnot(isTRUE(all.equal(1, sum(coord$mark))))

  # Map occupancy 
  occupancy <- terra::rasterize(coord, map, values = coord$mark)
  occupancy <- terra::classify(occupancy, cbind(NA, 0))
  occupancy <- terra::mask(occupancy, map)
  names(occupancy) <- "map_value"
  terra::plot(occupancy)
  
  # Each map requires 0.069128 MB of storage
  occupancy.tif <- tempfile(fileext = ".tif")
  terra::writeRaster(occupancy, occupancy.tif, overwrite = TRUE)
  file.size(occupancy.tif) / 1e6
  # patter::map_dens(map, it_states, .discretise = TRUE)
  
  # Validate occupancy map with via patter::map_pou()
  # > The map produced from the summarised Julia outputs
  #   should match that produced using all states from patter
  if (!is.null(it_states)) {
    # Use map_pou() to verify occupancy 
    occupancy_2 <- patter::map_pou(map, it_states, .plot = FALSE)$ud
    stopifnot(isTRUE(terra::all.equal(occupancy, occupancy_2, maxcell = terra::ncell(occupancy))))
    # Visual check
    pp <- par(mfrow = c(1, 2))
    terra::plot(occupancy)
    terra::plot(occupancy_2)
    par(pp)
  }
  
}

#### Map states + trajectory
# This is a quick validation analysis
# The map & simulated trajectory should align
if (analysis == "sim") {
  # Read paths 
  paths <- qs::qread(here_input_analysis("paths.qs"))
  path  <- paths[path_id == it$individual_id[1], ]
  # Plot map & paths
  terra::plot(occupancy)
  patter:::add_sp_path(path$x, path$y, length = 0.01, lwd = 0.1)
}


###########################
###########################
#### Define implementation strategy 

#### Define the number of particles for algorithms
# We run algorithms on {n_cpu}
# We record {n_record} particles
# We quadruple the required memory
# - We may start a new run while still finishing another (x2)
# - We need a safety buffer (x3)
# We use a maximum of {mb_available} memory in total
n_record     <- 2000
n_cpu        <- 20
mb_available <- 30e3
mb_required  <- p_mem(n_record, 22320, 4, n_cpu) * 3
p_batch(mb_required, mb_available)

#### Storage requirements for simulations (1500 particles):
# Breakdown (Arrow.write default, Arrow.write compress = ZstdCompressor(level = 9))
# * callstats   : 0.00189,  0.00209
# * diagnostics : unknown,  0.955378 (formerly 0.896738 without two area metrics)
# * jld2 files  : 1071.414, 1071.414 (same)
# * states      : 1874.882, 678.995
# * pou files   : unknown,  24.10011
# * map         : 0.068324, 0.068324 (see below)
# * Total       : 1878.502, 679.9622
#   (excluding jld2 files, which we delete)
# To store all outputs of the simulations, we would need ~142.7921 GB
(0.00209 + 0.955378 + 678.995 + 0.068324) * 210 / 1e3 

#### Storage requirements for real-world analysis (2000 particles): 
# Breakdown (Arrow.write compress = ZstdCompressor(level = 9))
# * callstats   : 0.00209  
# * diagnostics : unknown (formerly 1.071778)
# * jld2 files  : 1428.534
# * states      : 927.1699
# * pou files   : TO DO (unknown)
# * map         : 0.036435
# * Total       : 928.2802
#   (excluding jld2 files, which we delete)
# To store all outputs of the real world analysis, we would need ~2527.707 GB
# A) Guess based on simulation analysis with 1500 versus 2000 particles:
z <- 2000/1500
(0.00209 * z + 0.896738 * z + 678.995 * z + 0.068324) * 2723 / 1e3 
# B) Estimate based on real iteration[1, ]
(0.00209 + 1.071778 + 927.1699 + 0.036435) * 2723 / 1e3 

#### Define implementation strategy
# Run iterations in parallel 
# For each iteration
# * Run particle algorithms using batching & record:
#   - callstats
#   - diagnostics, including area spanned by smoothed distributions
#   - pou data.tables by batch, which we can use to collate maps in R
# * Unlink large particle files (R) 
# This minimises memory & disk space requirements for effective parallelisation
# - Particle outputs are only produced temporarily
# - We analyse outputs immediately without reading all particle into memory
# - We delete particle files on the fly after analysis 
# - Later in R, we analyse diagnostics & collate maps using smaller outputs

#### Summarise storage requirements (GB) for implementation strategy

## (1) Temporary storage requirements
# We need temporary storage for all CPUs producing:
# - callstats.feather
# - diagnostics.feather
# - fwd-{i}.jld2, bwd-{i}.jld2, smo-{i}.jld2 matrices
# - pou-{i}.feather files
(0.00209 + 0.955378 + 1071.414 * 3 + 24.10011) * n_cpu / 1e3 # sim, 64.78599 GB
(0.00209 + 0.955378 + 1428.534 * 3 + 24.10011) * n_cpu / 1e3 # real, 86.21319 GB

## (2) Ultimate Julia iteration requirements (GB)
# All Julia iterations ultimately require storage for:
# - callstats.feather
# - diagnostics.feather
# - pou-{i}.feather files
(0.00209 + 0.955378 + 24.10011) * 210 / 1e3  # sim, 5.262091 GB
(0.00209 + 0.955378 + 24.10011) * 2723 / 1e3 # real, 68.23178 GB

## (3) Ultimate iteration requirements after R mapping 
# These are as in (2) but we include occupancy.tif files
# (We could clean up pou-{i}.feather files)
(0.00209 + 0.955378 + 24.10011 + 0.068324) * 210 / 1e3  # sim, 5.276439 GB
(0.00209 + 0.955378 + 24.10011 + 0.036435) * 2723 / 1e3 # real, 68.331 GB


#### End of code. 
###########################
###########################