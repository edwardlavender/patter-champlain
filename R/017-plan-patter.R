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
it <- iteration[1, ]
list.files(it$folder_output)

#### Run algorithms 
# Run 001-run-algorithms.jl for example individual
# Check file sizes (MB)
if (FALSE) {
  dir_size(it$folder_output)
  file.size(it$file_callstats) / 1e6 
  file.size(it$file_diagnostics) / 1e6
  sum(file.size(list.files(it$folder_output, pattern = "\\.jld2$", full.names = TRUE))) / 1e6
}

#### Collate states 
# Run 001-run-algorithms.jl for example individual & check file size for it$file_states (MB)
if (FALSE) {
  file.size(it$file_states) / 1e6
}


###########################
###########################
#### Examine outputs

#### Read data 
it <- iteration[1, ]
it_callstats   <- arrow::read_feather(it$file_callstats)
it_diagnostics <- arrow::read_feather(it$file_diagnostics)
it_states      <- arrow::read_feather(it$file_states)
# Check file sizes
# > We should be able to read all diagnostic files into memory 
nrow(it_diagnostics) # 66960
lobstr::obj_size(it_diagnostics) 
nrow(it_states)

#### Examine computation time (SIA-LAVENDED, 50000/1500 particles, 1 thread, 10 batches)
## Check outputs
it_callstats
sum(it_callstats$time)
#
## Simulations: 
#
# For iteration[1, ]:
#
# timestamp           routine              n_particle n_iter   loglik convergence  time
# <dttm>              <chr>                     <int>  <dbl>    <dbl> <lgl>       <dbl>
# 1 2026-01-02 10:46:43 filter: forward           50000      1 -108537. TRUE         910.
# 2 2026-01-02 11:01:53 filter: backward          50000      1 -106402. TRUE         642.
# 3 2026-01-02 11:12:36 smoother: two-filter       1500    NaN     NaN  TRUE        3433.
#
# > Total : 4984.926 s (1.38 hours) for iteration[1, ]
# > ETA   : 0.60 days for 210 iterations on 20 cores
(4984.926 * 210)  / 60 / 60 / 24 / 20
#
## Real-world analysis
# * TO DO for it[1, ] 
# * >>8 days (ETA) for 2723 iterations on 20 cores
(4984.926 * 2723)  / 60 / 60 / 24 / 20

#### Examine convergence 
#
## For sim iteration[1, ]:
# percent_valid convergence
# <num>      <lgcl>
#   1:      99.98656        TRUE
#
# For real iteration[1, ]:
#
# TO DO
# 
it_callstats$convergence[1]
it_callstats$convergence[2]
it_diagnostics |>
  filter(routine == "smoother: two-filter") |> 
  summarise(percent_valid = length(which(!is.na(ess))) / n() * 100, 
            convergence = percent_valid > 75) |> 
  as.data.table()

#### Examine ESS
#
## For sim it[1, ]:
#
# routine                min   mean median    max     sd    IQR    MAD
# <chr>                <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
#   1 filter: backward      152. 14137.  8793. 49933. 13181. 17798. 9236. 
# 2 filter: forward       277. 14095.  8833. 49935. 13123. 17606. 9285. 
# 3 smoother: two-filter    1    214.   182.  1500    159.   111.   82.8
#
## For real it[1, ]
# 
# TO DO
#
it_diagnostics |> 
  group_by(routine) |> 
  reframe(utils.add::basic_stats(ess, na.rm = TRUE))

#### Map states (~7 s)
# Each map requires 0.069128 MB of storage
tic()
tmp.tif <- tempfile(fileext = ".tif")
tmp     <- map_pou(map, it_states)$ud
terra::writeRaster(tmp, tmp.tif, overwrite = TRUE)
file.size(tmp.tif) / 1e6
# map_dens(map, it_states, .discretise = TRUE)
toc()

#### Map states + trajectory
# This is a quick validation analysis
# The map & simulated trajectory should align
if (analysis == "sim") {
  # Read paths 
  paths <- qs::qread(here_input_analysis("paths.qs"))
  path  <- paths[path_id == it$individual_id[1], ]
  # Plot map & paths
  terra::plot(tmp)
  patter:::add_sp_path(path$x, path$y, length = 0.01, lwd = 0.1)
}


###########################
###########################
#### Define implementation strategy 

#### Define the number of particles for algorithms
# We run algorithms on {n_cpu}
# We record {n_record} particles
# We quadrouple the required memory
# - We may start a new run while still finishing another (x2)
# - We need a safety buffer (x3)
# We use a maximum of {mb_available} memory in total
n_record     <- 1500
n_cpu        <- 20
mb_available <- 25e3
mb_required  <- p_mem(n_record, 22320, 4, n_cpu) * 4
p_batch(mb_required, mb_available)

#### Storage requirements for simulations (1500 particles):
# Breakdown (Arrow.write default, Arrow.write compress = ZstdCompressor(level = 9))
# * callstats   : 0.00189,  0.00209
# * diagnostics : 3.550282, 0.896738
# * jld2 files  : 1071.414, 1071.414 (same)
# * states      : 1874.882, 678.995
# * map         : 0.068324, 0.068324 (see below)
# * Total       : 1878.502, 679.9622
#   (excluding jld2 files, which we delete)
# To store the outputs of the simulations, we need ~142.7921 GB
(0.00209 + 0.896738 + 678.995 + 0.068324) * 210 / 1e3 

#### Storage requirements for real-world analysis (2000 particles): 
# To store the outputs of the real world analysis, we would need 2468.654 GB
z <- 2000/1500
(0.00209 * z + 0.896738 * z + 678.995 * z + 0.068324) * 2723 / 1e3 

#### Define implementation strategy
# Handle iterations in batches e.g., 1:100, 101:200, ..., N
# For each batch:
# * Run particle algorithms in parallel (Julia)
# * Collate states (Julia)
# * Make maps & compute diagnostics (R)
# * Unlink large particle files (R) 
# This minimises disk space requirements
# * Each batch requires only ~200 GB of storage space: 1878.502 * 100 / 1e3 
# * All batches ultimately require only:
# - (0.00189 + 3.550282 + 0.068324) * 210 = 760 MB (simulations)
# - (0.00189 + 3.550282 + 0.068324) * 2723 = 9858.611 MB (real-world analysis)


#### End of code. 
###########################
###########################