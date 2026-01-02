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
library(patter)
library(proj.verse)
library(spatial.extensions)
library(tictoc)
library(truncdist)
files_source_r(here_src())

#### Load data
map       <- terra::rast(here_input("map.tif"))
champlain <- qreadvect(here_input("champlain-utm.qs"))
pars      <- qs::qread(here_input("pars-patter.qs"))


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
if (analysis == "real") {
  detections_raw <- qs::qread(here_input_analysis("detections-raw.qs"))
}


###########################
###########################
#### Define detection dataset (blocks)

#### Update detections with individual/time_id blocks
# Define individual_id, time_id blocks
detections <- 
  detections |> 
  group_by(individual_id) |> 
  mutate(time_id = lubridate::floor_date(timestamp, "months")) |>
  select(individual_id, time_id, timestamp, receiver_id) |>
  as.data.table()
# Check number of individuals in full dataset
length(unique(detections$individual_id))

#### Focus on individual/time (month) units with sufficient data
# NB: filter_detectionsd assumes monthly blocks
detections <- filter_detections(detections)

#### Checks
# Validate time_id assignment in detections
stopifnot(all(detections$time_id == 
                lubridate::floor_date(detections$timestamp, "months")))
# Check the number of days with detections meets minimum criteria
# * This is based on the threshold specified in filter_detections.R
ck <- 
  detections |> 
  group_by(paste(individual_id, time_id)) |> 
  summarise(ck = length(unique(lubridate::yday(timestamp)))) |> 
  pull(ck) |> 
  sort()
stopifnot(all(ck >= 14))
# Check the number of detections
# * We want to catch time series with 'too few' observations
# * What is 'too few' here is somewhat arbitrary
# * The goal is to catch potential mistakes in data processing
# * E.g., that would otherwise allow 1 or 2 row datasets forward for analysis
ck <- 
  detections |> 
  group_by(paste(individual_id, time_id)) |> 
  summarise(ck = n()) |> 
  pull(ck) |> 
  sort()
stopifnot(all(ck >= 50))

#### Summarise (real) detection dataset used for modelling 
# cf. raw data summary statistics (setup-data-detection.R)
if (TRUE) {
  nrow(detections)
  length(unique(detections$individual_id))
  range(detections$timestamp)
  int <- lubridate::interval(min(detections$timestamp),max(detections$timestamp))
  lubridate::time_length(int, "months")
  lubridate::time_length(int, "years")
  length(unique(detections$receiver_id))
  length(unique(paste(detections$individual_id, detections$time_id)))
  plot(detections$timestamp, detections$individual_id)
}

#### Visualise real detection dataset used for modelling 
# (optional) TO DO Move this code to appropriate synthesis script
# Plot raw time series (light grey)
# Add modelled time series, coloured by region as in map
if (analysis == "real") {
  
  #### Define colour scheme
  # Define scheme 
  cols <- tribble(
    ~region,            ~col,
    "Main Lake Central", "#6c81de", 
    "Main Lake North",   "#9834df", 
    "Main Lake South",   "#e41ea5", 
    "Malletts Bay",      "#df756d", 
    "Missisquoi Bay",    "#87e93b", 
    "Northeast Arm",     "#22e45f", 
    "South Lake",        "#0dcbd9")
  # Add to moorings
  moorings <- 
    moorings |> 
    mutate(region = terra::extract(champlain, cbind(receiver_x, receiver_y))$region, 
           col = cols$col[match(region, cols$region)]) |> 
    as.data.table()
  stopifnot(all(!is.na(moorings$col)))
  
  #### Process raw detections
  ids  <- unique(detections_raw$individual_id)
  draw <- copy(detections_raw)
  draw[, model := FALSE]
  draw[detections, on = .(individual_id, timestamp), model := TRUE]
  draw[, individual_id := factor(individual_id, levels = ids)]
  draw[, col := moorings$col[match(receiver_id, moorings$receiver_id)]]
  draw[, col := ifelse(model == TRUE, col, scales::alpha("dimgrey", 0.25))]
  
  #### Plot (~15 s)
  tic()
  range(draw$timestamp)
  png(here_fig(analysis, "detections.png"), 
      height = 9.69 * 1.75, width = 6.27 * 1.75, units = "in", res = 800)
  # Set parameters 
  pp <- par(oma = c(1.5, 1.5, 0, 0))
  xshift    <-  5 * 24 * 60 * 60
  cex.axis  <- 2
  cex.mtext <- 2.25
  # Blank plot
  plot(draw$timestamp, draw$individual_id, 
       type = "n",
       xlim = c(min(draw$timestamp) - xshift, max(draw$timestamp) + xshift),
       ylim = c(0, length(ids) + 1),
       xaxs = "i", yaxs = "i",
       xlab = "", ylab = "", 
       xaxt = "n", cex.axis = cex.axis, las = TRUE)
  # Add x grid (by month)
  months <- seq(lubridate::floor_date(min(draw$timestamp), "months"), 
                lubridate::floor_date(max(draw$timestamp), "months"),
                by = "months")
  sapply(months, \(month) abline(v = month, col = "lightgrey", lty = 3)) |> invisible()
  # Add y grid (by individual)
  sapply(unique(draw$individual_id), \(id) abline(h = id, col = "lightgrey")) |> invisible()
  # Add points (on top of grid)
  n <- nrow(draw)
  # n <- 1e5
  points(draw$timestamp[1:n], draw$individual_id[1:n],
         pch = 3, col = draw$col, las = TRUE)
  # Add axes
  xat <- as.POSIXct(paste0(rep(2014:2017, each = 2), c("-01-01", "-06-01")), tz = "UTC")
  axis(side = 1, at = xat, labels = format(xat, "%b-%y"), cex.axis = cex.axis)
  mtext(side = 1, "Time (month-year)", line = 4, cex = cex.mtext)
  mtext(side = 2, "Individual", line = 4, cex = cex.mtext)
  par(pp)
  dev.off()
  toc()
  
}


###########################
###########################
#### Define unitsets

#### Define unitsets
unitsets <- 
  detections |> 
  group_by(individual_id, time_id) |> 
  slice(1L) |> 
  ungroup() |> 
  arrange(individual_id, time_id) |> 
  mutate(unit_id = row_number()) |> 
  select(unit_id, individual_id, time_id, timestamp) |>
  as.data.table()
# Update detections with unit_id
detections <- 
  detections |> 
  left_join(unitsets[, .(unit_id, individual_id, time_id)], 
            by = c("individual_id", "time_id")) |> 
  select(unit_id, individual_id, time_id, timestamp, receiver_id) |> 
  as.data.table()

#### Checks
# Visually validate matching between unitsets & detections
unitsets[unit_id == 14, ]
detections[unit_id == 14, ]
# Validate all unit_ids present in each dataset
stopifnot(all(unitsets$unit_id %in% detections$unit_id) & 
            all(detections$unit_id %in% unitsets$unit_id))
# Validate matching between all unitsets and detections
cl_lapply(split(unitsets, seq_len(nrow(unitsets))), function(sim) {
  vdetections <- detections[unit_id == sim$unit_id, ]
  stopifnot(all(sim$unit_id == vdetections$unit_id))
  stopifnot(all(sim$individual_id == vdetections$individual_id))
  stopifnot(all(sim$time_id == vdetections$time_id))
})


###########################
###########################
#### Prepare iterations

#### Define iteration 
iteration <- 
  unitsets |> 
  cross_join(pars) |> 
  mutate(
    # Define input files 
    # * Some files depend on both unit_id & sensitivity parameters
    # * For convenience, we store all files in an {individual_id}/{unit_id}/{parameter_id} directory 
    folder_input        = file.path("data", "input", analysis, subanalysis, "runs", 
                                    individual_id, time_id, parameter_id), 
    file_timeline       = file.path(folder_input, "timeline.feather"),
    file_acoustics      = file.path(folder_input, "acoustics.feather"),
    file_containers_fwd = file.path(folder_input, "containers-fwd.feather"),
    file_containers_bwd = file.path(folder_input, "containers-bwd.feather"),
    # Define  output files (unit-specific & sensitivity specific)
    # * We use .feather to record outputs
    # * We can write these from Julia & read them into R correctly
    folder_output       = file.path("data", "output", analysis, subanalysis, "runs", 
                                    individual_id, time_id, parameter_id),
    file_states         = file.path(folder_output, "states.feather"),
    file_diagnostics    = file.path(folder_output, "diagnostics.feather"),
    file_callstats      = file.path(folder_output, "callstats.feather"),
    file_particles      = file.path(folder_output, "particles.feather"),
    # Add modelling columns
    n_batch             = 10L,
    n_particle_filter   = ifelse(analysis == "sim", 50000L, 75000L), 
    n_particle_smoother = ifelse(analysis == "sim", 1500L, 2000L), 
  ) |> 
  as.data.table()

#### Check nrow
# This must be feasible! 
nrow(iteration)

#### Build directories 
if (FALSE) {
  unlink(iteration$folder_input, recursive = TRUE)
  unlink(iteration$folder_output, recursive = TRUE)
}
dirs.create(iteration$folder_input)
dirs.create(iteration$folder_output)


###########################
###########################
#### Create iteration input files

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

#### Write files 
pbo <- pbapply::pboptions(nout = 2L)
cl_lapply(split(iteration, seq_len(nrow(iteration))), 
          .cl = 10L,
          .chunk = TRUE,
          .fun = function(d) {
  
  ## Define file_timeline
  dets     <- detections[unit_id == d$unit_id, ]
  timeline <- seq(dets$time_id[1], 
                  lubridate::ceiling_date(max(dets$timestamp), "months") - 60 * 2, 
                  by = "2 mins")
  stopifnot(length(timeline) > 20000 & length(timeline) < 30000)
  timeline <- data.table(timestamp = timeline)
  write_feather_compressed(timeline, d$file_timeline)
  
  ## Define file_acoustics
  # Define moorings, with detection probability parameters
  moors <- 
    moorings |> 
    lazy_dt(immutable = TRUE) |> 
    mutate(receiver_alpha = d$receiver_alpha, 
           receiver_beta = d$receiver_beta, 
           receiver_gamma = d$receiver_gamma) |> 
    as.data.table()
  # Define acoustics 
  accs <- assemble_acoustics(.timeline = timeline$timestamp, .detections = dets, .moorings = moors)
  write_feather_compressed(accs, d$file_acoustics)
  
  ## Define acoustic containers (file_containers_fwd, file_containers_bwd)
  containers <- assemble_acoustics_containers(.timeline = timeline$timestamp, 
                                              .acoustics = accs,
                                              .mobility = d$mobility, 
                                              .map = map)
  containers_fwd <- containers$forward
  containers_bwd <- containers$backward
  write_feather_compressed(containers_fwd, d$file_containers_fwd)
  write_feather_compressed(containers_bwd, d$file_containers_bwd)
  
  nothing()
})
pbapply::pboptions(pbo)

#### Check total size of input directories
# With write_feather_compressed():
# * 0.69476 MB per iteration
# * 184.3025 MB for simulations
# * 1060 for real-world analysis (estimated)
dir_size(iteration$folder_input[1])
dir_size(file.path("data", "input", analysis, subanalysis), recursive = TRUE)
2723 * 81.82257 / 210

#### Write iteration
qs::qsave(iteration, here_input_analysis("iteration.qs"))
write_feather_compressed(iteration, here_input_analysis("iteration.feather"))


#### End of code. 
###########################
###########################