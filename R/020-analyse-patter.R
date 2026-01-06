###########################
###########################
#### analysis-patter.R

#### Aims
# 1) This script provides generic analysis of patter outputs

#### Prerequisites
# 1) Run previous scripts


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
here_fig_analysis   <- switch_here_fig_analysis_subanalysis(analysis, subanalysis)
iteration           <- qs::qread(here_input_analysis("iteration.qs"))


###########################
###########################
#### Identify convergence 

#### Compute total run time
# Read callstats
iteration <- iteration[file.exists(file_callstats), ]
callstats <- lapply(iteration$index, function(i) {
  iteration$file_callstats[iteration$index == i] |> 
    arrow::read_feather() |> 
    mutate(index = i, .before = 1L) |> 
    cbind(iteration[index == i, .(individual_id, time_id, sensitivity)]) |> 
    as.data.table()
  }) |> rbindlist()
# Compute total run time
sum(callstats$time)

#### Compute total output size (MB, GB)
# Compute folder sizes
iteration[, folder_output_mb := 
            sapply(seq_len(nrow(iteration)),
                   \(i) dir_size(iteration$folder_output[i], recursive = TRUE))]
# Check total size (GB)
sum(iteration$folder_output_mb) / 1e3

#### Identify convergence
# Collate diagnostics
diagnostics <- cl_lapply(iteration$index, function(i) {
  iteration$file_diagnostics[iteration$index == i] |> 
    arrow::read_feather() |> 
    mutate(index = i, .before = 1L) |> 
    cbind(iteration[index == i, .(individual_id, time_id, sensitivity)]) |> 
    as.data.table()
}) |> rbindlist()
# Define convergence for an example individual
diagnostics |> 
  filter(index == 1L) |> 
  filter(routine == "smoother: two-filter") |> 
  summarise(prop = length(which(!is.na(ess))) / n())
# Define convergence 
convergence <- 
  diagnostics |> 
  lazy_dt() |> 
  group_by(index) |> 
  summarise(
    # Both filters converged if diagnostics smoother: two filter included
    pass_filter = any(routine %in% "smoother: two-filter"), 
    # Compute proportion of successful smoothing runs
    # (Do not use n() here which does not account for filtering)
            pass_smoother = 
              pick(routine, ess) |> 
              filter(routine == "smoother: two-filter") |>
              summarise(prop = length(which(!is.na(ess))) / 
                          length(which(routine == "smoother: two-filter"))) |> 
              pull(prop)
  ) |> 
  mutate(success = pass_filter & (pass_smoother >= 0.75)) |> 
  left_join(iteration[, .(index, individual_id, time_id, sensitivity)], by = "index") |> 
  as.data.table()
# Check convergence
convergence
table(convergence$success)
convergence[success == FALSE, ]

#### Filter by convergence (for subsequent steps)
# Get successful indices
successful_indices <- convergence$index[convergence$success == TRUE]
# Filter iteration
iteration <- iteration[index %in% successful_indices, ]
# Filter callstats
callstats <- callstats[index %in% successful_indices, ]
# Filter diagnostics
diagnostics <- diagnostics[index %in% successful_indices, ]
# > Subsequent analyses focus on successful runs only 


###########################
###########################
#### Summarise callstats

#### Compute computation time (mins)
# Summary statistics
callstats |> 
  group_by(routine) |> 
  reframe(utils.add::basic_stats(time / 60))
# Visualisation (~1 s)
tic()
png(here_fig_analysis("computation-time.png"), 
    height = 4, width = 5, units = "in", res = 800)
p <- 
  callstats |> 
  mutate(routine = stringr::str_to_sentence(routine)) |> 
  as_tibble() |> 
  ggplot(aes(routine, time / 60, fill = routine)) + 
  geom_violin() + 
  xlab("Routine") + ylab("Computation time (mins)") + 
  guides(fill = "none") +
  theme_bw()
print(p)
dev.off()
toc()


###########################
###########################
#### Summarise diagnostics

#### Summarise ESS (for filters & smoother)
# Summary statistics 
diagnostics |> 
  group_by(routine) |> 
  reframe(utils.add::basic_stats(ess, na.rm = TRUE))
# Visualisation (~14 s)
# TO DO Fix scientific notation here
# (This is ignored on the final panel)
tic()
png(here_fig_analysis("ess.png"), 
    height = 5, width = 10, units = "in", res = 800)
p <-
  diagnostics |> 
  mutate(routine = stringr::str_to_sentence(routine)) |> 
  as_tibble() |> 
  ggplot(aes(ess, fill = routine)) + 
  geom_density() + 
  xlab("Effective sample size") + ylab("Kernel density") + 
  guides(fill = "none") +
  facet_wrap(~routine, nrow = 1, scales = "free") +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0), labels = prettyGraphics::sci_notation) +
  theme_bw() + 
  theme(plot.margin = margin(t = 10, r = 20, b = 10, l = 20, unit = "pt"))
print(p)
dev.off()
toc()

#### Summarise area spanned by 95 % of the distribution
# TO DO, modify code above


###########################
###########################
#### Visualise maps

#### Visualise example maps
# TO DO


#### End of code. 
###########################
###########################