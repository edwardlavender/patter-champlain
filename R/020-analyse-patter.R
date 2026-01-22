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
library(glue)
library(ggplot2)
library(proj.verse)
library(tictoc)
files_source_r(here_src())

#### Load data
map <- terra::rast(here_input("map.tif"))


###########################
###########################
#### Select analysis

#### Define analysis 
# analysis <- "sim"
analysis <- "real"
subanalysis <- "main"

#### Define analysis-specific data
here_input_analysis  <- switch_here_input_analysis_subanalysis(analysis, subanalysis)
here_output_analysis <- switch_here_output_analysis_subanalysis(analysis, subanalysis)
here_fig_analysis    <- switch_here_fig_analysis_subanalysis(analysis, subanalysis)
iteration            <- qs::qread(here_input_analysis("iteration.qs"))


###########################
###########################
#### Analyse trials

# Using a sample real-world time series, we check convergence via run-filter.jl.
# This script only produces a file_callstats_filter file that records convergence
# (unlike the full run-algorithms.jl script). 
# We check those results here. 
# We further dig into convergence failures in refine-patter.R. 
# Later sections in this script are designed to analyse the outputs of the 
# full workflow i..e, run-algorithms.jl. 

if (FALSE) {
  
  #### Read callstats for the filter
  table(file.exists(iteration$file_callstats_filter))
  iteration <- iteration[file.exists(file_callstats_filter), ]
  callstats <- lapply(iteration$index, function(i) {
    iteration$file_callstats_filter[iteration$index == i] |> 
      arrow::read_feather() |> 
      mutate(index = i, .before = 1L) |> 
      cbind(iteration[index == i, .(individual_id, time_id, sensitivity, sensitivity_label)]) |> 
      as.data.table()
  }) |> 
    rbindlist()
  
  #### Check computation time
  # Estimate average run time  (mins)
  utils.add::basic_stats(callstats$time / 60)
  # Estimate total time (mins) shared across n_cpu
  n_cpu <- 100L
  sum(callstats$time) / 60 / n_cpu 
  
  #### Identify convergence failures
  failures <- 
    callstats |> 
    filter(routine == "filter: forward") |> 
    filter(convergence == FALSE) |> 
    as.data.table()
  
  #### Summarise failures
  # Count failures
  nrow(failures)
  # Proportion of failures
  nrow(failures) / nrow(iteration)
  # Count failures by sensitivity 
  failures |> 
    group_by(sensitivity) |>
    summarise(n())
  # Proportion of failures by sensitivity 
  table(failures$sensitivity) / table(iteration$sensitivity)

  #### Check failures
  # All failures 
  failures
  # Failures for main analysis
  failures[sensitivity == "best", .(index, individual_id, time_id, sensitivity, routine, n_particle, time, convergence)]
  
  #### (optional) Record datasets
  today <- as.Date(Sys.time())
  now   <- as.numeric(Sys.time())
  debug <- here_output_analysis("debug", paste0(today, "-", now))
  # dir.create(debug, recursive = TRUE)
  # qs::qsave(callstats, file.path(debug, "callstats.qs"))
  # qs::qsave(failures, file.path(debug, "failures.qs"))
  
  #### Record log 
  # Check movement model 
  iter <- iteration[1, ]
  glue(
    '
    # Movement model formulation:
    ModelMoveCXY(env, 
                {iter$mobility}, 
                truncated(Gamma({iter$shape}, {iter$scale}), upper = {iter$mobility}), 
                MixtureModel([truncated(Normal(0.0, {iter$phi}), -pi, pi), Uniform(-pi, pi)], [0.99, 0.01]))
    
    # Observation model formulation:
    truncated(logistic({iter$receiver_alpha} + {iter$receiver_beta} * distance), {iter$receiver_gamma})
    '
  )
  # Record settings & outcome
  glue(
    '
    
    # {today}, {now}
    
    * n_particle = XXX
    * n_move = XXX
    * Proportion best failures: {length(which(failures$sensitivity == "best")) / length(which(iteration$sensitivity == "best"))}
    * Run time: {sum(callstats$time) / 60} mins 
    '
    ) |> 
    cat(file = here_doc("debug-batch.txt"), append = TRUE)
  
  # Cleanup
  # unlink(iteration$file_callstats_filter)

}


###########################
###########################
#### Compute resource requirements

#### Compute total run time
# Read callstats
iteration <- iteration[file.exists(file_callstats), ]
callstats <- lapply(iteration$index, function(i) {
  iteration$file_callstats[iteration$index == i] |> 
    arrow::read_feather() |> 
    mutate(index = i, .before = 1L) |> 
    cbind(iteration[index == i, .(individual_id, time_id, sensitivity, sensitivity_label)]) |> 
    as.data.table()
  }) |> 
  rbindlist() |> 
  mutate(routine_label = stringr::str_to_sentence(routine), .after = routine) |> 
  mutate(routine_label = factor(routine_label, levels = c("Filter: forward", 
                                                          "Filter: backward", 
                                                          "Smoother: two-filter"))) |> 
  as.data.table()

# Compute total run time (days)
# > sim     : 22.98759 / 20              # 1.14 days on 20 cores
# > ETA real: 22.98759 / 210 * 2723 / 20 # 14 days on 20 cores
sum(callstats$time) / 60 / 60 / 24

#### Compute total output size (MB, GB)
# Compute folder sizes
iteration[, folder_output_mb := 
            sapply(seq_len(nrow(iteration)),
                   \(i) dir_size(iteration$folder_output[i], recursive = TRUE))]
# Check total size (GB)
# > sim: 4.190048 GB
sum(iteration$folder_output_mb) / 1e3


###########################
###########################
#### Identify convergence

#### Collate diagnostics
diagnostics <- 
  cl_lapply(iteration$index, function(i) {
    iteration$file_diagnostics[iteration$index == i] |> 
      arrow::read_feather() |> 
      mutate(index = i, .before = 1L) |> 
      cbind(iteration[index == i, .(individual_id, time_id, sensitivity)]) |> 
      as.data.table()
  }) |> 
  rbindlist()  |> 
  mutate(routine_label = stringr::str_to_sentence(routine), .after = routine) |> 
  mutate(routine_label = factor(routine_label, levels = c("Filter: forward", 
                                                          "Filter: backward", 
                                                          "Smoother: two-filter"))) |> 
  as.data.table()

#### Define convergence for an example individual
diagnostics |> 
  filter(index == 1L) |> 
  filter(routine == "smoother: two-filter") |> 
  summarise(prop = length(which(!is.na(ess))) / n())

#### Define convergence 
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

#### Check convergence
convergence
table(convergence$success)
table(convergence$success, convergence$sensitivity == "best")
table(convergence$sensitivity[convergence$success == FALSE])
convergence[success == FALSE, ]
utils.add::basic_stats(convergence$pass_smoother, na.rm = TRUE)
utils.add::basic_stats(convergence$pass_smoother[convergence$success], na.rm = TRUE)

#### Review convergence failures

## sim:
# index pass_filter pass_smoother success individual_id    time_id sensitivity
# <int>      <lgcl>         <num>  <lgcl>         <int>     <POSc>      <char>
# 1:     7       FALSE           NaN   FALSE             1 2025-01-01       ac(+)
# 2:    13       FALSE           NaN   FALSE             2 2025-01-01       ac(-)
# 3:    14       FALSE           NaN   FALSE             2 2025-01-01       ac(+)
# 4:    63        TRUE     0.4146057   FALSE             9 2025-01-01       ac(+)
# 5:   119       FALSE           NaN   FALSE            17 2025-01-01       ac(+)
# 6:   133       FALSE           NaN   FALSE            19 2025-01-01       ac(+)
# 7:   140       FALSE           NaN   FALSE            20 2025-01-01       ac(+)
# 8:   147       FALSE           NaN   FALSE            21 2025-01-01       ac(+)
# 9:   161       FALSE           NaN   FALSE            23 2025-01-01       ac(+)
# 10:   166       FALSE           NaN   FALSE            24 2025-01-01    angle(+)
# 11:   168       FALSE           NaN   FALSE            24 2025-01-01       ac(+)
# 12:   176       FALSE           NaN   FALSE            26 2025-01-01        best
# 13:   178       FALSE           NaN   FALSE            26 2025-01-01     step(+)
# 14:   179       FALSE           NaN   FALSE            26 2025-01-01    angle(-)
# 15:   182       FALSE           NaN   FALSE            26 2025-01-01       ac(+)
# 16:   195       FALSE           NaN   FALSE            28 2025-01-01       ac(-)
# 17:   196       FALSE           NaN   FALSE            28 2025-01-01       ac(+)

## real:
# index pass_filter pass_smoother success individual_id    time_id sensitivity
# <int>      <lgcl>         <num>  <lgcl>         <int>     <POSc>      <char>
#   1:   421       FALSE           NaN   FALSE         24334 2015-05-01        best
# 2:   442        TRUE     0.6628352   FALSE         24334 2016-02-01        best
# 3:   596       FALSE           NaN   FALSE         24339 2016-11-01        best
# 4:   603       FALSE           NaN   FALSE         24339 2017-05-01        best
# 5:   988       FALSE           NaN   FALSE         24352 2016-05-01        best
# 6:  1331        TRUE     0.6094982   FALSE         24370 2017-05-01        best
# 7:  1513       FALSE           NaN   FALSE         24378 2015-05-01        best
# 8:  1534       FALSE           NaN   FALSE         24378 2016-11-01        best
# 9:  1541       FALSE           NaN   FALSE         24380 2015-05-01        best
# 10:  1555       FALSE           NaN   FALSE         24380 2015-11-01        best
# 11:  1646       FALSE           NaN   FALSE         24383 2017-02-01        best
# 12:  1716       FALSE           NaN   FALSE         24385 2016-05-01        best
# 13:  1779       FALSE           NaN   FALSE         24385 2017-05-01        best
# 14:  1786       FALSE           NaN   FALSE         24386 2015-10-01        best
# 15:  1814       FALSE           NaN   FALSE         24386 2016-11-01        best
# 16:  1898       FALSE           NaN   FALSE         24387 2017-06-01        best
# 17:  2087       FALSE           NaN   FALSE         24393 2015-11-01        best
# 18:  2262       FALSE           NaN   FALSE         26792 2015-05-01        best
# 19:  2416       FALSE           NaN   FALSE         26803 2015-06-01        best
# 20:  2549       FALSE           NaN   FALSE         26808 2016-05-01        best

#### Examine convergence failures
# We know from setup-data-detection.R that there are some unlikely transitions 
# between receiver stations without detection. Do any of these account for the 
# real convergence failures above? Only one!
# For further examination, see refine-patter.R 
iteration[individual_id == 24321 & time_id == as.POSIXct("2017-06-01 00:00:00", tz = "UTC"), ]
iteration[individual_id == 24327 & time_id == as.POSIXct("2015-12-01 00:00:00", tz = "UTC"), ]
iteration[individual_id == 24391 & time_id == as.POSIXct("2015-05-01 00:00:00", tz = "UTC"), ]
iteration[individual_id == 24385 & time_id == as.POSIXct("2016-07-01 00:00:00", tz = "UTC"), ]
iteration[individual_id == 24385 & time_id == as.POSIXct("2017-05-01 00:00:00", tz = "UTC"), ] # included
iteration[individual_id == 24339 & time_id == as.POSIXct("2016-12-01 00:00:00", tz = "UTC"), ]

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

#### Summarise computation time (mins)
# Total computation time
callstats |> 
  group_by(index) |>
  summarise(time = sum(time)) |> 
  ungroup() |> 
  summarise(utils.add::basic_stats(time / 60))
# Summary statisics, by routine
callstats |>
  mutate(routine = if_else(grepl("^filter:", routine), "filter", routine)) |>
  group_by(routine) |>
  reframe(utils.add::basic_stats(time / 60))
# Summary statistics, by routine & sensitivity 
callstats |> 
  group_by(routine, sensitivity) |> 
  reframe(utils.add::basic_stats(time / 60)) |> 
  as.data.table()

#### Visualise total computation time
# > For simulations, total computation time varies from 150 - 190 mins (~3 hours)
tic()
png(here_fig_analysis("computation-time-total.png"), 
    height = 4, width = 5, units = "in", res = 800)
p <- 
  callstats |> 
  group_by(index) |> 
  mutate(time = sum(time) / 60) |> 
  slice(1L) |> 
  as_tibble() |> 
  ggplot(aes("Total", time)) + 
  geom_violin() + 
  xlab("Category") + ylab("Total computation time (mins)")
print(p)
dev.off()

#### Visualise total computation time by sensitivity
# This is simply to check whether any sensitivity runs 
# are associated with much longer computation times
# There are no substantial differences in computation time between runs
tic()
png(here_fig_analysis("computation-time-total-by-sensitivity.png"), 
    height = 4, width = 6, units = "in", res = 800)
p <- 
  callstats |> 
  group_by(index) |> 
  mutate(time = sum(time) / 60) |> 
  slice(1L) |> 
  as_tibble() |> 
  ggplot(aes(sensitivity_label, time, fill = sensitivity_label)) + 
  geom_violin() + 
  xlab("Sensitivity") + ylab("Total computation time (mins)") + 
  labs(fill = "Sensitivity")
print(p)
dev.off()

#### Visualise computation time by routine (~1 s)
tic()
png(here_fig_analysis("computation-time-by-routine.png"), 
    height = 4, width = 5, units = "in", res = 800)
p <- 
  callstats |> 
  as_tibble() |> 
  ggplot(aes(routine_label, time / 60, fill = routine_label)) + 
  geom_violin() + 
  xlab("Routine") + ylab("Computation time (mins)") + 
  guides(fill = "none") +
  theme_bw()
print(p)
dev.off()
toc()

#### As above including the total computation time as a category
png(here_fig_analysis("computation-time-by-routine-with-total.png"), 
    height = 2, width = 6, units = "in", res = 800)
rbind(
  callstats, 
  callstats |> 
    group_by(index) |> 
    mutate(time = sum(time), 
           routine_label = "Total") |> 
    slice(1L) 
  ) |> 
  mutate(routine_label = factor(routine_label, c("Filter: forward", 
                                                 "Filter: backward", 
                                                 "Smoother: two-filter", 
                                                 "Total"))) |> 
  as_tibble() |>
  ggplot(aes(routine_label, time / 60, fill = routine_label)) + 
  geom_violin() + 
  xlab("Routine") + ylab("Computation time (mins)") + 
  guides(fill = "none") +
  theme_bw()
dev.off()

#### Visualise computation time routine & sensitivity (~3 s)
tic()
png(here_fig_analysis("computation-time-by-routine-and-sensitivity.png"), 
    height = 4, width = 12, units = "in", res = 800)
p <- 
  callstats |> 
  ggplot(aes(sensitivity_label, time / 60, fill = sensitivity_label)) + 
  geom_violin() + 
  xlab("Sensitivity") + ylab("Computation time (mins)") + 
  labs(fill = "Sensitivity") + 
  facet_wrap(~routine_label) + 
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
tic()
png(here_fig_analysis("diagnostics-ess.png"), 
    height = 5, width = 10, units = "in", res = 800)
p <-
  diagnostics |> 
  ggplot(aes(ess, fill = routine_label)) + 
  geom_density() + 
  xlab("Effective sample size") + ylab("Relative kernel density") + 
  guides(fill = "none") +
  facet_wrap(~routine_label, nrow = 1, scales = "free") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  scale_y_continuous(expand = c(0, 0), 
                     labels = function(x) prettyGraphics::sci_notation(x, magnitude = 1L)) + 
  theme_bw() + 
  theme(
    plot.margin = margin(t = 10, r = 20, b = 10, l = 20, unit = "pt"), 
    # Hide y axis tick mark labels
    # This is necessary b/c sci_notation is ignored on the final panel
    axis.text.y = element_blank())
print(p)
dev.off()
toc()

#### Summarise area spanned by 95 % of the distribution
# Summary statistics (ncell_core)
diagnostics |> 
  filter(sensitivity == "best") |> 
  filter(routine == "smoother: two-filter") |> 
  summarise(utils.add::basic_stats(ncell_core, na.rm = TRUE))
# Summary statistics (ncell_home)
diagnostics |> 
  filter(sensitivity == "best") |> 
  filter(routine == "smoother: two-filter") |> 
  summarise(utils.add::basic_stats(ncell_home, na.rm = TRUE))
# cf. 31298 grid cells in lake (not NA)
terra::freq(map)
# Visualisation (~14 s)
tic()
png(here_fig_analysis("diagnostics-nell.png"), 
    height = 3, width = 6, units = "in", res = 800)
p <-
  diagnostics |> 
  lazy_dt() |> 
  filter(routine == "smoother: two-filter") |> 
  select(Core = "ncell_core", Home = "ncell_home") |>
  tidyr::pivot_longer(
    cols = c(Core, Home),
    names_to = "type",
    values_to = "ncell"
  ) |>
  as_tibble() |> 
  ggplot(aes(ncell, fill = type)) + 
  geom_density() + 
  xlab("Area (number of cells)") + ylab("Kernel density") + 
  guides(fill = "none") +
  facet_wrap(~type, nrow = 1, scales = "fixed") +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() + 
  theme(plot.margin = margin(t = 10, r = 20, b = 10, l = 20, unit = "pt"))
print(p)
dev.off()
toc()


###########################
###########################
#### Visualise maps

#### Visualise example maps
# TO DO


#### End of code. 
###########################
###########################