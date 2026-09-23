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
#### Select iterations

# For simulations, for comparison with previous work (Futia et al. 2024) and 
# the real-world analysis, we focus on iterations that generated detections
if (analysis == "sim") {
  nrow(iteration)
  iteration <- iteration[n_detections > 0L, ]
  nrow(iteration)
}

# For real-world analysis, focus on julia = TRUE
if (analysis == "real") {
  iteration <- iteration[julia == TRUE, ]
}

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
#### Process callstats 

#### Define iterations
iteration <- iteration[file.exists(file_callstats), ]

#### Read callstats
callstats <- lapply(iteration$index, function(i) {
  iteration$file_callstats[iteration$index == i] |> 
    arrow::read_feather() |> 
    mutate(index = i, .before = 1L) |> 
    cbind(iteration[index == i, .(individual_id, block_id, sensitivity, sensitivity_label)]) |> 
    as.data.table()
}) |> 
  rbindlist() |> 
  mutate(routine_label = stringr::str_to_sentence(routine), .after = routine) |> 
  mutate(routine_label = factor(routine_label, levels = c("Filter: forward", 
                                                          "Filter: backward", 
                                                          "Smoother: two-filter"))) |> 
  as.data.table()

#### Add columns for wahoo-champlain
iteration[, timeline_min := do.call(c, pbapply::pblapply(
  file_timeline,
  \(f) min(arrow::read_feather(f)$timestamp)
))]
iteration[, timeline_max := do.call(c, pbapply::pblapply(
  file_timeline,
  \(f) max(arrow::read_feather(f)$timestamp)
))]
callstats[, timeline_min := iteration$timeline_min[match(index, iteration$index)]]
callstats[, timeline_max := iteration$timeline_max[match(index, iteration$index)]]

#### Record callstats
qs::qsave(callstats, here_output_analysis("synthesis", "callstats.qs"))


###########################
###########################
#### Compute resource requirements

#### Compute total run time (days on 100 cl)
# sim : 0.6895771
# real: 3.300165
sum(callstats$time) / 60 / 60 / 24 / 100

#### Compute total output size by block (MB, GB)
# Compute folder sizes
iteration[, folder_output_block_mb := 
            sapply(seq_len(nrow(iteration)),
                   \(i) dir_size(iteration$folder_output_block[i], recursive = TRUE))]
# Check total size (GB)
# > sim: 12.9 GB
# > real: 56.20986 GB
sum(iteration$folder_output_block_mb) / 1e3

#### As above by chain (MB)
# > real: 22 MB
sum(sapply(unique(iteration$folder_output_chain), dir_size, recursive = TRUE))


###########################
###########################
#### Identify convergence

#### Collate diagnostics,  filtering by block timeline 
# * 1 cl: ~129 s for "real", 88672444 rows, 8.74 GB), 
# * 5 cl: no speed up
tic()
overwrite      <- FALSE
diagnostics.qs <- here_output_analysis("synthesis", "diagnostics.qs")
if (overwrite | !file.exists(diagnostics.qs)) {

  diagnostics <- 
    cl_lapply(iteration$index, .fun = function(i) {
      # Define iteration row for index
      it <- iteration[index == i, ]
      # Define block timeline
      block_timeline <- seq(it$block_start, it$block_end, by = "2 mins")
      # Read diagnostics, filtering by block timeline
      it$file_diagnostics |> 
        arrow::read_feather() |> 
        mutate(index = i, .before = 1L) |> 
        filter(timestamp %in% block_timeline) |> 
        group_by(routine) |> 
        arrange(timestamp, .by_group = TRUE) |> 
        mutate(timestep = 1:n()) |> 
        ungroup() |> 
        arrange(routine, timestamp) |>
        cbind(it[index == i, .(individual_id, block_id, sensitivity)]) |> 
        as.data.table()
    }) |> 
    rbindlist()  |> 
    lazy_dt(immutable = FALSE) |> 
    mutate(routine_label = stringr::str_to_sentence(routine), .after = routine) |> 
    mutate(routine_label = factor(routine_label, levels = c("Filter: forward", 
                                                            "Filter: backward", 
                                                            "Smoother: two-filter"))) |> 
    arrange(index, routine, timestamp) |> 
    as.data.table()
  
  qs::qsave(diagnostics, diagnostics.qs)

} else {
  diagnostics <- qs::qread(diagnostics.qs)
}
# lobstr::obj_size(diagnostics)
toc()

#### Define convergence for an example individual
diagnostics |> 
  filter(index == 8L) |> 
  filter(routine == "smoother: two-filter") |> 
  summarise(prop = length(which(!is.na(ess))) / n())

#### Define convergence 
# Defile pass_filter_fwd, pass_filter_bwd
convergence_filter <- 
  callstats |>
  select(index, routine, convergence) |>
  filter(routine %in% c("filter: forward", "filter: backward")) |>
  mutate(routine = recode(routine,
                          "filter: forward"  = "pass_filter_fwd",
                          "filter: backward" = "pass_filter_bwd")) |>
  tidyr::pivot_wider(names_from  = routine,
                     values_from = convergence) |>
  mutate(across(starts_with("pass_filter_"), ~ tidyr::replace_na(.x, FALSE))) |>
  as.data.table()
# Collect pass_filter_bwd, pass_filter_bwd and pass_smoother
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
  ungroup() |> 
  mutate(success = pass_filter & (pass_smoother >= 0.75)) |> 
  left_join(convergence_filter, by = "index") |> 
  select("index", "pass_filter_fwd", "pass_filter_bwd", "pass_filter", "pass_smoother", "success") |> 
  left_join(iteration[, .(index, individual_id, chain_id, block_id, sensitivity)], by = "index") |> 
  # Define success_chain = TRUE if all blocks for that chain succeeded
  group_by(individual_id, sensitivity, chain_id) |> 
  mutate(success_chain = all(success)) |> 
  ungroup() |> 
  as.data.table()

#### Check convergence by block, with comments on the "real" analysis
# Review convergence (sense check)
callstats[convergence == FALSE, ]
convergence[success == FALSE, ]
# Review filter issues (82 % pass rate)
# * 1243 / (267 +  1243) pass filter
table(convergence$pass_filter)
# For the filter, increasing the number of particles helps
table(callstats$convergence, callstats$n_particle)
# Review smoothing (85 % pass rate for successful filter runs)
# * 1067 / (1067 + 176)
utils.add::basic_stats(convergence$pass_smoother, na.rm = TRUE)
utils.add::basic_stats(convergence$pass_smoother[convergence$success], na.rm = TRUE)
table(convergence$pass_smoother[convergence$pass_filter == TRUE] > 0.75)
hist(convergence$pass_smoother)
# Overall convergence rate is 71 %
# * 1067 /( 1067 + 443)
table(convergence$success)
table(convergence$success, convergence$sensitivity == "best")
table(convergence$success, convergence$sensitivity != "best")
table(convergence$sensitivity[convergence$success == FALSE])

#### Check convergence by chain
# For the real analysis, we have 1510 individual/month-year blocks
# We have 512 individual/season-year blocks
# We have complete convergence for only 242/512 blocks (47 % success)
convergence_chains <- 
  convergence |> 
  group_by(individual_id, sensitivity, chain_id) |>
  summarise(success_chain = all(success)) |> 
  ungroup() |>
  as.data.table()
convergence_chains |> 
  count(success_chain)

#### Review convergence failures
convergence[success == FALSE, ]

#### Examine convergence failures
# We know from setup-data-detection.R that there are some unlikely transitions 
# between receiver stations without detection. Do any of these account for the 
# real convergence failures above? Only one!
# For further examination, see refine-patter.R 
if (analysis == "real") {
  iteration[individual_id == 24321 & block_start == as.POSIXct("2017-06-01 00:00:00", tz = "UTC"), ]
  iteration[individual_id == 24327 & block_start == as.POSIXct("2015-12-01 00:00:00", tz = "UTC"), ]
  iteration[individual_id == 24391 & block_start == as.POSIXct("2015-05-01 00:00:00", tz = "UTC"), ]
  iteration[individual_id == 24385 & block_start == as.POSIXct("2016-07-01 00:00:00", tz = "UTC"), ]
  iteration[individual_id == 24385 & block_start == as.POSIXct("2017-05-01 00:00:00", tz = "UTC"), ] # included
  iteration[individual_id == 24339 & block_start == as.POSIXct("2016-12-01 00:00:00", tz = "UTC"), ]
}

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

#### Update callstats$convergence
# Set callstats$convergence = TRUE to avoid confusion
# (callstats$convergence may be FALSE for some indicies b/c 
# for smoothing runs b/c/ convergence there is defined by Patter.jl 
# not using the criteria in this project)
callstats[, convergence := TRUE]


###########################
###########################
#### Summarise callstats

#### Summarise computation time (mins), by block
# Total computation time
# (Note that for the real-world analysis, the number of particles varies)
callstats |> 
  group_by(index) |>
  summarise(time = sum(time)) |> 
  ungroup() |> 
  summarise(utils.add::basic_stats(time / 60))
# Summary statistics, by routine
callstats |>
  mutate(routine = if_else(grepl("^filter:", routine), "filter", routine)) |>
  group_by(routine) |>
  reframe(utils.add::basic_stats(time / 60))
# Summary statistics, by routine & sensitivity 
callstats |> 
  group_by(routine, sensitivity) |> 
  reframe(utils.add::basic_stats(time / 60)) |> 
  as.data.table()
# Cf. the number of time steps for successful algorithm runs (~6 s)
# 20160 -> 58572 time steps
iteration[, nt := pbapply::pbsapply(iteration$file_timeline, 
                                    \(f) arrow::read_feather(f) |> nrow())]
callstats[, nt := iteration$nt[match(index, iteration$index)]]
utils.add::basic_stats(callstats$nt)

#### Summarise computation time, by time step and particle
#
# A) Total computation time per time step (secs)
callstats |> 
  group_by(index) |> 
  mutate(time = sum(time)) |> 
  slice(1L) |> 
  ungroup() |> 
  summarise(time_per_t = time / nt) |> 
  reframe(utils.add::basic_stats(time_per_t, p = NULL))
# min  mean median   max    sd   IQR   MAD
# <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>
# 0.918  1.36   1.31  2.20 0.241 0.385 0.251
#
# B) Computation time per time step by routine (secs)
callstats |> 
  mutate(routine = if_else(grepl("^filter:", routine), "filter", routine)) |>
  mutate(time_per_t = time / nt) |> 
  group_by(routine) |> 
  reframe(utils.add::basic_stats(time_per_t, p = NULL))
# routine                 min  mean median   max     sd    IQR    MAD
# <chr>                 <dbl> <dbl>  <dbl> <dbl>  <dbl>  <dbl>  <dbl>
# 1 filter               0.0732 0.201  0.193 0.543 0.0475 0.0620 0.0450
# 2 smoother: two-filter 0.656  0.961  0.917 1.57  0.171  0.281  0.181 
#
# C) Computation time per time step per particle by routine
callstats |> 
  mutate(routine = if_else(grepl("^filter:", routine), "filter", routine)) |>
  mutate(time_per_t_per_particle = time / nt / n_particle) |> 
  group_by(routine) |> 
  reframe(utils.add::basic_stats(time_per_t_per_particle, p = NULL))
# routine                     min       mean     median        max          sd        IQR         MAD
# <chr>                     <dbl>      <dbl>      <dbl>      <dbl>       <dbl>      <dbl>       <dbl>
# 1 filter               0.00000146 0.00000397 0.00000385 0.00000717 0.000000849 0.00000122 0.000000883
# 2 smoother: two-filter 0.000262   0.000385   0.000367   0.000627   0.0000683   0.000113   0.0000722 

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
#### Summarise diagnostics, by block

#### Summarise ESS (for filters & smoother)
# Summary statistics 
diagnostics |> 
  group_by(routine) |> 
  reframe(utils.add::basic_stats(ess, na.rm = TRUE)) 
# As above but grouping both filter runs
diagnostics[, routine_simple := 
              if_else(routine %in% c("filter: forward", "filter: backward"),
                      "filter", "smoother")]
diagnostics |> 
  group_by(routine_simple) |> 
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
# cf. 31296 grid cells in lake (not NA)
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
#### Record datasets

qs::qsave(convergence, here_output_analysis("synthesis", "convergence.qs"))
qs::qsave(convergence_chains, here_output_analysis("synthesis", "convergence-chains.qs"))


#### End of code. 
###########################
###########################