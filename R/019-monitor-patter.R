###########################
###########################
#### monitor-patter.R

#### Aims
# 1) This script monitors particle algorithm progress

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
here_input_analysis  <- switch_here_input_analysis_subanalysis(analysis, subanalysis)
here_output_analysis <- switch_here_output_analysis_subanalysis(analysis, subanalysis)
iteration            <- qs::qread(here_input_analysis("iteration.qs"))


###########################
###########################
#### Monitor progress

#### Monitor progress
# table(file.exists(iteration$file_callstats_filter))
table(file.exists(iteration$file_callstats))
table(file.exists(iteration$file_diagnostics))

#### Monitor number of Julia processes running
# Display for user user "lavended" up to N = 30 processes
if (Sys.info()[["nodename"]] == "siam-linux20") {
  system("top -b -n 1 | awk 'NR==1 || $2==\"lavended\"' | head -n 30")
}

#### Check rows with missing files
# iteration[!file.exists(file_callstats_filter), .(index, unit_id, individual_id, block_id, sensitivity)]
iteration[!file.exists(file_callstats), .(index, unit_id, individual_id, block_id, sensitivity)]

#### Check for errors
# Read logfiles 
logfiles <- list.files(here_output_analysis("logs"), full.names = TRUE, pattern = "\\.log$")
logfiles <- gtools::mixedsort(logfiles, decreasing = TRUE)
rows <- 
  logfiles |>
  basename() |>
  stringr::str_extract("\\d+(?=\\.log$)") |> 
  as.integer()
logs <- lapply(logfiles, readLines)
names(logs) <- as.character(rows)
# Check for missing logs
table(iteration$index %in% rows)
# Read logs
logs
# Search for 'error', 'fail' or simular
logtxt <- do.call(paste, lapply(logs, function(l) paste0(l, collapse = ", ")))
stringr::str_detect(tolower(logtxt), "error")
stringr::str_detect(tolower(logtxt), "fail")
stringr::str_detect(tolower(logtxt), "failure")


#### End of code. 
###########################
###########################