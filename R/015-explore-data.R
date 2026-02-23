###########################
###########################
#### explore-data.R

#### Aims
# 1) Explore detection datasets

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
regions    <- terra::rast(here_input("regions.tif"))
regions_cs <- qs::qread(here_input("regions-colour-scheme.qs"))
fish       <- qs::qread(here_input("fish.qs"))


###########################
###########################
#### Select analysis

#### Define analysis 
# analysis <- "sim"
analysis <- "real"
subanalysis <- "main"

#### Define analysis-specific data
here_input_analysis <- switch_here_input_analysis_subanalysis(analysis, subanalysis)
iteration           <- qs::qread(here_input_analysis("iteration.qs"))
detections          <- qs::qread(here_input_analysis("detections.qs"))
moorings            <- qs::qread(here_input_analysis("moorings.qs"))


###########################
###########################
#### Summarise detection dataset used for modelling 
# cf. raw data summary statistics (setup-data-detection.R)

# Number of observations
nrow(detections)
# Number of individuals
length(unique(detections$individual_id))
# Time ranges (in months & years)
range(detections$timestamp)
int <- lubridate::interval(min(detections$timestamp),max(detections$timestamp))
lubridate::time_length(int, "months")
lubridate::time_length(int, "years")
# Number of receivers with detections
length(unique(detections$receiver_id))
# Number of time_id blocks with detections
detections[, time_id := lubridate::floor_date(timestamp, "months")]
length(unique(paste(detections$individual_id, detections$time_id)))


###########################
###########################
#### Abacus plot 

#### Visualise real detection dataset used for modelling 
# Plot raw time series (light grey)
# Add modelled time series, coloured by region as in map

overwrite      <- TRUE
detections.png <- here_fig(analysis, "detections.png")

if (analysis == "real" & (overwrite | !file.exists(detections.png))) {
  
  #### Define mooring regions
  moorings <- 
    moorings |> 
    mutate(region = terra::extract(regions, cbind(receiver_x, receiver_y))$map_value, 
           col = regions_cs$col[match(region, regions_cs$region)]) |> 
    as.data.table()
  stopifnot(all(!is.na(moorings$col)))
  
  #### Define raw detection time series
  detections_raw <- 
    here_input_analysis("detections-raw.qs") |> 
    qs::qread() |> 
    select("individual_id", "timestamp", "receiver_id") |> 
    as.data.table()
  
  #### Define modelled detection time series (~20 s)
  # Isolate detections within period of interest
  # Note that additional detections may have been included outwidth the modelled period
  # (to support model initialisation)
  detections_mod <- detections[timestamp >= as.POSIXct("2014-12-01 00:00:00"), ]
  
  #### Collate detections (raw, modelled)
  head(detections_raw)
  head(detections_mod)
  ids  <- unique(detections_raw$individual_id)
  draw <- copy(detections_raw)
  draw[, model := FALSE]
  draw[detections_mod, on = .(individual_id, timestamp), model := TRUE]
  draw[, individual_id := factor(individual_id, levels = ids)]
  draw[, col := moorings$col[match(receiver_id, moorings$receiver_id)]]
  draw[, col := ifelse(model == TRUE, col, scales::alpha("dimgrey", 0.25))]
  
  #### Plot (~15 s)
  tic()
  range(draw$timestamp)
  png(detections.png, 
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
#### Review data availability 

if (analysis == "real") {
  
  # Review data availability by spawning site/season
  # (Check the mean number of days with observations per individual for each site/season)
  detections |> 
    mutate(site = fish$site[match(individual_id, fish$individual_id)], 
           season = season_factor(time_id)) |> 
    group_by(individual_id, time_id, site, season) |> 
    summarise(n = length(unique(lubridate::floor_date(timestamp, "days")))) |> 
    ungroup() |> 
    group_by(site, season) |> 
    summarise(utils.add::basic_stats(n)) |> 
    as.data.table()
  
  # Count the number of individual/month time series per site/season
  # (There may be times of year e.g., summer when data are poorer)
  categories <- 
    rbind(
      # Number of time series in full dataset
      detections |> 
        mutate(site = fish$site[match(individual_id, fish$individual_id)], 
               season = season_factor(time_id)) |> 
        group_by(site, season) |> 
        summarise(n = n_distinct(paste(individual_id, time_id))) |> 
        mutate(dataset = "full") |> 
        as.data.table(),
      # Number of time series in modelled dataset
      detections |> 
        mutate(site = fish$site[match(individual_id, fish$individual_id)], 
               season = season_factor(time_id)) |> 
        group_by(site, season) |> 
        summarise(n = n_distinct(paste(individual_id, time_id))) |> 
        mutate(dataset = "model") |> 
        as.data.table()
    )
  
  # Visualise number of time series available versus modelled
  p <- 
    ggplot(categories, aes(x = dataset, y = n)) +
    geom_bar(stat = "identity") +
    ylab("Number of time series") + 
    facet_grid(~site ~ season)
  plotly::ggplotly(p)
  
  # Visualise % of time series modelled 
  p <- 
    categories |>
    tidyr::pivot_wider(names_from = dataset,
                       values_from = n) |>
    mutate(percent = 100 * model / full) |> 
    as_tibble() |> 
    ggplot(aes(x = season, y = percent)) +
    geom_col() +
    ylab("Percentage of time series modelled") + 
    facet_wrap(~site)
  plotly::ggplotly(p)
  
}


#### End of code.
###########################
###########################