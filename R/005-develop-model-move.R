###########################
###########################
#### develop-model-move.R

#### Aims
# 1) Develop a movement model for lake trout

#### Prerequisites
# We consider the following data sources:
# 1) Detection data collected by M. Futia. 
# 2) Movement data provided by Blanchfield et al. (2023)
# 3) Calibration equations developed by Reeve et al. (2024)


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Set global options
Sys.setenv("JULIA_SESSION" = FALSE)
set.seed(1L)

#### Load essential packages
library(proj.verse)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(tictoc)
library(truncdist)
files_source_r(here_src())

#### Load data
map         <- terra::rast(here_input("map.tif"))
fish        <- qs::qread(here_input_real("fish.qs"))
moorings    <- qs::qread(here_input("moorings.qs"))
detections  <- qs::qread(here_input_real("detections.qs"))
blanchfield <- qs::qread(here_data("supp", "model-move", "blanchfield.qs"))


###########################
###########################
#### Explore shapes
# This code explores candidate shapes for a distribution of step lengths

#### Speeds
# 0.1 m/s -> 12 m/2 min,    18 m/3 min
# 0.2 m/s -> 24 m/2 min,    36 m/3 min
# 0.3 m/s -> 36 m/2 min,    54 m/3 min
# 0.4 m/s -> 48 m/2 min,    72 m/3 min
# 0.9 m/ms -> 108 m/2 min,  180 m/3 min

#### Define mobility 
mobility <- 200

#### Normal 
# This distribution probably permits overly low step lengths 
# But otherwise covers a broad range of possible cruising speeds
curve(dtrunc(x, "norm", a = 0, b = mobility, 24, 20), from = 0, to = mobility)

#### Gamma
# Reduce shape to shift to left (scale-dependent)
# Reduce rate to widen distribution
curve(dtrunc(x, "gamma", a = 0, b = mobility, 5, 0.25), from = 0, to = mobility)
curve(dtrunc(x, "gamma", a = 0, b = mobility, 3, 0.15), from = 0, to = mobility)

#### Gamma (ggplot2)
# Define example parameters
p1 <- c(2.8, 1 / 0.05)
# Shrink/expand distribution e.g., by 120 % while maintaining the same mode 
p2 <- gamma_rescale(p1[1], p1[2], fact = 1.2)
p3 <- gamma_rescale(p1[1], p1[2], fact = 0.8)
# Visualise Gamma distributions
ggplot(data.frame(x = c(0, mobility)), aes(x = x)) +
  # stat_dtruncgamma(mobility, shape = 4.5, scale = 1/0.1) + 
  # stat_dtruncgamma(mobility, shape = 3.2, scale = 1/0.06, col = "blue") + 
  stat_dtruncgamma(mobility, shape = 4.5, scale = 1/0.10, 
                   col = "grey", linetype = 2, linewidth = 1.5) + 
  stat_dtruncgamma(mobility, shape = p1[1], scale = p1[2], col = "black") + 
  stat_dtruncgamma(mobility, shape = p2[1], scale = p2[2], col = "green") + 
  stat_dtruncgamma(mobility, shape = p3[1], scale = p3[2], col = "red") + 
  scale_x_continuous(breaks = seq(0, mobility, by = 10))  + 
  theme_bw()

#### Log normal
# Reduce sdlog to broaden distribution
# This parameterisation peaks too early
curve(dtrunc(x, "lnorm", a = 0, b = mobility, 3, 1), from = 0, to = mobility)

#### Cauchy
# This distribution also permits overly low step lengths
# But has a longer tail to the right
# Increase scale to widen distribution 
curve(dtrunc(x, "cauchy", a = 0, b = mobility, 20, 10), from = 0, to = mobility)


###########################
###########################
#### Accelerometry analyses of movement speeds

#### Define accelerometry data.table
if (FALSE) {
  
  ## (A) Simulate example acceleration measurements (m/s^2)
  n_id    <- 50L     # number of individuals
  n_obs   <- 10000L  # number of observations per individual
  n_id * n_obs       # total number of observations (~500,000)
  # Simulate A measurements
  # * Use an approximate R
  blanchfield <- CJ(id = 1:n_id, index = 1:n_obs)
  blanchfield[, A := runif(nrow(blanchfield), 0, 2.5)]
  blanchfield[, log10_A := log10(A)]
  blanchfield[, A := NULL]
  blanchfield[, index := NULL]
  head(blanchfield)
  
  # id     log10_A
  # <int>   <num>
  #  1    -0.17798130
  #  1    -0.03137243
  
} else {
  
  ## (B) Use real-world measurements
  blanchfield <- 
    blanchfield |> 
    mutate(log10_A = log10(accel)) |> 
    select(id = individual_id, timestamp, log10_A) |> 
    filter(!is.infinite(log10_A)) |>
    as.data.table()
  
}

#### Define acceleration ~ speed model estimates (Reeve et al., 2024)
# Fixed effects
alpha    <- 0.15353 # intercept
beta     <- 0.36788 # gradient (effect of mean acceleration)
# Uncertainties
alpha_se <- 0.01529 # alpha
beta_se  <- 0.02338 # beta
corr     <- 0.224   # alpha-beta correlation
sigma_id <- 0.05566 # random variation among individuals
sigma    <- 0.04809 # random noise

#### Simulate speeds (BL/s)
# Simulate speeds ignoring all uncertainty
cl  <- 6L
SS1 <- sim_SS(blanchfield,
              alpha    = alpha,
              alpha_se = 0,
              beta     = beta,
              beta_se  = 0,
              corr     = NA,
              sigma_id = 0,
              sigma    = 0,
              n_sim    = 1L)
# Simulate speeds ignoring individual variation
SS2 <- sim_SS(blanchfield,
              alpha    = alpha,
              alpha_se = alpha_se,
              beta     = beta,
              beta_se  = beta_se,
              corr     = corr,
              sigma_id = 0,
              sigma    = sigma,
              cl       = cl)
# Simulate speeds ignoring coefficient uncertainty
SS3 <- sim_SS(blanchfield,
              alpha    = alpha,
              alpha_se = 0,
              beta     = beta,
              beta_se  = 0,
              corr     = corr,
              sigma_id = sigma_id,
              sigma    = sigma,
              cl       = cl)
# Simulate speeds, accounting for uncertainties
SS4 <- sim_SS(blanchfield,
              alpha    = alpha,
              alpha_se = alpha_se,
              beta     = beta,
              beta_se  = beta_se,
              corr     = corr,
              sigma_id = sigma_id,
              sigma    = sigma,
              cl       = cl)

#### Repeat simulation for SS4 with thinning strategies
# This is used to evaluate robustness under autocorrelation
SS5 <- cl_lapply(c("1 hour", "2 hours", "12 hours", "full"), function(thin) {
  # (optional) Thin time series 
  if (thin != "full") {
    blanchfield_thin <- 
      blanchfield |> 
      mutate(timestamp = lubridate::round_date(timestamp, "days")) |> 
      group_by(id, timestamp) |> 
      slice(1L) |>
      as.data.table() 
  } else {
    blanchfield_thin <- copy(blanchfield)
  }
  # Simulate swimm speeds
  SS <- sim_SS(blanchfield_thin,
               alpha    = alpha,
               alpha_se = alpha_se,
               beta     = beta,
               beta_se  = beta_se,
               corr     = corr,
               sigma_id = sigma_id,
               sigma    = sigma,
               cl       = cl)
    data.table(grp = thin, SS = as.vector(SS))
}) |> rbindlist()
# Summarise the number of observations per group
SS5 |> 
  group_by(grp) |> 
  summarise(n())

#### Repeat simulation with random subsampling (50 %)
SS6 <- lapply(c(0.25, 0.5, 0.75), function(prop) {
  SS <- sim_SS(blanchfield[sample.int(floor(prop * .N)), ],
               alpha    = alpha,
               alpha_se = alpha_se,
               beta     = beta,
               beta_se  = beta_se,
               corr     = corr,
               sigma_id = sigma_id,
               sigma    = sigma,
               cl       = cl)
  data.table(grp = prop, SS = as.vector(SS))
}) |> rbindlist()

#### Visualise speeds (BL/s)
# (A) Visualise speeds, accounting for different sources of uncertainty 
# > Accounting for all sources of uncertainty makes a minimal difference 
xmax <- ceiling(max(c(SS1, SS2, SS3, SS4)))
xlim <- c(0, xmax)
plot(density(SS1), xlim = xlim, col = "red")  # SS1 : ignore all
lines(density(SS2), col = "orange")           # SS2: ignore ID
lines(density(SS3), col = "green")            # SS3: ignore coef
lines(density(SS4), col = "black")            # SS4: full
# (B) Visualise speeds, accounting for autocorrelation
# > The distributions are highly robust to autocorrelation
plot(density(SS4), xlim = xlim)
cols <- rainbow(SS5[, uniqueN(grp)])
SS5[, lines(density(SS), col = cols[.GRP]), by = grp]
legend("topright", legend = SS5[, unique(grp)], lty = 1, col = cols)
# (C) Visualise speeds, for random subsets
# > The results are robust to the data structure (individuals, time steps)
plot(density(SS4), xlim = xlim)
cols <- rainbow(SS6[, uniqueN(grp)])
SS6[, lines(density(SS), col = cols[.GRP]), by = grp]
legend("topright", legend = SS6[, unique(grp)], lty = 1, col = cols)

#### Compute speeds (m/s or m/ 120 s)
# We have generated a distribution of swimming speeds (BL/s)
# We use body sizes to translate this into a distribution in m/s
# TO DO
# * Review whether we need to use fork length 
s <- 120 
xlim <- c(0, max(SS4 * max(fish$len)) * s)
dmin <- density(SS4 * min(fish$len) * s)
dmax <- density(SS4 * max(fish$len) * s)
plot(dmin)
lines(dmax, col = "dimgrey")

#### Model speeds
# Plot 'observed' distributions
plot(dmin, lwd = 0.5)
lines(dmax, col = "dimgrey", lwd = 0.5)
# Plot initial 'best-guess' model
# plot_dbn("gamma", xlim = xlim, add = TRUE, pars = list(shape = 2.33, scale = 27), col = "blue")
# Add data-driven models 
# plot_dbn("norm", xlim = xlim, add = TRUE, pars = list(mean = 60, sd = 50))
# plot_dbn("cauchy", xlim = xlim, add = TRUE, pars = list(location = 60, scale = 50), col = "red")
plot_dbn("gamma", xlim = xlim, add = TRUE, pars = list(shape = 3, scale = 30), col = "blue")

#### Use quantiles to inform mobility
quantile(SS4, 0.99)


###########################
###########################
#### Detection analyses of mobility 

#### Estimation method
# Using detection data, we can estimate mobility
# We compute three speed metrics & for each we compute the min/mean/max value
# Metrics:
# * speed_min: minimum speed, based on movement between nearest container edges 
# * speed_avg: average speed, based on movement between receivers
# * speed_max: maximum speed, based on movement between farthest container edge
# Interpretation:
# * min speeds are 'inappropriately low' b/c indirect routes between receivers are taken
# * maximum values for each metric can inform us about mobility 

#### Estimate mobility (~4 s)
if (requireNamespace("flapper", quietly = TRUE)) {
  tic()
  # Define LCP grid (resolution in x and y must be exactly identical)
  bb <- terra::ext(map)[1:4]
  bb <- plyr::round_any(bb, 100)
  grid <- terra::rast(terra::ext(bb), res = 100)
  grid <- terra::resample(map, grid)
  # Prepare observations
  detections[, fct := individual_id]
  msp <- sp::SpatialPoints(moorings[, c("receiver_x", "receiver_y")], sp::CRS(terra::crs(map)))
  msp <- sp::SpatialPointsDataFrame(msp, data.frame(receiver_id = moorings$receiver_id))
  # Compute speeds
  # TO DO 
  # * Repeat with best-guess detection range
  mvt <- flapper::get_mvt_mobility_from_acoustics(data = detections, 
                                                  fct = "individual_id", 
                                                  moorings = msp, 
                                                  detection_range = 6000, 
                                                  calc_distance = "lcp", 
                                                  bathy = raster::raster(grid),
                                                  step = 120,
                                                  transmission_interval = 160)
  toc()
}

#### Results 
# These results are broadly robust to detection_range, step and transmission_interval
# --------------------------------------
#   Estimates (m/s)-----------------------
#   variable min mean  max
# 1 speed_min_ms   0 0.03 0.40
# 2 speed_avg_ms   0 0.10 0.90
# 3 speed_max_ms   0 0.18 1.64
# --------------------------------------
#   Estimates (m/step)--------------------
#   variable  min  mean    max
# 1 speed_min_mstep 0.00  3.92  47.71
# 2 speed_avg_mstep 0.07 12.52 108.29
# 3 speed_max_mstep 0.12 21.12 197.28
# --------------------------------------

#### Conclusions 
# Even if we assume 'minimum distance' movements, maximum speeds of 48 m/2 min are apparent.
# (But travelled  speeds are likely to be higher than this value 
# ... as this assumes the minimum possible travel distance)
# If we assume larger movements, speeds may be up to 109 m/2 min.
# This is a reasonable middle-of-the-road estimate for maximum movement speeds. 
# In the most extreme case (unlikely), movement speeds up to 198 m/2 min are apparent.
# (This is the upper bound for movement speeds suggested by the data.)
# These results are consistent for different transmission intervals
# (But the movement model needs to consider the effect of random transmission)


###########################
###########################
#### Analyse turning angles

# TO DO
# M. Futia to add code. 


###########################
###########################
#### Synthesise parameters

#### Define parameters
# TO DO
pars_model_move <- list(shape = 3.0, scale = 30.0, mobility = 350, phi = 1.2)
qs::qsave(pars_model_move, here_input("pars-model-move.qs"))

#### Build validity map(s)
mobility <- pars_model_move$mobility
vmap <- patter:::spatVmap(.map = map, .mobility = mobility, .plot = TRUE)
terra::writeRaster(vmap, here_input_real("vmap.tif"), overwrite = TRUE)


#### End of code. 
###########################
###########################