###########################
###########################
#### develop-model-move-blanchfield.R

#### Aims
# 1) Develop a movement model for lake trout with analyses of Blanchfield et al. (2023) data

#### Prerequisites
# 1) Process Blanchfield et al. (2023) data


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
files_source_r(here_src())

#### Load data
fish        <- qs::qread(here_input("fish.qs"))
blanchfield <- qs::qread(here_data("supp", "model-move", "blanchfield.qs"))


###########################
###########################
#### Simulate swimming speeds

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
cl  <- 10L
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

#### Record simulation outputs
qs::qsave(SS4, here_data("supp", "model-obs", "SS4.qs"))


###########################
###########################
#### Visualise speeds (BL/s)

#### (A) Visualise speeds, accounting for different sources of uncertainty 
# Accounting for all sources of uncertainty makes a minimal difference 
xmax <- ceiling(max(c(SS1, SS2, SS3, SS4)))
xlim <- c(0, xmax)
plot(density(SS1), xlim = xlim, col = "red")  # SS1 : ignore all
lines(density(SS2), col = "orange")           # SS2: ignore ID
lines(density(SS3), col = "green")            # SS3: ignore coef
lines(density(SS4), col = "black")            # SS4: full

#### (B) Visualise speeds, accounting for autocorrelation
# The distributions are highly robust to autocorrelation
plot(density(SS4), xlim = xlim)
cols <- rainbow(SS5[, uniqueN(grp)])
SS5[, lines(density(SS), col = cols[.GRP]), by = grp]
legend("topright", legend = SS5[, unique(grp)], lty = 1, col = cols)

#### (C) Visualise speeds, for random subsets
# The results are robust to the data structure (individuals, time steps)
plot(density(SS4), xlim = xlim)
cols <- rainbow(SS6[, uniqueN(grp)])
SS6[, lines(density(SS), col = cols[.GRP]), by = grp]
legend("topright", legend = SS6[, unique(grp)], lty = 1, col = cols)

#### Compute speeds (m/s or m/ 120 s)
# We have generated a distribution of swimming speeds (BL/s)
# We use body sizes to translate this into a distribution in m/s
s <- 1 
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
plot_dbn("gamma", xlim = xlim, add = TRUE, pars = list(shape = 3.25, scale = 25), col = "blue")

#### Summary statistics
# Summary statistics (BL/s)
mean(SS4)
sd(SS4)
quantile(SS4, 0.99)
# Summary statistics (m/s or m per two min) for small fish 
s <- 120
mean(SS4 * min(fish$len)) * s           # 0.4671949, 56.06339 (m per 2 min)
sd(SS4 * min(fish$len)) * s             # 0.1740966, 20.8916
quantile(SS4 * min(fish$len) * s, 0.99) # 1.078395,  129.4074
# Summary statistics (m/s or m per two min) for BIG fish 
mean(SS4 * max(fish$len)) * s           # 0.7029355, 84.35226 (m per 2 min)
sd(SS4 * max(fish$len)) * s             # 0.2619436, 31.43323
quantile(SS4 * max(fish$len) * s, 0.99) # 1.622539,  194.7047


#### End of code. 
###########################
###########################