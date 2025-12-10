###########################
###########################
#### develop-model-move.R

#### Aims
# 1) Develop a movement model for lake trout

#### Prerequisites
# We consider the following data sources:
# 1) Movement data provided by Blanchfield et al. (2023)
#    & Calibration equations developed by Reeve et al. (2024)
# 2) VPS data collected by M. Futia.
# 3) Detection data collected by M. Futia. 


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
library(truncdist)
files_source_r(here_src())

#### Load data
map       <- terra::rast(here_input("map.tif"))
fish      <- qs::qread(here_input("fish.qs"))
SS4       <- qs::qread(here_data("supp", "model-obs", "SS4.qs"))
vps_step  <- qs::qread(here_data("supp", "model-move", "futia-step.qs"))
vps_angle <- qs::qread(here_data("supp", "model-move", "futia-angle.qs"))


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
#### Synthesise parameters

#### Define 'best-guess' parameters (list)
pars_model_move_best <- list(shape = 3.25, scale = 25.0, mobility = 216, phi = 0.4)

#### Define restrictive/flexible parameters 
# Define parameter uncertainty 
adj      <- 0.25
inflate  <- 1 + adj
deflate  <- 1 - adj
pars_adj <- list(inflate = inflate, deflate = deflate)
# Collect 'best-guess' parameters
mobility <- pars_model_move_best$mobility
shape    <- pars_model_move_best$shape
scale    <- pars_model_move_best$scale
phi      <- pars_model_move_best$phi
# Define more restrictive/flexible models
restrictive <- gamma_rescale(shape, scale, fact = deflate)
flexible    <- gamma_rescale(shape, scale, fact = inflate)

#### Collect movement parameters in data.table (best, restrictive, flexible)
pars_model_move_full <- data.table(mobility = c(mobility, mobility * deflate, mobility * inflate),
                                   shape = c(shape, restrictive[1], flexible[1]),
                                   scale = c(scale,  restrictive[2], flexible[2]),
                                   phi = c(phi, phi * deflate, phi * inflate))


###########################
###########################
#### Publication-quality visualisation of movement model

#### Isolate VPS datasets
# We compute density below 
drummond_step <- vps_step[site == "Drummond", ]
thunder_step  <- vps_step[site == "Thunder", ]

#### For accelerometry, precompute observed step length densities  (slow)
# Define speeds (values) for 120 s
s    <- 120
vmin <- SS4 * min(fish$len) * s
vmax <- SS4 * max(fish$len) * s
# Round values to add as a rug (round for speed)
# vrug <- unique(c(plyr::round_any(vmin, 5), plyr::round_any(vmax, 5)))
# Compute densities
dmin <- density(vmin, from = 0)
dmax <- density(vmax, from = 0)

#### Step lengths
png(here_fig("model-move-step.png"), 
    height = 4, width = 4, units = "in", res = 800)
pp <- par(mgp = c(3, 0.7, 0))
# Set graphical parameters
xlim <- c(0, 500)
# ylim <- c(0, 0.025)
ylim <- c(0, 0.06)
plot(dmin, 
     type = "n",
     xlim = xlim, ylim = ylim,
     xlab = "", ylab = "", main = "",
     axes = FALSE)
# Observed distributions for VPS analyses 
add_poly(density(drummond_step$step_length, from = 0), col = scales::alpha("lightblue", 1))
add_poly(density(thunder_step$step_length, from = 0), col = scales::alpha("blue", 0.5))
# Observed distributions for small and large fish 
add_poly(dmin, col = scales::alpha("lightgrey", 0.75))
add_poly(dmax, col = scales::alpha("dimgrey", 0.5))
# Best model (step-length)
plot_dbn("gamma", 
         xlim = c(0, 300), 
         pars = list(shape = pars_model_move_full$shape[1],
                     scale = pars_model_move_full$scale[1]), 
         upper = pars_model_move_full$mobility[1],
         add = TRUE, lwd = 2)
# Restrictive model (step-length)
plot_dbn("gamma", 
         xlim = c(0, 300), 
         pars = list(shape = pars_model_move_full$shape[2], 
                     scale = pars_model_move_full$scale[2]), 
         upper = pars_model_move_full$mobility[2],
         add = TRUE, col = "red", lty = 1, lwd = 1)
# Flexible model (step-length)
plot_dbn("gamma", 
         xlim = c(0, 300), 
         pars = list(shape = pars_model_move_full$shape[3], 
                     scale = pars_model_move_full$scale[3]), 
         upper = pars_model_move_full$mobility[3],
         add = TRUE, col = "darkgreen", lty = 1, lwd = 1)
# Mark mobility
mark_mobility(pars_model_move_full$mobility[1])
mark_mobility(pars_model_move_full$mobility[2], col = "red")
mark_mobility(pars_model_move_full$mobility[3], col = "darkgreen")
# Add axes (m/s, m per two min, density)
axis(side = 1, c(xlim[1], xlim[2]), labels = c("", ""), lwd.tick = 0, pos = ylim[1])
axis(side = 1, (0:4) * s, labels = 0:4, pos = ylim[1])
# axis(side = 1, seq(xlim[1], xlim[2], by = 100), pos = -0.005) 
axis(side = 1, seq(xlim[1], xlim[2], by = 100), pos = -0.0125) 
axis(side = 2, ylim, labels = FALSE, lwd.ticks = 0, pos = xlim[1])
# axis(side = 2, c(0, 0.01, 0.02), pos = xlim[1], las = TRUE)
axis(side = 2, c(0, 0.02, 0.04, 0.06), pos = xlim[1], las = TRUE)
par(pp)
dev.off()

#### Turning angle
# (For speed, for this plot we just amend code from patter-flapper)
# Define data
drummond_angle <- vps_angle[site == "Drummond", ]
thunder_angle  <- vps_angle[site == "Thunder", ]
# Make plot 
png(here_fig("model-move-turning-angle.png"), 
    height = 4, width = 4, units = "in", res = 800)
pp <- par(mgp = c(3, 0.7, 0))
x <- seq(-pi, pi, length.out = 1e5)
y <- dmix(x, 0, pars_model_move_full$phi[1])
ylim <- c(0, 1.2)
plot(x, y,
     ylim = ylim,
     xlab = "", ylab = "",
     type = "n", 
     axes = FALSE)
# Add observed distributions from VPS
add_poly(density(drummond_angle$turn_angle, from = -pi, to = pi), col = scales::alpha("lightblue", 1))
add_poly(density(thunder_angle$turn_angle, from = -pi, to = pi), col = scales::alpha("blue", 0.5))
# Add model
lines(x, y, lwd = 2)
y <- dmix(x, 0, pars_model_move_full$phi[2])
lines(x, y, col = "red", lty = 1, lwd = 1)
y <- dmix(x, 0, pars_model_move_full$phi[3])
lines(x, y, col = "darkgreen", lty = 1, lwd = 1)
axis(side = 1, at = c(-pi, pi), labels = FALSE, lwd.tick = 0, pos = 0)
axis(side = 1, 
     at = c(-pi, -pi/2, 0, pi/2, pi), 
     labels = c(expression(-pi), expression(-pi/2), expression(0), expression(pi/2), expression(pi)), 
     pos = 0)
yat <- seq(ylim[1], ylim[2], by = 0.4)
axis(side = 2, prettyGraphics:::add_lagging_point_zero(yat), las = TRUE, pos = -pi)
par(pp)
dev.off()


###########################
###########################
#### Write parameters to file

# Parameters
qs::qsave(pars_adj, here_input("pars-adj.qs"))
qs::qsave(pars_model_move_best, here_input("pars-model-move-best.qs"))
qs::qsave(pars_model_move_full, here_input("pars-model-move-full.qs"))

# vmaps
dirs.create(here_input("vmap", pars_model_move_full$mobility))
pp <- par(mfrow = c(1, nrow(pars_model_move_full)))
lapply(split(pars_model_move_full, seq_len(nrow(pars_model_move_full))), function(d) {
  vmap <- patter:::spatVmap(.map = map, .mobility = d$mobility, .plot = TRUE)
  terra::writeRaster(vmap, 
                     here_input("vmap", d$mobility, "vmap.tif"), 
                     overwrite = TRUE)
})
par(pp)


#### End of code. 
###########################
###########################