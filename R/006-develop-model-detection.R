###########################
###########################
#### setup-model-detection.R

#### Aims
# 1) Setup the detection probability model 

#### Prerequisites
# 1) Process Klinard et al. (2019) dataset


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Set global options
Sys.setenv("JULIA_SESSION" = FALSE)

#### Load essential packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(DHARMa)
library(ggplot2)
library(mgcv)
library(mgcViz)
library(proj.verse)
files_source_r(here_src())

#### Load data
klinard <- qs::qread(here_data("supp", "model-obs", "klinard.qs"))


###########################
###########################
#### Explore curves
# This code explores possible 

#### Initial 'best-guess' models 
# Summer model 
curve(plogis(1.8 + -0.005 * x), from = 0, to = 6000, ylim = c(0, 1))
# Winter model
curve(plogis(3 + -0.0025 * x), from = 0, to = 6000, col = "blue", add = TRUE)
# Compromise model
curve(plogis(2.5 + -0.003 * x), from = 0, to = 6000, col = "darkgreen", add = TRUE)

#### ggplot representation
alpha          <- 2.5    # intercept
beta           <- -0.003 # rate of decline (larger values, nearer 0, increase steepness)
receiver_gamma <- 7000
ggplot(data.frame(x = c(0, receiver_gamma)), aes(x = x)) +
  stat_function(fun = function(x) {
    prob <- plogis(alpha + beta * x)
    prob[x > receiver_gamma] <- 0
    prob
  }) +
  ylim(0, 1) + 
  scale_x_continuous(breaks = seq(0, receiver_gamma, by = 500))  + 
  theme_bw()


###########################
###########################
#### Analyse Klinard et al. (2019) datasets

#### Model detection probability

# detections ~ B(n, p)
# p = logistic(dist) or p = s(dist)

# Unweighted GLM
m1 <- glm(cbind(success, failure) ~ dist,
          data = klinard, family = binomial())

# Weighted GLM
m2 <- glm(cbind(success, failure) ~ dist,
          data = klinard, family = binomial(), weights = w)

# Unweighted GAM
m3 <- gam(cbind(success, failure) ~ s(dist), 
          data = klinard, family = binomial)

# Weighted GAM 
m4 <- gam(cbind(success, failure) ~ s(dist), 
          data = klinard, family = binomial, weights = w)

#### Extract GLM coefficients
# Model 2 is our prefered model (weighted GLM)
equatiomatic::extract_eq(m2)
(receiver_alpha <- coef(m2)[1]) # 1.885708
(receiver_beta  <- coef(m2)[2]) # -0.001613148
dbinom(1, size = 1, prob = plogis(receiver_alpha + receiver_beta * 8000))
dbinom(1, size = 1, prob = plogis(receiver_alpha + receiver_beta * 8001))

#### Visualise models
## (A) Compute predictions
nd   <- data.frame(dist = seq(0, 1e4, length.out = 1e3L))
fit  <- data.frame(dist = nd$dist, 
                   y0 = plogis(2.5 + -0.003 * nd$dist), # initial guess, 
                   y1 = predict(m1, newdata = nd, type = "response"), 
                   y2 = predict(m2, newdata = nd, type = "response"), 
                   y3 = predict(m3, newdata = nd, type = "response"),
                   y4 = predict(m4, newdata = nd, type = "response"))
## (B) Visualise models
# The GLMs are more 'generous' than our initial guess
# The GLMs and GAMs match well
# The GAMs better capture low-probability detections at higher distances
# The GAMs behave more poorly beyond range of data
# Weighted/unweighted models are similar
ggplot(klinard, aes(x = dist, y = prop)) +
  geom_bin_2d(bins = 50) +
  scale_fill_viridis_c(name = "Count") +
  geom_point(shape = ".") + 
  geom_line(data = fit, aes(x = dist, y = y0),
            lwd = 1.5, color = "grey", inherit.aes = FALSE) +
  geom_line(data = fit, aes(x = dist, y = y1),
            lwd = 1.5, color = "red", inherit.aes = FALSE) +
  geom_line(data = fit, aes(x = dist, y = y2), 
            lwd = 1.5, color = "darkred", inherit.aes = FALSE) +
  geom_line(data = fit, aes(x = dist, y = y3), 
            lwd = 1.5, color = "skyblue", inherit.aes = FALSE) +
  geom_line(data = fit, aes(x = dist, y = y4), 
            lwd = 1.5, color = "blue", inherit.aes = FALSE) +
  scale_x_continuous(limits = c(0, 1e4), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
  labs(x = "Distance", y = "Detection Probability") +
  theme_bw()

#### Examine residuals
# Residual diagnostics are poor
# But MLE parameter estimates look reasonable (above)
# and uncertainty quantification is not the aim here
r1 <- simulateResiduals(m1)
r2 <- simulateResiduals(m3) 
plot(r1)
plot(r2)

#### Examine receiver_gamma
kmax <- 
  klinard |> 
  group_by(transmitter_id) |> 
  mutate(max_dist = max(dist)) |> 
  slice(1L) |>
  select(transmitter_id, dB, max_dist) |> 
  as.data.table()
# Adjusted max detection ranges for 147 dB tag
# (Amplitude Distance Law)
kmax$max_dist * 10^((147 - kmax$dB) / 20)


###########################
###########################
#### Record parameters

#### Define 'best-guess' parameters (list)
pars_model_obs_best <- list(receiver_alpha = receiver_alpha, 
                            receiver_beta  = receiver_beta, 
                            receiver_gamma = 8000)

#### Define restrictive/flexible parameter combinations
# We assume these are known
# To minimise computation time, we only explore the effects of uncertainty in movement
# We find this more interesting

#### Collect all parameters (data.table)
pars_model_obs_full <- as.data.table(pars_model_obs_best)


###########################
###########################
#### Publication-quality plot

png(here_fig("model-obs.png"), 
    height = 4, width = 6, units = "in", res = 800)
gg <- 
  ggplot(klinard, aes(x = dist, y = prop)) +
  geom_bin_2d(bins = 50) +
  scale_fill_viridis_c(name = "Count", direction = -1, alpha = 0.95, 
                       guide     = guide_colorbar(
                         # draw a frame around the bar
                         frame.colour    = "black",
                         frame.linewidth = 0.5,
                         # draw ticks and labels
                         ticks           = TRUE,
                         ticks.colour    = "black",
                         ticks.linewidth = 0.5,
                         # size of the bar
                         barwidth        = unit(0.5, "cm"),
                         barheight       = unit(4,   "cm"),
                         # put title on top, labels beneath
                         title.position  = "top",
                         label.position  = "right"
                       )) +
  geom_point(shape = ".") + 
  # geom_line(data = fit, aes(x = dist, y = y0),
  #           lwd = 1.5, color = "grey", inherit.aes = FALSE) +
  # geom_line(data = fit, aes(x = dist, y = y1),
  #           lwd = 1.5, color = "red", inherit.aes = FALSE) +
  geom_line(data = fit, aes(x = dist, y = y2),
            lwd = 1.25, color = "black", inherit.aes = FALSE) +
  # geom_line(data = fit, aes(x = dist, y = y3), 
  #           lwd = 1.5, color = "skyblue", inherit.aes = FALSE) +
  geom_line(data = fit, aes(x = dist, y = y4), 
            lwd = 1.25, color = "dimgrey", inherit.aes = FALSE) +
  scale_x_continuous(limits = c(0, pars_model_obs_best$receiver_gamma), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
  labs(x = "Distance", y = "Detection probability") +
  theme_bw() + 
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.x.top  = element_blank(),
    axis.text.y.right = element_blank(),
    axis.title.x = element_text(size = 14, colour = "black", margin = margin(t = 7.5)),
    axis.title.y = element_text(size = 14, colour = "black", margin = margin(r = 7.5)),
    axis.text.x  = element_text(size = 12, colour = "black"),
    axis.text.y  = element_text(size = 12, colour = "black"), 
    axis.ticks.length = unit(0.3, "cm"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.title = element_text(color = "black"),
    legend.text  = element_text(color = "black")
  )
print(gg)
dev.off()
print(gg)

###########################
###########################
#### Write to file

qs::qsave(pars_model_obs_best, here_input("pars-model-obs-best.qs"))
qs::qsave(pars_model_obs_full, here_input("pars-model-obs-full.qs"))


#### End of code. 
###########################
###########################