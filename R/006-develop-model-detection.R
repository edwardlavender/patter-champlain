###########################
###########################
#### setup-model-detection.R

#### Aims
# 1) Setup the detection probability model 

#### Prerequisites
# 1) TO DO


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())

#### Set global options
Sys.setenv("JULIA_SESSION" = FALSE)

#### Load essential packages
library(DHARMa)
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
m1 <- glm(cbind(success, failure) ~ dist,
          data = klinard, family = binomial("cloglog"))
m2 <- gam(cbind(success, failure) ~ s(dist), 
          data = klinard, family = binomial)

#### GLM coefficients
(receiver_alpha <- coef(m1)[1]) # 2.904267
(receiver_beta  <- coef(m1)[2]) # -0.001546906

#### Visualise models
# The GLM is more 'generous' than our initial guess
# The GLM and GAM match well
# The GAM better captures low-probability detections at higher distances
# The GAM behaves more poorly beyond range of data
dist <- seq(0, 1e4, length.out = 1e3L)
nd   <- data.frame(dist = dist)
y0   <- plogis(2.5 + -0.003 * dist) # initial guess
y1   <- predict(m1, newdata = nd, type = "response")
y2   <- predict(m2, newdata = nd, type = "response")
fit  <- data.frame(dist = dist, y0 = y0, y1 = y1, y2 = y2)
ggplot(klinard, aes(x = dist, y = prop)) +
  geom_bin_2d(bins = 50) +
  scale_fill_viridis_c(name = "Count") +
  geom_point(shape = ".") + 
  geom_line(data = fit, aes(x = dist, y = y0),
            lwd = 1.5, color = "grey", inherit.aes = FALSE) +
  geom_line(data = fit, aes(x = dist, y = y1),
            lwd = 1.5, color = "blue", inherit.aes = FALSE) +
  geom_line(data = fit, aes(x = dist, y = y2), 
            lwd = 1.5, color = "red", inherit.aes = FALSE) +
  scale_x_continuous(limits = c(0, 1e4), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
  labs(x = "Distance", y = "Detection Probability") +
  theme_bw()

#### Examine residuals
# Residual diagnostics are poor
# But MLE parameter estimates look reasonable
r1 <- simulateResiduals(m1)
r2 <- simulateResiduals(m2) 
plot(r1)
plot(r2)


###########################
###########################
#### Record parameters

#### Define 'best-guess' parameters (list)
pars_model_obs_best <- list(receiver_alpha = receiver_alpha, 
                            receiver_beta  = receiver_beta, 
                            receiver_gamma = 7500)

#### Define restrictive/flexible parameter combinations
# We assume these are known
# To minimise computation time, we only explore the effects of uncertainty in movement
# We find this more interesting

#### Collect all parameters (data.table)
pars_model_obs_full <- as.data.table(pars_model_obs_best)

#### Write to file
qs::qsave(pars_model_obs_best, here_input("pars-model-obs-best.qs"))
qs::qsave(pars_model_obs_full, here_input("pars-model-obs-full.qs"))


#### End of code. 
###########################
###########################