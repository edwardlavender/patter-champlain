###########################
###########################
#### setup-model-detection.R

#### Aims
# 1) Setup the detection probability model 

#### Prerequisites
# 1) Process detection probability dataset


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
dcounts  <- qs::qread(here_data("supp", "model-obs", "range-testing.qs"))
pars_adj <- qs::qread(here_input("pars-adj.qs"))


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
#### Analyse range-testing datasets

#### Model detection probability

# detections ~ B(n, p)
# p = logistic(dist * study) or p = s(dist * study)

# GLM
m1 <- glm(cbind(success, failure) ~ dist * study,
          data = dcounts, family = binomial())

# GAM
# * Model I in Pedersen et al. (2019)
# * Complete flexibility: different smoothness, different penalties 
m2 <- gam(cbind(success, failure) ~ study + s(dist, by = study, k = 5), 
          data = dcounts, family = binomial)

# SCAM (enforce monotonic decline)
# m_scam <- scam::scam(cbind(success, failure) ~ study + 
#                  s(dist, by = study, k = 10, bs = "mpd"), 
#                data = dcounts, family = binomial)


#### Extract example GLM coefficients
# Model 1 is our main model (GLM)
# Here we extract intercept & distance coefficient for one of the studies
equatiomatic::extract_eq(m1)
(receiver_alpha <- coef(m1)[1]) # 1.885708
(receiver_beta  <- coef(m1)[2]) # -0.001613148
dbinom(1, size = 1, prob = plogis(receiver_alpha + receiver_beta * 8000))
dbinom(1, size = 1, prob = plogis(receiver_alpha + receiver_beta * 8001))

#### Compute predictions
# Define data.table of studies & distances
# * Note that GAMs behave poorly beyond the range of the data
dists <- seq(0, 1e4, by = 1)
nd <- lapply(unique(dcounts$study), function(s) {
  data.table(study = s, dist = seq(0, max(dcounts$dist[dcounts$study == s]), by = 1))
  # data.table(study = s, dist = dists)
}) |> rbindlist()
# Generate predictions, for each model
pred <- lapply(1:2, function(i) {
    ms <- list(Best = m1, GAM = m2)
    m  <- ms[[i]]
    p <- predict(m, newdata = nd, se.fit = TRUE, type = "link")
    p <- prettyGraphics::list_CIs(p, inv_link = m$family$linkinv, plot_suggestions = FALSE)
    cbind(nd, 
          data.table(model = names(ms)[i], 
                     fit = as.numeric(p$fit), 
                     lwr = as.numeric(p$lowerCI), 
                     upr = as.numeric(p$upperCI)))
  }) |> 
  rbindlist() |>
  arrange(study, model, dist) |> 
  as.data.table()

#### Visualise GLMs
# The GLMs are more 'generous' than our initial guess
# The GLMs and GAMs match well
# For K153, the GAMs better capture low-probability detections at higher distances
# Otherwise, GAMs behave more pooly at the edges of the data
ggplot(dcounts, aes(x = dist, y = prop)) +
  geom_bin_2d(bins = 50) +
  scale_fill_viridis_c(name = "Count") +
  geom_point(shape = ".") + 
  geom_line(data = pred[model == "Best", ], 
            aes(x = dist, y = fit, colour = study, group = study), 
            lwd = 1.25, inherit.aes = FALSE) +
  # Add preliminary best-guess
  # * After a preliminary GLM analysis based on the Klinard et al. data
  #   we generated the following line. We found this was too generous
  #   (receivers 'blocked' movements between sequential detections). 
  geom_line(data = data.table(dist = dists, 
                              fit = plogis(1.885708 - 0.001613148 * dists)), 
            aes(x = dist, y = fit),
            col = "black", lwd = 1.5, inherit.aes = FALSE) + 
  scale_x_continuous(limits = c(0, 1e4), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
  # facet_wrap(~model) + 
  labs(x = "Distance", y = "Detection Probability") +
  theme_bw()

# Green line
# * start black -> detection pr was too high
# * we're trying to use one model for all receivers (limitation)
# * red line may be biased low b/c of range test design

#### Visualise GLMs versus GAMS
head(dcounts)
head(pred)
ggplot() +
  geom_bin_2d(aes(x = dist, y = prop), data = dcounts, bins = 50) +
  scale_fill_viridis_c(name = "Count") +
  geom_point(shape = ".") + 
  geom_line(data = pred,
            aes(x = dist, y = fit, colour = model, group = model), 
            lwd = 1.25, inherit.aes = FALSE) + 
  facet_wrap(~study)

#### Examine residuals
# Residual diagnostics are poor
# But MLE parameter estimates look reasonable (above)
# and uncertainty quantification is not the aim here
r1 <- simulateResiduals(m1)
r2 <- simulateResiduals(m2) 
plot(r1)
plot(r2)

#### Examine receiver_gamma
# * Max detection range may be affected by study design
# * NB Futia testing does not cover entire Lake Champlain 
dcounts |> 
  group_by(study) |> 
  mutate(max_dist = max(dist)) |> 
  slice(1L) |>
  select(study, transmitter_id, dB, max_dist) |> 
  as.data.table()

#### Compute balance of observations
dcounts |> 
  group_by(study) |> 
  summarise(sum(success + failure))


###########################
###########################
#### Record parameters

#### Define 'best-guess' parameters (list)
# Define parameters
a <- receiver_alpha
b <- receiver_beta
g <- 8000.0  # use 7000.0
# Collate in list
pars_model_obs_best <- list(receiver_alpha = a, 
                            receiver_beta  = b, 
                            receiver_gamma = g)

#### Collect all parameters (data.table)
# We use the same degree of uncertainty as for the movement model
# For restrictive model: deflate alpha, inflate beta
# For flexible model: inflate alpha, deflate beta
inflate <- pars_adj$inflate
deflate <- pars_adj$deflate
pars_model_obs_full <- data.table(receiver_alpha = c(a, a * deflate, a * inflate), 
                                  receiver_beta = c(b, b * inflate, b * deflate), 
                                  receiver_gamma = c(g, g * deflate, g * inflate))


###########################
###########################
#### Publication-quality plot

#### Define datasets
# Best-guess (based on weighted GLM): y2
# GAM (comparison)                  : y4
# * Define above 
# Restrictive model                 : y5
# Flexible model                    : y6
head(fit)
p <- pars_model_obs_full
fit[dist > p$receiver_gamma[1], y2 := 0]
fit[, y5 := plogis(p$receiver_alpha[2] + p$receiver_beta[2] * dist)]
fit[dist > p$receiver_gamma[2], y5 := NA]
fit[, y6 := plogis(p$receiver_alpha[3] + p$receiver_beta[3] * dist)]
fit[dist > p$receiver_gamma[3], y5 := NA]
rm(p)

#### Make plot
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
  # Weighted GAM (baseline)
  geom_line(data = fit, aes(x = dist, y = y4), 
            lwd = 1, colour = "dimgrey", inherit.aes = FALSE) +
  # Best model (weighted GLM, truncated)
  geom_line(data = fit, aes(x = dist, y = y2),
            lwd = 1.25, colour = "black", inherit.aes = FALSE) +
  # Restrictive and flexible models (truncated)
  geom_line(data = fit, aes(x = dist, y = y5), 
            lwd = 0.75, colour = "red", inherit.aes = FALSE) +
  geom_line(data = fit, aes(x = dist, y = y6), 
            lwd = 0.75, colour = "darkgreen", inherit.aes = FALSE) +
  # Axes
  scale_x_continuous(limits = c(0, max(pars_model_obs_full$receiver_gamma)), expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1.025), clip = "off") +
  # receiver_gamma (added after axes)
  annotate(
    "segment",
    x     = pars_model_obs_full$receiver_gamma,
    xend  = pars_model_obs_full$receiver_gamma,
    y     = -0.04, 
    yend  = -0.0075,
    arrow = arrow(length = unit(0.15, "cm")),
    colour = c("black", "red", "darkgreen")
  ) +
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