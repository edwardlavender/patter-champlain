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

#### Define studies
# Define studies
studies <- sort(unique(dcounts$study))
# Check balance of observations
dcounts |> 
  group_by(study) |> 
  summarise(sum(success + failure))

#### Model detection probability
# detections ~ B(n, p)
# p = logistic(dist)

models  <- cl_lapply(studies, function(s) {
  
  # Define data 
  # s <- studies[1]
  data <- as.data.frame(dcounts[study == s, ])
  
  # GLM 
  m1 <- glm(cbind(success, failure) ~ dist,
      data = data, 
      family = binomial())
  
  # GAM
  m2 <- gam(cbind(success, failure) ~ dist, 
            data = dcounts, family = binomial)
  
  list(GLM = m1, GAM = m2)
})
names(models) <- studies

# SCAM (enforce monotonic decline)
# m_scam <- scam::scam(cbind(success, failure) ~ study + 
#                  s(dist, by = study, k = 10, bs = "mpd"), 
#                data = dcounts, family = binomial)


#### Extract example GLM coefficients
# Here we extract intercept & distance coefficient for one of the studies
m1 <- models[[1]]$GLM
equatiomatic::extract_eq(m1)
(receiver_alpha <- coef(m1)[1]) 
(receiver_beta  <- coef(m1)[2])
dbinom(1, size = 1, prob = plogis(receiver_alpha + receiver_beta * 7000))
dbinom(1, size = 1, prob = plogis(receiver_alpha + receiver_beta * 7001))

#### Compute predictions
# Define data.table of studies & distances
# * Note that GAMs behave poorly beyond the range of the data
dists <- seq(0, 1e4, by = 1)
nd <- data.table(dist = dists)
# Generate predictions, for each model
pred_empirical <- 
  lapply(studies, function(s) {
    lapply(c("GLM", "GAM"), function(m) {
      # s = studies[[1]]; m <- "GLM"
      mod  <- models[[s]][[m]]
      p <- predict(mod, newdata = nd, se.fit = TRUE, type = "link")
      p <- prettyGraphics::list_CIs(p, inv_link = mod$family$linkinv, plot_suggestions = FALSE)
      cbind(nd, 
            data.table(study = s, 
                       model = m,
                       fit = as.numeric(p$fit), 
                       lwr = as.numeric(p$lowerCI), 
                       upr = as.numeric(p$upperCI)))
    }) |> rbindlist()
  }) |> 
  rbindlist() |>
  arrange(study, model, dist) |> 
  as.data.table()

#### Visualise GLMs
ggplot(dcounts, aes(x = dist, y = prop)) +
  geom_bin_2d(bins = 50) +
  scale_fill_viridis_c(name = "Count") +
  geom_point(shape = ".") + 
  geom_line(data = pred_empirical[model == "GLM", ], 
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

# > Black line -> detection pr was too high
# > We're trying to use one model for all receivers (limitation)
# > Red line may be biased low b/c of range test design
# > Green lines are more reasonable choices

#### Visualise GLMs versus GAMS
head(dcounts)
head(pred_empirical)
ggplot() +
  geom_bin_2d(aes(x = dist, y = prop), data = dcounts, bins = 50) +
  scale_fill_viridis_c(name = "Count") +
  geom_point(shape = ".") + 
  geom_line(data = pred_empirical,
            aes(x = dist, y = fit, colour = model, group = model), 
            lwd = 1.25, inherit.aes = FALSE) + 
  facet_wrap(~study)

#### Examine residuals
# Residual diagnostics are poor
# But MLE parameter estimates look reasonable (above)
# and uncertainty quantification is not the aim here
if (FALSE) {
  lapply(models, function(mods) {
    m1 <- mods[["GLM"]]
    m2 <- mods[["GAM"]]
    r1 <- simulateResiduals(m1)
    r2 <- simulateResiduals(m2) 
    readline("See GLM...")
    plot(r1)
    readline("See GAM...")
    plot(r2)
    invisible(NULL)
  }) |> invisible()
}

#### Examine estimated 50 % detection range
pred_empirical |> 
  filter(model == "GLM") |> 
  group_by(study) |> 
  slice(which.min(abs(0.5 - fit)))
  

#### Examine empirical maximum detection range
# * Max detection range may be affected by study design
# * NB Futia testing does not cover entire Lake Champlain 
dcounts |> 
  group_by(study) |> 
  mutate(max_dist = max(dist)) |> 
  slice(1L) |>
  select(study, transmitter_id, dB, max_dist) |> 
  as.data.table()


###########################
###########################
#### Record parameters

#### Define 'best-guess' parameters (list)
# Pull out coefficients for suitable model
coefs <- coef(models[["F151"]][["GLM"]])
# Define parameters
a <- as.numeric(coefs[1])
b <- as.numeric(coefs[2])
g <- 7000
# Collate parameter list
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

#### Define detection probability curves
pm <- copy(pars_model_obs_full)
pm[, sensitivity_label := c("Best", "Restricted", "Flexible")]
pred_patter <- 
  lapply(split(pm, seq_len(nrow(pm))), function(d) {
  dist <- nd$dist
  fit <- trunclogis(d$receiver_alpha, d$receiver_beta, d$receiver_gamma, dist)
  data.table(sensitivity_label = d$sensitivity_label, 
             dist = dist, 
             fit = fit)
}) |> rbindlist()


###########################
###########################
#### Publication-quality plot

#### Make plot
png(here_fig("model", "model-obs", "detection-probability.png"), 
    height = 5, width = 10, units = "in", res = 800)
gg <- 
  ggplot(dcounts, aes(x = dist, y = prop)) +
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
  # geom_point(shape = ".") + 
  # Add best, restrictive and flexible models
  geom_line(data = pred_patter, aes(x = dist, y = fit, 
                                      colour = sensitivity_label, group = sensitivity_label),
            lwd = 1.75,  inherit.aes = FALSE) +
  # Add GLMs
  geom_line(data = pred_empirical[model == "GLM", ], 
            aes(x = dist, y = fit, colour = study, group = study), 
            lwd = 0.75, linetype = 2, inherit.aes = FALSE) +
  # Colour lines 
  scale_colour_manual(values = c(
      "F146" = "lightblue",
      "F151" = "skyblue",
      "F152" = "blue",
      "K145" = "mediumpurple1",
      "K153" = "purple3",
      "Best"        = "black",
      "Restricted"  = "red",
      "Flexible"    = "darkgreen"
    )
  ) + 
  # Axes
  scale_x_continuous(limits = c(0, max(pars_model_obs_full$receiver_gamma)), expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1.025), clip = "off") +
  # receiver_gamma (added after axes)
  annotate(
    "segment",
    x     = pars_model_obs_full$receiver_gamma,
    xend  = pars_model_obs_full$receiver_gamma,
    y     = -0.05, 
    yend  = -0.0075,
    arrow = arrow(length = unit(0.2, "cm")),
    colour = c("black", "red", "darkgreen"), 
    linewidth = 1
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