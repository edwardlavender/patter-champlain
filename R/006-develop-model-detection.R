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
library(proj.verse)
files_source_r(here_src())

#### Load data
# TO DO


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
#### Update moorings

# TO DO
moorings[, receiver_alpha := 2.25]
moorings[, receiver_beta := -0.0022]
moorings[, receiver_gamma := 7000]


#### End of code. 
###########################
###########################