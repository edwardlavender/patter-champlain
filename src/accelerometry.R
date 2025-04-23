# Calculate SS (BL/s), given effects
# * Constants alpha, beta
# * Vectors zeta (random intercept) and eps (random error)
# * Vector log10_A observations
calc_SS <- function(alpha, beta, zeta, eps, log10_A) {
  log10_SS <- alpha + beta * log10_A + zeta + eps
  10^log10_SS
}

# Simulate speeds (BL/s)
# * This function returns a matrix
# * For each observation (row), we have multiple simulated SS
sim_SS <- function(data,
                   alpha, alpha_se,
                   beta, beta_se,
                   corr,
                   sigma_id, sigma,
                   n_sim = 250L,
                   cl = 1L) {
  
  #### Prepare data
  n    <- nrow(data)
  data <- data[, list(id, log10_A)]
  
  #### Simulate alpha/beta for each simulation
  # (We do this outside the loop below for speed)
  if (n_sim == 1L) {
    fixed   <- matrix(c(alpha, beta), ncol = 2)
  } else {
    cov     <- corr * alpha_se * beta_se
    cov_mat <- matrix(c(alpha_se^2, cov, cov, beta_se^2), nrow = 2)
    fixed   <- MASS::mvrnorm(n = n_sim, mu = c(alpha, beta), Sigma = cov_mat)
  }
  
  #### Simulate a list of SS values with one element for n_sim
  sims <- cl_lapply(seq_len(n_sim), .cl = cl, .fun = function(i) {
    
    # Extract simulated fixed effects
    sim_alpha <- fixed[i, 1]
    sim_beta  <- fixed[i, 2]
    
    # Simulate random effects & noise
    sim_data <-
      data |>
      lazy_dt() |>
      # Simulate individual-specific random effects
      group_by(id) |>
      mutate(zeta = rnorm(1L, mean = 0, sd = sigma_id)) |>
      ungroup() |>
      # Simulate random noise
      mutate(eps = rnorm(n, mean = 0, sd = sigma)) |>
      as.data.table()
    
    # Calc SS using simulated parameters
    calc_SS(alpha   = sim_alpha,
            beta    = sim_beta,
            zeta    = sim_data$zeta,
            eps     = sim_data$eps,
            log10_A = data$log10_A)
    
  })
  
  # Convert to matrix
  # * Each row is an observation
  # * Each column is a realisation
  do.call(cbind, sims)
  
}