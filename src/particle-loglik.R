# Computes the log-likelihood given input parameters
pf_filter_loglik <- function(.theta, .sim) {
  
  #### Check inputs
  stopifnot(length(.theta) == 2L)
  shape <- .theta[1]
  scale <- .theta[2]
  
  #### Handle impossible parameters
  if (shape < 0 || scale < 0) {
    return(-Inf)
  }
  
  #### Prepare filter arguments
  args <- constructor_ac_sim(.sim = .sim, .datasets = list(), .verbose = TRUE)
  
  #### Update args
  args          <- args$forward
  args$.collect <- TRUE
  args$.batch   <- NULL
  
  #### Run filter
  # (Running pf_filter directly is simpler than estimate_coord_particle
  fwd  <- do.call(patter::pf_filter, args)
  
  #### Return log-lik
  # (-Inf is returned for convergence failures)
  fwd$callstats$loglik
  
}

# Wrapper function for iterative applications
pf_filter_loglik_optim <- function(.sim) {
  stopifnot(nrow(.sim) == 1L)
  proj.build::check_names(.sim, c("shape", "scale"))
  # pf_filter_loglik(.theta = c(.sim$shape, .sim$scale), .sim = .sim)
  optim(par = c(.sim$shape, .sim$scale),
        fn = pf_filter_loglik, .sim = .sim,
        hessian = TRUE,
        control = list(fnscale = -1))
}