model_obs_champlain <- function(.moorings, .pars, .as_ModelObs = TRUE) {
  proj.build::check_names(.pars, c("receiver_alpha", "receiver_beta", "receiver_gamma"))
  moorings <- copy(.moorings)
  moorings[, receiver_alpha := .pars$receiver_alpha]
  moorings[, receiver_beta := .pars$receiver_beta]
  moorings[, receiver_gamma := .pars$receiver_gamma]
  if (!.as_ModelObs) {
    # Return moorings (includes receiver_id column)
    return(moorings)
  } else {
    # Return ModelObs object (includes sensor_id column)
    return(model_obs_acoustic_logis_trunc_los(moorings))
  }
}