using Patter

function run_particle_filter(; iter, 
                               env_init, 
                               timeline, 
                               state, 
                               model_move, 
                               datasets, 
                               model_obs_types, 
                               yobs, 
                               n_particle, 
                               n_record,
                               t_resample, 
                               direction, 
                               batch)

  # Initialise filter 
  xinit = simulate_states_init(map            = env_init,
                               timeline        = timeline,
                               state_type      = state,
                               xinit           = nothing,
                               model_move      = model_move,
                               datasets        = datasets,
                               model_obs_types = model_obs_types,
                               n_particle      = n_particle,
                               direction       = direction,
                               output          = "Vector")

  # Run the forward filter
  out = particle_filter(timeline   = timeline,
                        xinit      = xinit,
                        yobs       = yobs,
                        model_move = model_move,
                        n_move     = iter.n_move,
                        n_resample = iter.n_resample,
                        t_resample = t_resample,
                        n_record   = n_record,
                        direction  = direction, 
                        batch      = batch,
                        progress   = Patter.progress_control(enabled = isinteractive()),
                        verbose    = isinteractive());

  return out

end 