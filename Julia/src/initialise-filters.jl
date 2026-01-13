import Rasters
import LibGEOS
using DataFrames
using Dates 
using Distributions
using JLD2
using Patter

# To initialise the forward particle filter we:
# (1) Sample particles uniformally within the portion of the study area compatible with initial observations, if detection at t = 1
# (2) Sample particles by running the filter backwads in time from first detection container to the start of the timeline, accounting for non detections
#     - This method improves convergence when there is a gap between the start of the timeline and the first detection
#     - Particles in areas without receivers have better survival properties (in South Lake and North of study area) b/c non-detection barriers between basins
#     - This causes convergence issues in particle filter
#     - This method leverages information in non-detections to define more informative starting locations 

function initialise_forward_filter(iter, env_init, timeline, state, model_move, datasets_fwd, model_obs_types, acoustics)

    # Isolate first detection
    detections  = acoustics[acoustics.obs .== 1, :]
    detections0 = detections[1, :]

    # If there is a detection at the first time step:
    # sample particles uniformally within the portion of the study area compatible with initial observations (default method)
    xinit = nothing
    if minimum(timeline) == detections0.timestamp

        xinit = simulate_states_init(map             = env_init,
                                     timeline        = timeline,
                                     state_type      = state,
                                     xinit           = nothing,
                                     model_move      = model_move,
                                     datasets        = datasets_fwd,
                                     model_obs_types = model_obs_types,
                                     n_particle      = iter.n_particle_filter,
                                     direction       = "forward",
                                     output          = "Vector")

    # Otherwise, run filter backwards in time from first detection to start of timeline
    else

        # Define timeline0 from timeline[1] to time of first detection
        timeline0 = collect(minimum(timeline):Minute(2):detections0.timestamp)

        # Define map0 for first detection
        container0 = LibGEOS.buffer(LibGEOS.Point(detections0.receiver_x, detections0.receiver_y),
                                    detections0.receiver_gamma,
                                    1000)
        env0       = Rasters.mask(env_init, with = container0)

        # Define initial states, accounting for line of sight 
        coord0 = Patter.spatSample(x = env0, size = iter.n_particle_filter, drop_missing = true)
        los    = [in_line_of_sight(env, detections0.receiver_x, detections0.receiver_y, coord0.x[i], coord0.y[i]) for i in 1:nrow(coord0)]
        coord0 = coord0[los .== true, :]
        if nrow(coord0) != iter.n_particle_filter
            coord0 = coord0[sample(1:nrow(coord0), iter.n_particle_filter, replace = true), :]
        end
        (nrow(coord0) == iter.n_particle_filter) || error("coord0 does not contain iter.n_particle_filter rows!")
        xinit0 = Patter.states_init(StateCXY, coord0)
        xinit0 = Patter.julia_get_xinit(StateCXY, xinit0)

        # Define dataset (non detections)
        nondetections       = acoustics[acoustics.timestamp .<= detections0.timestamp, :]
        nondetections.obs  .= 0
        yobs0               = assemble_yobs(datasets = [nondetections],
                                            model_obs_types = [ModelObsAcousticLogisTruncLos])

        # Define batches
        # We use b batches to store values v for N particles * T/b time steps (where n = n_particle_smoother)
        # I.e., v = N * T/b 
        # To maintain this memory usage, we need b = N * T / v batches here
        v = iter.n_particle_smoother * ceil(Int, (length(timeline) / iter.n_batch))
        b = iter.n_particle_smoother * length(timeline) / v
        b = ceil(Int, min(b, 9, length(timeline0)))
        n_batch0     = ceil(Int, length(timeline0) / ceil(Int, length(timeline) / iter.n_batch))
        bwd_batches0 = [joinpath(iter.folder_output, "bwd0-$i.jld2") for i in 1:n_batch0]

        # Run filter backwards from first detection to start of timeline, assuming only non-detections
        # * We only record iter.n_particle_smoother particles to maintain memory use 
        bwd0 = particle_filter(timeline   = timeline0,
                               xinit      = xinit0,
                               yobs       = yobs0,
                               model_move = model_move,
                               n_move     = 1000,
                               n_resample = iter.n_resample,
                               n_record   = iter.n_particle_smoother, # iter.n_particle_filter,
                               direction  = "backward",
                               batch      = bwd_batches0,
                               progress   = Patter.progress_control(enabled = isinteractive()),
                               verbose    = isinteractive())

        # Take the states at timestep[1] & resample for xinit
        # (rows = particles, columns = timesteps)
        # xinit = bwd0.states[:, 1]
        @load bwd_batches0[1] xbwd
        xinit = xbwd[:, 1]
        xbwd = nothing 
        foreach(f -> rm(f; force = true), bwd_batches0)
        xinit = xinit[sample(1:length(xinit), iter.n_particle_filter, replace = true)]
        xinit = [StateCXY(1.0, xinit[i].x, xinit[i].y, rand(Uniform(0.0, 2π))) for i in 1:iter.n_particle_filter]

    end

    return xinit

end

# To initialise the backward filter, we resample the final particles from the forward filter
function initialise_backward_filter(iter, fwd_batches)
    # Read ending particles from forward filter
    @load fwd_batches[end] xfwd
    fwd_states_end = xfwd[:, end]
    xfwd = nothing 
    # Define coordinates DataFrame (contains n_record rows)
    coord0 = DataFrame(map_value = [fwd_states_end[i].map_value for i in 1:iter.n_particle_smoother],
                       x         = [fwd_states_end[i].x for i in 1:iter.n_particle_smoother], 
                       y         = [fwd_states_end[i].y for i in 1:iter.n_particle_smoother])
    # Resample DataFrame so we have n_forward_filter rows
    coord0 = coord0[sample(1:nrow(coord0), iter.n_particle_filter, replace = true), :]
    # Add heading (for StateCXY)
    coord0.heading = rand(Uniform(0.0, 2π), iter.n_particle_filter)
    # Coerce to Vector{State}
    xinit = Patter.julia_get_xinit(StateCXY, coord0)
    return xinit
end 

