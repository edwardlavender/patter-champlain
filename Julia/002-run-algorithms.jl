###########################
###########################
#### run-algorithms.jl

#### Aims
# 1) Run the particle algorithms

#### Prerequisites
# 1) This code can be run interactively or via bash from SIA-LAVENDED or siam-linux20
#    (It is mainly designed for parallel implementation via bash on siam-linux20)


###########################
###########################
#### Set up

#### Load base packages
import Base.Threads
import Pkg

#### Load local packages
Pkg.activate(".")
using BenchmarkTools
import Arrow
import GeoArrays
import Random
using DataFrames
using Dates
using Distributions
using JLD2
using Patter
# using Plots

#### Load source files
include("./src/utils.jl")
include("./src/initialise-filters.jl")
include("./src/observation-model.jl")

#### Load datasets (map, iteration)
# Load map & iteration 
# analysis  = "sim"
analysis  = "real"
subanalysis = "main"
env       = GeoArrays.read(joinpath("data", "input", "map.tif"));
env_init  = Patter.rast(joinpath("data", "input", "map.tif"));
iteration = DataFrame(Arrow.Table(joinpath("data", "input", analysis, subanalysis, "iteration.feather")))
# iteration = iteration[iteration.index .∈ Ref([7, 13, 14, 119, 133, 140, 147, 161, 166, 168, 176, 178, 179, 182, 195, 196]), :];

#### Select iteration 
# Set column types as needed
iteration.n_move              = Int.(iteration.n_move);
iteration.n_resample          = Float64.(iteration.n_resample);
iteration.n_particle_filter   = Int.(iteration.n_particle_filter);
iteration.n_particle_smoother = Int.(iteration.n_particle_smoother);

# Select row 
if isinteractive()
    row = 1
else
    if length(ARGS) < 1
        error("Row index not provided")
    end
    row = try
        parse(Int, ARGS[1])
    catch
        error("Row index is NA")
    end
    if !(1 ≤ row ≤ size(iteration, 1))
        error("Invalid row index.")
    end
end
iter = iteration[row, :]

#### Define local settings
# Check JULIA_NUM_THREADS = 1
Threads.nthreads()
if (Threads.nthreads() != 1) && !isinteractive()
  error("JULIA_NUM_THREADS must be 1 for parallelisation!")
end 
# Set seed
Random.seed!(123);
start = now();


###########################
###########################
#### Define state-space model 


###########################
#### Define movement model

#### Define state 
state = StateCXY

#### Define movement model
model_move = ModelMoveCXY(env, 
                          iter.mobility, 
                          truncated(Gamma(iter.shape, iter.scale), upper = iter.mobility), 
                          MixtureModel([truncated(Normal(0.0, iter.phi), -pi, pi), Uniform(-pi, pi)], [0.99, 0.01]))


###########################
#### Define observation model 

#### Load timeline 
# Define timeline 
timeline           = DataFrame(Arrow.Table(iter.file_timeline))
timeline.timestamp = DateTime.(timeline.timestamp)
timeline           = timeline.timestamp
# Define time steps & timesteps_by_batch
timesteps          = collect(1:length(timeline))
timesteps_by_batch = Patter.split_indices(timesteps, Int(iter.n_batch))

#### Load acoustic observations & containers
acoustics      = DataFrame(Arrow.Table(iter.file_acoustics))
containers_fwd = DataFrame(Arrow.Table(iter.file_containers_fwd))
containers_bwd = DataFrame(Arrow.Table(iter.file_containers_bwd))

#### Process columns
# Timestamps
acoustics.timestamp      = DateTime.(acoustics.timestamp);
containers_fwd.timestamp = DateTime.(containers_fwd.timestamp);
containers_bwd.timestamp = DateTime.(containers_bwd.timestamp);
# Acoustic columns
acoustics.sensor_id      = Int.(acoustics.sensor_id);
acoustics.obs            = Int.(acoustics.obs);
acoustics.receiver_x     = Float64.(acoustics.receiver_x);
acoustics.receiver_y     = Float64.(acoustics.receiver_y);
acoustics.receiver_alpha = Float64.(acoustics.receiver_alpha);
acoustics.receiver_beta  = Float64.(acoustics.receiver_beta);
acoustics.receiver_gamma = Float64.(acoustics.receiver_gamma);
# containers_fwd
containers_fwd.obs        = Int.(containers_fwd.obs);
containers_fwd.sensor_id  = Int.(containers_fwd.sensor_id);
containers_fwd.centroid_x = Float64.(containers_fwd.centroid_x);
containers_fwd.centroid_y = Float64.(containers_fwd.centroid_y);
containers_fwd.radius     = Float64.(containers_fwd.radius);
# containers_bwd
containers_bwd.obs        = Int.(containers_bwd.obs);
containers_bwd.sensor_id  = Int.(containers_bwd.sensor_id);
containers_bwd.centroid_x = Float64.(containers_bwd.centroid_x);
containers_bwd.centroid_y = Float64.(containers_bwd.centroid_y);
containers_bwd.radius     = Float64.(containers_bwd.radius);
# # > This is necessary to avoid issues with data imported via Arrow
# > This is not necessary for data imported from CSV

#### Assemble datasets 
# Collate datasets & associated `ModelObs` instances into a typed dictionary 
datasets_fwd    = [acoustics, containers_fwd];
datasets_bwd    = [acoustics, containers_bwd];
model_obs_types = [ModelObsAcousticLogisTruncLos, ModelObsContainer];
yobs_fwd        = assemble_yobs(datasets = datasets_fwd,
                                model_obs_types = model_obs_types);
yobs_bwd        = assemble_yobs(datasets = datasets_bwd,
                                model_obs_types = model_obs_types);


###########################
###########################
#### Run particle algorithms 


###########################
#### Set up algorithms 

# Define batches 
fwd_batches = [joinpath(iter.folder_output, "fwd-$i.jld2") for i in 1:iter.n_batch]
bwd_batches = [joinpath(iter.folder_output, "bwd-$i.jld2") for i in 1:iter.n_batch]
smo_batches = [joinpath(iter.folder_output, "smo-$i.jld2") for i in 1:iter.n_batch]
pou_batches = [joinpath(iter.folder_output, "pou-$i.feather") for i in 1:iter.n_batch]

# Define output objects
fwd = bwd = smo = nothing 

# Define duration placeholders
td_fwd = td_bwd = td_smo = NaN


###########################
#### Forward filter 

#### Simulate initial states for the forward filter
xinit = initialise_forward_filter(iter, env_init, timeline, state, model_move, datasets_fwd, model_obs_types, acoustics)
# Plots.plot(env)
# scatter!([xinit[i].x for i in 1:iter.n_particle_filter], [xinit[i].y for i in 1:iter.n_particle_filter], markersize = 0.01)

#### Run the forward filter
t1_fwd = now()
fwd = particle_filter(timeline   = timeline,
                      xinit      = xinit,
                      yobs       = yobs_fwd,
                      model_move = model_move,
                      n_move     = iter.n_move,
                      n_resample = iter.n_resample,
                      n_record   = iter.n_particle_smoother,
                      direction  = "forward", 
                      batch      = fwd_batches,
                      progress   = Patter.progress_control(enabled = isinteractive()),
                      verbose    = isinteractive());
t2_fwd = now()
td_fwd = diffsecs(t2_fwd, t1_fwd)

#### Collect outputs
# (To conserve disk space, we do not record states)
callstats               = fwd.callstats
diagnostics             = fwd.diagnostics 
diagnostics.routine    .= callstats.routine
diagnostics.ncell_core .= NaN
diagnostics.ncell_home .= NaN


###########################
#### (2) Backward filter 

convergence = fwd.callstats.convergence[1]
if convergence

  #### Simulate initial states for the backward filter
  # Use states from forward filter
  xinit = initialise_backward_filter(iter, fwd_batches)
  # Plots.plot(env)
  # scatter!([xinit[i].x for i in 1:iter.n_particle_filter], [xinit[i].y for i in 1:iter.n_particle_filter], markersize = 0.1)

  #### Run the backward filter
  t1_bwd = now()
  bwd = particle_filter(timeline   = timeline,
                        xinit      = xinit,
                        yobs       = yobs_bwd,
                        model_move = model_move,
                        n_move     = iter.n_move,
                        n_resample = iter.n_resample,
                        n_record   = iter.n_particle_smoother,
                        direction  = "backward", 
                        batch      = bwd_batches,
                        progress   = Patter.progress_control(enabled = isinteractive()),
                        verbose    = isinteractive());
  t2_bwd = now()
  td_bwd = diffsecs(t2_bwd, t1_bwd)
  
  #### Collect outputs
  append!(callstats, bwd.callstats)
  bwd.diagnostics.routine    .= bwd.callstats.routine
  bwd.diagnostics.ncell_core .= NaN
  bwd.diagnostics.ncell_home .= NaN
  append!(diagnostics, bwd.diagnostics)
  convergence = bwd.callstats.convergence[1]

end 


###########################
#### Smoother

if convergence

  #### Define Monte Carlo settings 
  # Set n_sim = 0 and cache = nothing for unrestricted models (n_move = 1)
  vmap = nothing
  cache = false 
  n_sim = 0
  if iter.n_move > 1
    # vmap = GeoArrays.read(joinpath("data", "input", "vmap", string(Int(iter.mobility[1])), "vmap.tif"))
    vmap = GeoArrays.read(iter.file_vmap)
    cache = true 
    n_sim = 30
  end 

  #### Run smoother
  t1_smo = now()
  smo = particle_smoother_two_filter(timeline   = timeline,
                                     xfwd       = fwd_batches,
                                     xbwd       = bwd_batches,
                                     model_move = model_move,
                                     vmap       = vmap,
                                     n_particle = iter.n_particle_smoother,
                                     n_sim      = n_sim, 
                                     cache      = cache, 
                                     batch      = smo_batches, 
                                     progress   = Patter.progress_control(enabled = isinteractive()), 
                                     verbose    = isinteractive());
  t2_smo = now()
  td_smo = diffsecs(t2_smo, t1_smo)

  #### Collate smoother/summary outputs
  # This is implemented below 
                                     
end 

#### Clean up fwd_batches and bwd_batches
# This is implemented at the earliest possible stage
# We remove smo batch files after collating results 
foreach(f -> rm(f; force = true), fwd_batches)
foreach(f -> rm(f; force = true), bwd_batches)


###########################
###########################
#### Collate outputs 

if convergence

  #### (1) Collate callstats 
  callstats.n_iter .= Float64.(callstats.n_iter)
  append!(callstats, smo.callstats)

  #### (2) Collate state-dependent outputs by batch 
  # This is implemented batch-wise to manage memory
  areas_by_batch = Vector{DataFrame}(undef, iter.n_batch)
  for i in eachindex(smo_batches)
    
    # Define DataFrame of states for batch from Matrix (rows: particles; columns: time steps)
    @load smo_batches[i] xsmo
    smo_states = Patter.r_get_states(xsmo, timesteps_by_batch[i], timeline[timesteps_by_batch[i]])
    
    # (A) Compute grid cell weights, for mapping 
    # * This is a summarised DataFrame, with one row for each timestep & grid cell with the weight 
    # * We write this straight to file to keep memory usage low (over all batches)
    # * We collate DataFrames in R (summing weights over all time steps) for mapping 
    coord = map_marks(env, smo_states)
    Arrow.write(pou_batches[i], coord; compress = Arrow.ZstdCompressor(level = 9))

    # (B) Compute area (ncell) spanned by 50 % and 95 % of the probability mass, for diagnostics
    areas_by_batch[i] = map_uncertainty(coord)

  end

  #### Collate diagnostics over all batches, including area (ncell) spanned by 50 % and 95 % of the distribution 
  ncells                      = vcat(areas_by_batch...)
  smo.diagnostics.routine    .= smo.callstats.routine
  smo.diagnostics.ncell_core .= ncells.ncell_core
  smo.diagnostics.ncell_home .= ncells.ncell_home
  append!(diagnostics, smo.diagnostics)

  #### Collate smoothed DataFrame of states
  # This is too memory intensive for parallel applications
  # smo_batches = smo_batches[isfile.(smo_batches)]
  # if length(smo_batches) > 0
  #   # Collate states Matrix in Julia 
  #   smo_states = hcat([f["xsmo"] for f in map(jldopen, smo_batches)]...)
  #   # Convert to DataFrame 
  #   smo_states_df = Patter.r_get_states(smo_states, collect(1:length(timeline)), timeline)
  # end

end 


###########################
###########################
#### Record results

#### Cleanup smoothed batches
# This is implemented at the earliest possible stage
foreach(f -> rm(f; force = true), smo_batches)

#### Write outputs 

# (1) Write callstats 
Arrow.write(iter.file_callstats, callstats; compress = Arrow.ZstdCompressor(level = 9))

# (2) Write diagnostics 
Arrow.write(iter.file_diagnostics, diagnostics; compress = Arrow.ZstdCompressor(level = 9))

# (3) Write states (if convergence)
# (This is not currently implemented) 

# (4) Write map map_marks (if convergence)
# (This is implemented above)

# readdir(iter.folder_output, join = true)


#### End of code. 
###########################
###########################