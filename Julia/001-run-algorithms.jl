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
import CSV
import GeoArrays
import Parquet
import Random
using DataFrames
using Dates
using Distributions
using Patter

#### Load source files
include("./src/utils.jl")

#### Load datasets (map, iteration)
# Load map & iteration 
analysis  = "sim"
# analysis  = "real"
subanalysis = "main"
env       = GeoArrays.read(joinpath("data", "input", "map.tif"));
env_init  = Patter.rast(joinpath("data", "input", "map.tif"));
iteration = DataFrame(Arrow.Table(joinpath("data", "input", analysis, subanalysis, "iteration.feather")))

#### (optional) Use test settings
if false
  # Focus on a few iterations
  iteration = iteration[1:4, :]
  # Reduce batches & particle numbers
  iteration.n_batch .= 3
  iteration.n_particle_filter .= 20000
  iteration.n_particle_smoother .= 500
  # Clean up old files
  files = filter(f -> endswith(f, ".jld2"), readdir(ititerationer.folder_output; join=true))
  if length(files) > 0
    rm.(files; force=true)
  end
end

#### Select iteration 
# Set column types as needed
iteration.n_particle_filter = Int.(iteration.n_particle_filter);
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

#### Load timeline #
timeline = DataFrame(Arrow.Table(iter.file_timeline))
timeline.timestamp = DateTime.(timeline.timestamp)
timeline = timeline.timestamp
# timeline = CSV.read(iter.file_timeline, DataFrame, dateformat = "yyyy-mm-dd H:M:S")
# timeline.timestamp = DateTime.(timeline.timestamp)
# timeline = timeline.timestamp

#### Load acoustic observations & containers
# acoustics      = CSV.read(iter.file_acoustics, DataFrame, dateformat = "yyyy-mm-dd H:M:S")
# containers_fwd = CSV.read(iter.file_containers_fwd, DataFrame, dateformat = "yyyy-mm-dd H:M:S")
# containers_bwd = CSV.read(iter.file_containers_bwd, DataFrame, dateformat = "yyyy-mm-dd H:M:S")
acoustics = DataFrame(Arrow.Table(iter.file_acoustics))
containers_fwd = DataFrame(Arrow.Table(iter.file_containers_fwd))
containers_bwd = DataFrame(Arrow.Table(iter.file_containers_bwd))

#### Process columns
# Timestamps
acoustics.timestamp = DateTime.(acoustics.timestamp);
containers_fwd.timestamp = DateTime.(containers_fwd.timestamp);
containers_bwd.timestamp = DateTime.(containers_bwd.timestamp);
# Acoustic columns
acoustics.sensor_id = Int.(acoustics.sensor_id);
acoustics.obs = Int.(acoustics.obs);
acoustics.receiver_x = Float64.(acoustics.receiver_x);
acoustics.receiver_y = Float64.(acoustics.receiver_y);
acoustics.receiver_alpha = Float64.(acoustics.receiver_alpha);
acoustics.receiver_beta = Float64.(acoustics.receiver_beta);
acoustics.receiver_gamma = Float64.(acoustics.receiver_gamma);
# containers_fwd
containers_fwd
containers_fwd.obs = Int.(containers_fwd.obs);
containers_fwd.sensor_id = Int.(containers_fwd.sensor_id);
containers_fwd.centroid_x = Float64.(containers_fwd.centroid_x);
containers_fwd.centroid_y = Float64.(containers_fwd.centroid_y);
containers_fwd.radius = Float64.(containers_fwd.radius);
# containers_bwd
containers_bwd.obs = Int.(containers_bwd.obs);
containers_bwd.sensor_id = Int.(containers_bwd.sensor_id);
containers_bwd.centroid_x = Float64.(containers_bwd.centroid_x);
containers_bwd.centroid_y = Float64.(containers_bwd.centroid_y);
containers_bwd.radius = Float64.(containers_bwd.radius);
# # > This is necessary to avoid issues with data imported via Arrow
# > This is not necessary for data imported from CSV

#### Assemble datasets 
# Collate datasets & associated `ModelObs` instances into a typed dictionary 
datasets_fwd    = [acoustics, containers_fwd];
datasets_bwd    = [acoustics, containers_bwd];
model_obs_types = [ModelObsAcousticLogisTrunc, ModelObsContainer];
yobs_fwd        = assemble_yobs(datasets = datasets_fwd,
                                    model_obs_types = model_obs_types);
yobs_bwd        = assemble_yobs(datasets = datasets_bwd,
                                    model_obs_types = model_obs_types);


###########################
###########################
#### Run particle algorithms 

#### Set up algorithms 
# Define batches 
fwd_batches = [joinpath(iter.folder_output, "fwd-$i.jld2") for i in 1:iter.n_batch]
bwd_batches = [joinpath(iter.folder_output, "bwd-$i.jld2") for i in 1:iter.n_batch]
smo_batches = [joinpath(iter.folder_output, "smo-$i.jld2") for i in 1:iter.n_batch]
# Define output objects
fwd = bwd = smo = nothing 
# Define duration placeholders
td_fwd = td_bwd = td_smo = NaN

#### (1) Forward filter 

## Simulate initial states for the forward filter
xinit = simulate_states_init(map             = env_init, 
                             timeline        = timeline, 
                             state_type      = state,
                             xinit           = nothing, 
                             model_move      = model_move, 
                             datasets        = datasets_fwd,
                             model_obs_types = model_obs_types,
                             n_particle      = iter.n_particle_filter, 
                             direction       = "forward", 
                             output          = "Vector");

## Run the forward filter
t1_fwd = now()
fwd = particle_filter(timeline   = timeline,
                      xinit      = xinit,
                      yobs       = yobs_fwd,
                      model_move = model_move,
                      n_move     = 1,
                      n_record   = iter.n_particle_smoother,
                      direction  = "forward", 
                      batch      = fwd_batches,
                      progress   = Patter.progress_control(enabled = isinteractive()),
                      verbose    = isinteractive());
t2_fwd = now()
td_fwd = diffsecs(t2_fwd, t1_fwd)

## Collect outputs
# (To conserve disk space, we do not record states)
callstats            = fwd.callstats
diagnostics          = fwd.diagnostics 
diagnostics.routine .= callstats.routine

#### (2) Backward filter 
convergence = fwd.callstats.convergence[1] 
if convergence

  ## Simulate initial states for the backward filter
  xinit = simulate_states_init(map           = env_init, 
                             timeline        = timeline, 
                             state_type      = state,
                             xinit           = nothing, 
                             model_move      = model_move, 
                             datasets        = datasets_bwd,
                             model_obs_types = model_obs_types,
                             n_particle      = iter.n_particle_filter, 
                             direction       = "backward", 
                             output          = "Vector");

  ## Run the backward filter
  t1_bwd = now()
  bwd = particle_filter(timeline   = timeline,
                        xinit      = xinit,
                        yobs       = yobs_bwd,
                        model_move = model_move,
                        n_move     = 1,
                        n_record   = iter.n_particle_smoother,
                        direction  = "backward", 
                        batch      = bwd_batches,
                        progress   = Patter.progress_control(enabled = isinteractive()),
                        verbose    = isinteractive());
  t2_bwd = now()
  td_bwd = diffsecs(t2_bwd, t1_bwd)
  
  ## Collect outputs
  # (To conserve disk space, we do not record states)
  append!(callstats, bwd.callstats)
  bwd.diagnostics.routine .= bwd.callstats.routine
  append!(diagnostics, bwd.diagnostics)
  convergence = bwd.callstats.convergence[1]

end 

#### (3) Run smoother
# Set n_sim = 0 and cache = nothing for unrestricted models (n_move = 1)

if convergence

  ## Run smoother
  t1_smo = now()
  smo = particle_smoother_two_filter(timeline   = timeline,
                                     xfwd       = fwd_batches,
                                     xbwd       = bwd_batches,
                                     model_move = model_move,
                                     vmap       = nothing,
                                     n_particle = iter.n_particle_smoother,
                                     n_sim      = 0, 
                                     cache      = false, 
                                     batch      = smo_batches, 
                                     progress   = Patter.progress_control(enabled = isinteractive()), 
                                     verbose    = isinteractive());
  t2_smo = now()
  td_smo = diffsecs(t2_smo, t1_smo)

  ## Collate outputs
  callstats.n_iter .= Float64.(callstats.n_iter)
  append!(callstats, smo.callstats)
  smo.diagnostics.routine .= smo.callstats.routine
  append!(diagnostics, smo.diagnostics)

  ## Define smoothed DataFrame of states
  # * This requires bringing all smoothed files into memory
  # * This is not feasible
  # * Hence, this is implemented separately via Patter.pf_particles()
                                     
end 
                                   
#### (4) Record results
# We record all results for which the smoother was run 
# (i.e., for which the forward and backward filters converged)
if convergence

  # Write particles
  # (This is not currently implemented)
  # Parquet.write(iter_file_states, smo_df);
  
  # Write diagnostics 
  Arrow.write(iter.file_diagnostics, diagnostics; compress = Arrow.ZstdCompressor(level = 9))

  # Write callstats 
  Arrow.write(iter.file_callstats, callstats; compress = Arrow.ZstdCompressor(level = 9))

end 

#### (5) Cleanup batches
# We remove fwd_{i}.jld2 & bwd_{i}.jld2 files
# We only remove smo_{i}.jld2 after collating states (later)
foreach(f -> rm(f; force = true), fwd_batches)
foreach(f -> rm(f; force = true), bwd_batches)
readdir(iter.folder_output, join = true)


#### End of code. 
###########################
###########################