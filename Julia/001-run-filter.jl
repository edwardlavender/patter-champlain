###########################
###########################
#### run-algorithms.jl

#### Aims
# 1) Run the forward filter
#    - This code simply runs the particle filter forwards
#    - This is a stripped down version of 002-run-algorithms.R
#    - It is used to trial settings & quickly identify convergence issues
#      for diagnosis (via trial-patter.R) and 

#### Prerequisites
# 1) This code can be run interactively or via bash from SIA-LAVENDED or siam-linux20
#    (It is mainly designed for parallel implementation via bash on siam-linux20)
# 2) Each run of this script requires ~0.3 % memory on siam-linux20
#    So we can comfortably use 0.3 * 100 CPUs 


###########################
###########################
#### Set up

#### Load base packages
import Base.Threads
import Pkg

#### Load local packages
Pkg.activate(".")
import Arrow
import GeoArrays
import Random
using DataFrames
using Dates
using Distributions
using Patter

#### Load source files
include("./src/utils.jl")

#### Load datasets (map, iteration)
# Load map & iteration 
# analysis  = "sim"
analysis  = "real"
subanalysis = "main"
env       = GeoArrays.read(joinpath("data", "input", "map.tif"));
env_init  = Patter.rast(joinpath("data", "input", "map.tif"));
iteration = DataFrame(Arrow.Table(joinpath("data", "input", analysis, subanalysis, "iteration-1.feather")))

#### Select iteration 
# Set column types as needed
iteration.n_particle_filter = Int.(iteration.n_particle_filter);
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

#### Load acoustic observations & containers
acoustics      = DataFrame(Arrow.Table(iter.file_acoustics))
containers_fwd = DataFrame(Arrow.Table(iter.file_containers_fwd))

#### Process columns
# Timestamps
acoustics.timestamp      = DateTime.(acoustics.timestamp);
containers_fwd.timestamp = DateTime.(containers_fwd.timestamp);
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

#### Assemble datasets 
# Collate datasets & associated `ModelObs` instances into a typed dictionary 
datasets_fwd    = [acoustics, containers_fwd];
model_obs_types = [ModelObsAcousticLogisTrunc, ModelObsContainer];
yobs_fwd        = assemble_yobs(datasets = datasets_fwd,
                                model_obs_types = model_obs_types);


###########################
###########################
#### Run filter 

#### Simulate initial states for the forward filter
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

#### Run the forward filter
# We are only interested in convergence success/failure, 
# so we only record 1 particle in memory at each time step (no batches)
# For the subset of convergence failures, we will dig into the particle dynamics
# in a separate script (we can't do this for all runs as it requires
# too much disk space). 
fwd = particle_filter(timeline   = timeline,
                      xinit      = xinit,
                      yobs       = yobs_fwd,
                      model_move = model_move,
                      n_move     = 1,
                      n_record   = 1,
                      direction  = "forward", 
                      progress   = Patter.progress_control(enabled = isinteractive()),
                      verbose    = isinteractive());

#### Write to file 
Arrow.write(iter.file_callstats_filter, fwd.callstats; compress = Arrow.ZstdCompressor(level = 9))


# readdir(iter.folder_output, join = true)


#### End of code. 
###########################
###########################