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
using BenchmarkTools
# using Revise
import Pkg

#### Load local packages
Pkg.activate(".")
import CSV
using DataFrames
import Dates
import GeoArrays
import Parquet
using Patter
import Random 

#### Load source files
include("./src/utils.jl")

#### Load datasets (map, iteration)
# Load map & iteration 
analysis  = "sim"
# analysis  = "real"
subanalysis = "main"
env       = GeoArrays.read(joinpath("data", "input", "map.tif"));
iteration = CSV.read(joinpath("data", "input", analysis, subanalysis, "iteration.csv"), DataFrame)

#### Select iteration 
row = parse(Int, ARGS[1])
iter = iteration[row, :]

#### Define local settings
# Check JULIA_NUM_THREADS = 1
if Threads.nthreads() != 1
  error("JULIA_NUM_THREADS must be 1 for parallelisation!")
end 
# Set seed
Random.seed!(123);
start = now();


###########################
###########################
#### Prepare data 

#### Examine iterations
for (n,t) in zip(names(iter), eltype.(eachcol(iter)))
    println(n, ": ", t)
end

#### Process column types 
# iter.D         = Vector(iter.D);
# iter.D         = Float64.(iter.D);
iter.sensitivity = string.(iter.sensitivity);


###########################
###########################
#### Define state-space model 


###########################
#### Define movement model

#### Define state 
state = StateXY

#### Define movement model
model_move = ModelMoveCXY(env, 
                          iter.mobility, 
                          truncated(Gamma(iter.k, iter.theta), upper = iter.mobility), 
                          MixtureModel([truncated(Normal(0.0, iter.sigma), -pi, pi), Uniform(-pi, pi)], [0.99, 0.01]))


###########################
#### Define observation model 

#### Load timeline 
timeline = CSV.read(iter.file_timeline,
                    DataFrame, 
                    dateformat = "yyyy-mm-dd H:M:S");

#### Load acoustic observations
# Acoustics 
acoustics = CSV.read(joinpath("data", "input", analysis, subanalysis, "acoustics.csv"), DataFrame);
first(acoustics, 6)
# Forward containers 
containers_fwd = CSV.read(joinpath("data", "input", analysis, subanalysis, "containers-fwd.csv"), DataFrame);
first(containers_fwd, 6)
# Backward containers 
containers_bwd = CSV.read(joinpath("data", "input", analysis, subanalysis, "containers-bwd.csv"), DataFrame);
first(containers_bwd, 6)

#### Process time stamps
acoustics.timestamp      = DateTime.(acoustics.timestamp, "yyyy-mm-dd HH:MM:SS");
containers_fwd.timestamp = DateTime.(containers_fwd.timestamp, "yyyy-mm-dd HH:MM:SS");
containers_bwd.timestamp = DateTime.(containers_bwd.timestamp, "yyyy-mm-dd HH:MM:SS");
archival.timestamp       = DateTime.(archival.timestamp, "yyyy-mm-dd HH:MM:SS");

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
fwd_batches = [joinpath(iter.folder_patter, "fwd-{i}.jld2") for i in 1:iter.n_batch]
bwd_batches = [joinpath(iter.folder_patter, "bwd-{i}.jld2") for i in 1:iter.n_batch]
smo_batches = [joinpath(iter.folder_patter, "smo-{i}.jld2") for i in 1:iter.n_batch]
# Define output objects
fwd = bwd = smo = nothing 
# Define duration placeholders
td_fwd = td_bwd = td_smo = NaN

#### (1) Forward filter 

## Simulate initial states for the forward filter
xinit = simulate_states_init(map             = env, 
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
                      n_record   = iter.n_particle_smoother
                      direction  = "forward", 
                      batch      = fwd_batches,
                      progress   = (),
                      verbose    = false);
t2_fwd = now()
td_fwd = diffsecs(t2_fwd, t1_fwd)

## Collect outputs
# (To conserve disk space, we do not record states)
callstats            = fwd.callstats
diagnostics          = fwd.diagnostics 
diagnostics.routine .= callstats.routine

#### (2) Backward filter 
convergence = fwd.diagnostics.convergence 
if convergence

  ## Simulate initial states for the backward filter
  xinit = simulate_states_init(map             = env, 
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
  bwd_batches = [joinpath("tmp", "fwd-{i}.jld2") for i in 1:iter.n_batch]
  bwd = particle_filter(timeline   = timeline,
                        xinit      = xinit,
                        yobs       = yobs_bwd,
                        model_move = model_move,
                        n_move     = 1,
                        n_record   = iter.n_particle_smoother
                        direction  = "backward", 
                        batch      = bwd_batches,
                        progress   = (),
                        verbose    = false);
  t2_bwd = now()
  td_bwd = diffsecs(t2_bwd, t1_bwd)
  
  ## Collect outputs
  # (To conserve disk space, we do not record states)
  append!(callstats, bwd.callstats)
  bwd.diagnostics.routine .= bwd.callstats.routine
  append!(diagnostics, bwd.diagnostics)
  
  end 

#### (3) Run smoother
# Set n_sim = 0 and cache = nothing for unrestricted models (n_move = 1)
convergence = bwd.diagnostics.convergence 
if convergence

  ## Run smoother
  t1_smo = now()
  smo = particle_smoother_two_filter(timeline   = timeline,
                                     xfwd       = fwd_batches,
                                     xbwd       = bwd.batches,
                                     model_move = model_move,
                                     vmap       = nothing,
                                     n_particle = iter.n_particle_smoother,
                                     n_sim      = 0, 
                                     cache      = nothing, 
                                     batch      = smo_batches, 
                                     progress   = (), 
                                     verbose    = false);
  t2_smo = now()
  td_smo = diffsecs(t2_smo, t1_smo)

  ## Collate outputs
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
if convergence

  # Write particles
  # (This is not currently implemented)
  # Parquet.write(iter_file_states, smo_df);
  
  # Write diagnostics 
  Parquet.write(iter.file_diagnostics, diagnostics);
  
  # Write callstats
  # * key = {individual_id}-{time_id}-{sensitivity}
  # * time_fwd, time_bwd and time_smo collectively define computation time 
  # - If NaN -> convergence failure 
  # - In R, we select rows without convergence failures that also successed with smoothing 
  callstats = DataFrame(key       = iter.key, 
                        timestamp = start, 
                        time_fwd  = td_fwd
                        time_bwd  = td_bwd
                        time_smo  = td_smo, 
                        # Add smoothing callstats
                        
                        );
  CSV.write(iter.file_callstats, callstats);

end 

#### (5) Cleanup batches
# We remove fwd_{i}.jld2 & bwd_{i}.jld2 files
# We only remove smo_{i}.jld2 after collating states (later)
foreach(f -> rm(f; force = true), fwd_batches)
foreach(f -> rm(f; force = true), bwd_batches)


#### End of code. 
###########################
###########################