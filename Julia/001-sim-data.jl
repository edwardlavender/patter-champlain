###########################
###########################
#### sim-data.jl

#### Aims
# 1) Simulate data

#### Prerequisites
# 1) NA


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
include("./src/observation-model.jl")

#### Load datasets
env  = GeoArrays.read(joinpath("data", "input", "map.tif"));
iter = DataFrame(Arrow.Table(joinpath("data", "input", "sim", "main", "iter.feather")));

#### Define local settings
Random.seed!(123);


###########################
###########################
#### Run simulation

#### Define iteration settings
iter = iter[1, :]

#### Define timeline 
timeline           = DataFrame(Arrow.Table(iter.file_timeline))
timeline.timestamp = DateTime.(timeline.timestamp)
timeline           = timeline.timestamp

#### Define starting points
xinit = DataFrame(Arrow.Table(iter.file_xinit));
xinit = [StateCXY(xinit.map_value[i], xinit.x[i], xinit.y[i], xinit.heading[i]) for i in 1:nrow(xinit)];

#### Define movement model
model_move = ModelMoveCXY(env, 
                          iter.mobility, 
                          truncated(Gamma(iter.shape, iter.scale), upper = iter.mobility), 
                          MixtureModel([truncated(Normal(0.0, iter.phi), -pi, pi), Uniform(-pi, pi)], [0.99, 0.01]))

#### Simulate movement paths (~3 s)
@time paths = simulate_path_walk(xinit = xinit, model_move = model_move, timeline = timeline);
@time paths_df = Patter.r_get_states(paths, collect(1:length(timeline)), timeline)

#### Simulate observations (~22 s, 500 s)
# Define model_obs 
moorings = DataFrame(Arrow.Table(iter.file_moorings));
model_obs = Patter.julia_get_model_obs([moorings], [ModelObsAcousticLogisTruncLos]);
# Simulate acoustic observations 
# * Acoustics is a Dict with one element for each individual (path)
# * Each individual's element is a Dict with one element for each time step
# * Each timestep is a Vector with one element for each receiver 
@time acoustics = simulate_yobs(paths = paths, model_obs = model_obs, timeline = timeline);
# Count the number of detections per individual
n_detections = DataFrame((individual = id,
                         count = sum(receiver[1]
                            for timesteps in values(individuals)
                            for receiver  in timesteps))
                            for (id, individuals) in acoustics)
sort!(n_detections, :count)
# Build DataFrame 
@time acoustics_list = Patter.r_get_dataset(acoustics, ModelObsAcousticLogisTruncLos)
for i in eachindex(acoustics_list)
    acoustics_list[i].path_id .= i
end 
acoustics_df = vcat(acoustics_list...)

#### Write to file (~6 s)
Arrow.write(iter.file_paths_raw, paths_df; compress = Arrow.ZstdCompressor(level = 9))
Arrow.write(iter.file_acoustics_raw, acoustics_df; compress = Arrow.ZstdCompressor(level = 9))


#### End of code. 
###########################
###########################