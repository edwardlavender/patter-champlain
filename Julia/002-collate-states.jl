###########################
###########################
#### collate-states.jl

#### Aims
# 1) Collate smoothed states

#### Prerequisites
# 1) This code should be run in parallel via bash on siam-linux20
#    (or on the machine where all particle outputs live)


###########################
###########################
#### Set up

#### Load base packages
import Base.Threads
import Pkg

#### Load local packages
Pkg.activate(".")
import CSV
using DataFrames
import Dates
import Patter 
import Parquet

#### Load datasets (map, iteration)
# Load map & iteration 
analysis  = "sim"
# analysis  = "real"
subanalysis = "main"
iteration = CSV.read(joinpath("data", "input", analysis, subanalysis, "iteration.csv"), DataFrame)

#### Select iteration 
row = parse(Int, ARGS[1])
iter = iteration[row, :]

#### Define local settings
if Threads.nthreads() != 1
  error("JULIA_NUM_THREADS must be 1 for parallelisation!")
end 


###########################
###########################
#### Collate states

# Define timeline 
timeline = CSV.read(iter.file_timeline,
                    DataFrame, 
                    dateformat = "yyyy-mm-dd H:M:S");

# Define smo-{i}.jld2 files
batch <- smo_batches = [joinpath(iter.folder_patter, "smo-{i}.jld2") for i in 1:iter.n_batch]

# Collate states Matrix in Julia 
smo_states = hcat([f["xsmo"] for f in map(jldopen, bbb)]...)

# Convert to DataFrame 
smo_states_df = Patter.r_get_states(smo_states, collect(1:length(timeline)), timeline);
  
# Collate 'pf_particles' object
# * This is implemented in R 

# Write output to file
Parquet.write(iter.file_states, smo_states_df)
  
# Cleanup batches
foreach(f -> rm(f; force = true), smo_batches)


#### End of code. 
###########################
###########################