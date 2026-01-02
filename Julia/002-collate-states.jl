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
import Arrow
import CSV
import Patter 
using DataFrames
using Dates
using JLD2

#### Load datasets (map, iteration)
# Load map & iteration 
analysis  = "sim"
# analysis  = "real"
subanalysis = "main"
iteration = DataFrame(Arrow.Table(joinpath("data", "input", analysis, subanalysis, "iteration.feather")))

#### Select iteration 
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
Threads.nthreads()
if (Threads.nthreads() != 1) && !isinteractive()
  error("JULIA_NUM_THREADS must be 1 for parallelisation!")
end


###########################
###########################
#### Collate states

# Define timeline 
timeline = DataFrame(Arrow.Table(iter.file_timeline))
timeline.timestamp = DateTime.(timeline.timestamp)
timeline = timeline.timestamp

# Define smo-{i}.jld2 files
# > This returns String[] if no files exist 
smo_batches = [joinpath(iter.folder_output, "smod-$i.jld2") for i in 1:iter.n_batch]
smo_batches = smo_batches[isfile.(smo_batches)]

# Collate states, if smo-$i.jld2 files have been produced 
# (i.e., if the two filters converged)
if length(smo_batches) > 0

  # Collate states Matrix in Julia 
  smo_states = hcat([f["xsmo"] for f in map(jldopen, smo_batches)]...)

  # Convert to DataFrame 
  smo_states_df = Patter.r_get_states(smo_states, collect(1:length(timeline)), timeline);
  
  # Write output to file
  Arrow.write(iter.file_states, smo_states_df; compress = Arrow.ZstdCompressor(level = 9))

end 
  
# Cleanup batches
foreach(f -> rm(f; force = true), smo_batches)
# readdir(iter.folder_output, join = true)


#### End of code. 
###########################
###########################