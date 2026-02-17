import Arrow
import GeoArrays
import JLD2
using DataFrames
using DataFramesMeta
using Dates

# Compute time difference (s) between t2 and t1 
# * `t2` and `t1` may be DateTime or Vector{DateTime}
# function diffsecs(t2::Union{Dates.DateTime,Vector{Dates.DateTime}},
#     t1::Union{Dates.DateTime,Vector{Dates.DateTime}})
#     Dates.value.(t2 .- t1) ./ 1000
# end

# Compute the 'map_marks' DataFrame
# > Given the Map (env) GeoArray and a DataFrame of sampled particles
#   this function computes the total weight in each grid cell per time step
function map_marks(env::GeoArrays.GeoArray, coord::DataFrame)

    # Define grid cells and coordinates
    coord.id = [GeoArrays.indices(env, [coord.x[i], coord.y[i]]) for i in 1:nrow(coord)]
    xy       = GeoArrays.coords.(Ref(env), coord.id)
    coord.x  = getindex.(xy, 1)
    coord.y  = getindex.(xy, 2)
    @select!(coord, :timestep, :id, :x, :y)

    # Compute grid cell weights by timestep
    coord = @chain coord begin
        @groupby(:timestep)
        @transform(:mark = 1 / length(:timestep))
        # For each time step, compute the total weight of each grid cell 
        @groupby(:timestep, :x, :y)
        @combine(:mark = sum(:mark))
    end

    coord

end

# Compute map grid cell area (m^2) on UTM grid 
# import LinearAlgebra
# A = abs(LinearAlgebra.det(env.f.linear))
# A = prod(abs.(step.(dims(env_init))))

# Compute the area (ncells) spanned by 50 and 95 % of the probability mass by time step
function map_uncertainty(coord::DataFrame)

    # By timestep, compute the minimum number of cells containing 50 or 95 % of the probability mass
    ncells = DataFrames.combine(groupby(coord, :timestep)) do grp
        cmark = cumsum(sort(grp.mark, rev=true))
        (; ncell_core = findfirst(>=(0.5), cmark),
            ncell_home = findfirst(>=(0.95), cmark))
    end

    return ncells

end

# Collate smoothed states for a specified iteration over all batches
# > Smoothed states are written to file as iter.file_states
function write_file_states(iter)
    
    # Define timeline 
    timeline           = DataFrame(Arrow.Table(iter.file_timeline))
    timeline.timestamp = DateTime.(timeline.timestamp)
    timeline           = timeline.timestamp

    # Define smo-{i}.jld2 files
    # > This returns String[] if no files exist 
    smo_batches = [joinpath(iter.folder_output, "smo-$i.jld2") for i in 1:iter.n_batch]
    smo_batches = smo_batches[isfile.(smo_batches)]

    # Collate states, if smo-$i.jld2 files have been produced 
    # (i.e., if the two filters converged)
    if length(smo_batches) > 0
        # Collate states Matrix in Julia 
        smo_states = hcat([f["xsmo"] for f in map(JLD2.jldopen, smo_batches)]...)
        # Convert to DataFrame 
        smo_states_df = Patter.r_get_states(smo_states, collect(1:length(timeline)), timeline)
        # Write output to file
        Arrow.write(iter.file_states, smo_states_df; compress = Arrow.ZstdCompressor(level = 9))
    end

    # Cleanup batches
    # foreach(f -> rm(f; force=true), smo_batches)

    return nothing
end 