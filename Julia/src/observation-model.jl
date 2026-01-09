import GeoArrays
using Distributions
using Patter
import Patter.distance
using LogExpFunctions: logistic, log1pexp

# Define line of sight map
los_map = GeoArrays.read(joinpath("data", "input", "map.tif"));

# Identify whether two points are connected by a direct line of sight (true/false)
# * Draw a linear line with n points between the two locations
# * Lookup values on the map
# * If any NaN values, line of sight is 0.0; otherwise 1.0
function in_line_of_sight(env::GeoArrays.GeoArray, x0::Real, y0::Real, x1::Real, y1::Real, n::Int=5)
    xs = range(x0, x1, length=n)
    ys = range(y0, y1, length=n)
    !any(isnan.([Patter.extract(env, xs[i], ys[i]) for i in eachindex(xs)])) + 0.0
end

# Define ModelObsAcousticLogisTrunc accounting for line of sight (Los) 
struct ModelObsAcousticLogisTruncLos <: Patter.ModelObs
    sensor_id::Int64
    # Receiver coordinates
    receiver_x::Float64
    receiver_y::Float64
    # Detection probability parameters
    receiver_alpha::Float64
    receiver_beta::Float64
    receiver_gamma::Float64
end

# Log probability method
function Patter.logpdf_obs(state::State, model_obs::ModelObsAcousticLogisTruncLos, t::Int64, obs::Int64, env::GeoArrays.GeoArray = los_map)

    # Compute distance between the particle and receiver
    dist = distance(state.x, state.y, model_obs.receiver_x, model_obs.receiver_y)

    # Evaluate line of sight, if needed
    los = 1.0
    if dist < model_obs.receiver_gamma
        los = in_line_of_sight(env, state.x, state.y, model_obs.receiver_x, model_obs.receiver_y)
    end 

    # Compute log probability 
    η = model_obs.receiver_alpha + model_obs.receiver_beta * dist
    if obs == 1
        if dist > model_obs.receiver_gamma || los == 0.0
            return -Inf
        else
            return -log1pexp(-η)
        end
    elseif obs == 0
        if dist > model_obs.receiver_gamma || los == 0.0
            return 0.0
        else
            return -log1pexp(η)
        end
    else
        error("Acoustic observations should be coded as 0 or 1.")
    end

end

# Simulate method 
function Patter.simulate_obs(state::State, model_obs::ModelObsAcousticLogisTruncLos, t::Int64, env = los_map)
    prob = exp(logpdf_obs(state::State, model_obs::ModelObsAcousticLogisTruncLos, t::Int64, obs::Int64, env))
    rand(Bernoulli(prob)) + 0
end