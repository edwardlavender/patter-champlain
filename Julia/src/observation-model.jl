import GeoArrays
import Patter.distance

using Distributions
using Patter
using LogExpFunctions: logistic, log1pexp

# Define line of sight map
const los_map = GeoArrays.read(joinpath("data", "input", "map.tif"));

# Identify whether two points are connected by a direct line of sight (true/false)
# * Draw a linear line with 3 points between the two locations
# * Lookup middle value on the map
# * If NaN value, line of sight is 0.0; otherwise 1.0
function in_line_of_sight(env::GeoArrays.GeoArray, x0::Float64, y0::Float64, x1::Float64, y1::Float64)::Bool
    xm = (x0 + x1) * 0.5
    ym = (y0 + y1) * 0.5
    return !isnan(Patter.extract(env, xm, ym))
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
    dist = distance(state.x, state.y, model_obs.receiver_x, model_obs.receiver_y)
    if dist > model_obs.receiver_gamma
        return obs == 1 ? -Inf : 0.0
    end
    if !in_line_of_sight(env, state.x, state.y, model_obs.receiver_x, model_obs.receiver_y)
        return obs == 1 ? -Inf : 0.0
    end 
    η = model_obs.receiver_alpha + model_obs.receiver_beta * dist
    return obs == 1 ? -log1pexp(-η) : -log1pexp(η)
end

# Simulate method 
function Patter.simulate_obs(state::State, model_obs::ModelObsAcousticLogisTruncLos, t::Int64, env::GeoArrays.GeoArray = los_map)
    # Compute probability of a detection from location (state.x, state.y)
    obs = 1
    prob = exp(Patter.logpdf_obs(state, model_obs, t, obs, env))
    # Simulate detection/non-detection
    return rand(Bernoulli(prob)) + 0
end