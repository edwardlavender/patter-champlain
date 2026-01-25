import GeoArrays
import Patter.distance

using Distributions
using Patter
using LogExpFunctions: logistic, log1pexp

# Define line of sight map
los_map = GeoArrays.read(joinpath("data", "input", "map.tif"));

# Identify whether two points are connected by a direct line of sight (true/false)
# * Draw a linear line with 3 points between the two locations
# * Lookup middle value on the map
# * If NaN value, line of sight is 0.0; otherwise 1.0
function in_line_of_sight(env::GeoArrays.GeoArray, x0::Float64, y0::Float64, x1::Float64, y1::Float64)::Bool
    xm = (x0 + x1) * 0.5
    ym = (y0 + y1) * 0.5
    !isnan(Patter.extract(env, xm, ym))
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
    # (Implementing this when obs == 1 does not seem to boost speed)
    if dist <= model_obs.receiver_gamma
        los = in_line_of_sight(env, state.x, state.y, model_obs.receiver_x, model_obs.receiver_y)
    end 

    # Compute probability of detection
    # - Allow the model to 'flicker' between the standard model and steeper model with low probability 
    # - This model recognises that there are moments in time when detection probability is much lower
    # - The properties of the steeper model are currently hard-coded (F146)
    # - TO DO Review this in due course & align with sensitivity analysis, if needed
    flicker = ifelse(rand() < 0.95, false, true)
    if (flicker)
         η = 0.635159329 + -0.002724118 * dist
    else 
        η = model_obs.receiver_alpha + model_obs.receiver_beta * dist
    end 

    # Compute log probability 
    # η = model_obs.receiver_alpha + model_obs.receiver_beta * dist
    if obs == 1
        if dist > model_obs.receiver_gamma  || los == false
            return -Inf
        else 
            return -log1pexp(-η)
        end 
    elseif obs == 0
        if dist > model_obs.receiver_gamma || los == false
            return 0.0
        else
            return -log1pexp(η)
        end
    end 

end

# Simulate method 
function Patter.simulate_obs(state::State, model_obs::ModelObsAcousticLogisTruncLos, t::Int64, env = los_map)
    # Compute probability of a detection from location (state.x, state.y)
    obs = 1
    prob = exp(Patter.logpdf_obs(state, model_obs, t, obs, env))
    # Simulate detection/non-detection
    rand(Bernoulli(prob)) + 0
end