###########################
###########################
#### run-algorithms.jl

#### Aims
# 1) Test selected routines 

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
import GeoArrays
import Patter
import Patter: distance
using LogExpFunctions: logistic, log1pexp
using Plots
using Test
include("../src/observation-model.jl")

#### Load data 
los_map = GeoArrays.read(joinpath("data", "input", "map.tif"));


###########################
###########################
#### Test line of sight 

# Test two points that are not in line of sight
x0 = 1786904
y0 = 948148
x1 = 1780253
y1 = 946763.2
Plots.plot(los_map)
scatter!([x0], [y0])
scatter!([x1], [y1])
@test in_line_of_sight(los_map, x0, y0, x1, y1) == 0.0

# Test two points that are in line of sight 
x0 = 1784188
y0 = 933737.8
x1 = 1787791
y1 = 929026.5
Plots.plot(env)
scatter!([x0], [y0])
scatter!([x1], [y1])
@test in_line_of_sight(los_map, x0, y0, x1, y1) == 1.0


###########################
###########################
#### Test logpdf_obs() logic 

# Define _logpdf_obs internal function 
# * At the moment of detection, individuals must be (a) within detection range and (b) line of sight
# * At the moment of non detection, individuals are more likely to be beyond detection range and/or not in line of sight (but may be within both)
function _logpdf_obs(obs, model_obs, dist::Real, los::Real)
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

#### Example
# Define example distance and model_obs object
dist = 50
model_obs = (receiver_alpha=1.205216, receiver_beta=-0.001672085, receiver_gamma=7000)
# (1) If obs == 1, probability @ 50 metres if line of sight should be close to 1
exp(_logpdf_obs(1.0, model_obs, dist, 1.0))
exp(-log1pexp(-(model_obs.receiver_alpha + model_obs.receiver_beta * dist)))
# (2) If obs == 1, probability @ 50 metres if _not_ line of sight (should be 0)
exp(_logpdf_obs(1.0, model_obs, dist, 0.0))

# (3) obs == 0, probability @ 50 m if line of sight should be low
exp(_logpdf_obs(0.0, model_obs, dist, 1.0))
exp(-log1pexp(model_obs.receiver_alpha + model_obs.receiver_beta * dist))
# (4) obs == 0, probability @ 50 m if no line of sight should be high (1.0) 
exp(_logpdf_obs(0.0, model_obs, dist, 0.0))



#### End of code. 
###########################
###########################