import Dates

# Time differences: compute the time difference (s)
# * `t2` and `t1` may be DateTime or Vector{DateTime}
function diffsecs(t2::Union{Dates.DateTime, Vector{Dates.DateTime}}, 
                  t1::Union{Dates.DateTime, Vector{Dates.DateTime}}) 
    Dates.value.(t2 .- t1) ./ 1000
end 