# # Standard utility 

# intertemporal aggregator
ϕ(::Log, c, v, β) = c^(1 - β) * v^β 
ϕ(u::Power, c, v, β) = ((1 - β) * c^u.m + β * v^u.m)^u.inv
ϕ(u, c, v, β) = ϕ(get_ies(u), c, v, β) 


# certainty equivalent aggregator 
function ce(::Log, π̄, v̄) 
    s = zero(eltype(π̄))
    @turbo for i in eachindex(π̄, v̄)
        s += π̄[i] * log(v̄[i])
    end 
    return exp(s)
end 

function ce(u::Power, π̄, v̄) 
    s = zero(eltype(π̄))
    @turbo for i in eachindex(π̄, v̄)
        s += π̄[i] * (v̄[i] ^ u.m)
    end 
    return s ^ u.inv
end 

ce(u, π̄, v̄) = ce(get_ra(u), π̄, v̄) 


# # Disutility of labor methods

# Fixed Labor 
n_of_w(v::FixedLabor, w) = v.n  # labor supply 
disutility(::FixedLabor{R}, n) where {R} = zero(R)  # disutility
disutility_given_w(::FixedLabor{R}, w) where {R} = zero(R)  # disutility
labor_income_given_w(v::FixedLabor, w) = w * v.n  # labor income


# GHH 
n_of_w(v::GHH, w) = (w / v.θ)^v.ν  # labor supply 
disutility(v::GHH, n) = v.θ * n^(1 + 1/v.ν) / (1 + 1/v.ν) 
disutility_given_w(v::GHH, w) = v.θ * (w / v.θ)^(v.ν + 1) / (1 + 1/v.ν)  
labor_income_given_w(v::GHH, w) = w^(1 + v.ν) / (v.θ^v.ν)  # labor income


# Consumption from budget constraint
get_c(h, w, other_net_income) = labor_income_given_w(get_v(h), w) + other_net_income

# GHH inner value 
get_x(h, w, other_net_income) = get_c(h, w, other_net_income) - disutility_given_w(get_v(h), w)