# # Standard utility 

# Assumes time separability where c is consumption net of the disutility_given_w 
# of labor

outside(u::AbstractUtility, c, v, β) = u(c) + β * v
inside(::AbstractUtility, x) = x
inverse_inside(::AbstractUtility, x) = x 
ss_value(u::AbstractUtility, c, β) = u(c) / (1 - β) # stationary value of constant consumption

# # EZ Methods -- non separability

function outside(f::EZ{<:Power, T1, T2}, c, v, β) where {T1, T2}
    ((1 - β) * c^f.ies.m + β * v^f.ies.m)^f.ies.n
end 

function outside(::EZ{<:Log, T1, T2}, c, v, β) where {T1, T2}
    c^(1 - β) * v^β 
end 

inside(::EZ{T1, <:Log, T2}, x) where {T1, T2} = log(x) 
inside(f::EZ{T1, <:Power, T2}, x) where {T1, T2} = x^f.ra.m

inverse_inside(::EZ{T1, <:Log, T2}, x) where {T1, T2} = exp(x) 
inverse_inside(f::EZ{T1, <:Power, T2}, x) where {T1, T2} = x^f.ra.n

ss_value(::EZ, c, β) = c  # stationary value with constant consumption 

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