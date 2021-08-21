# # Standard utility

# intertemporal aggregator
ϕ(::Log, c, v, β) = c^(1 - β) * v^β
ϕ(u::Power, c, v, β) = ((1 - β) * c^u.m + β * v^u.m)^u.inv
ϕ(u, c, v, β) = ϕ(intertemporal_aggregator(u), c, v, β)


# certainty equivalent aggregator
function ce(::Log, p⃗, v⃗)
    s = zero(eltype(p⃗))
    @turbo for i in eachindex(p⃗, v⃗)
        s += p⃗[i] * log(v⃗[i])
    end
    return exp(s)
end

function ce(u::Power, p⃗, v⃗)
    s = zero(eltype(p⃗))
    @turbo for i in eachindex(p⃗, v⃗)
        s += p⃗[i] * (v⃗[i] ^ u.m)
    end
    return s ^ u.inv
end

ce(u, p⃗, v⃗) = ce(ce_aggregator(u), p⃗, v⃗)


# # Disutility of labor methods

# Fixed Labor
labor(v::FixedLabor, w) = v.n  # labor supply
labor_income(v::FixedLabor, w) = w * v.n  # labor income
disutility(::FixedLabor{R}, n) where {R} = zero(R)  # disutility
disutility_given_w(::FixedLabor{R}, w) where {R} = zero(R)  # disutility


# GHH
labor(v::GHH, w) = (w / v.θ)^v.ν  # labor supply
labor_income(v::GHH, w) = w^(1 + v.ν) / (v.θ^v.ν)  # labor income
disutility(v::GHH, n) = v.θ * n^(1 + 1/v.ν) / (1 + 1/v.ν)
disutility_given_w(v::GHH, w) = v.θ * (w / v.θ)^(v.ν + 1) / (1 + 1/v.ν)


# Household consumption given w and other income
consumption(h, w, other_net_income) = labor_income(get_v(h), w) + other_net_income

# GHH inner value
net_consumption(h, w, other_net_income) = consumption(h, w, other_net_income) - disutility_given_w(get_v(h), w)
