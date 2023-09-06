

# # Disutility of labor methods

# Fixed Labor
labor(v::FixedLabor; w) = v.n  # labor supply
labor_income(v::FixedLabor; w) = w * v.n  # labor income
disutility(::FixedLabor{R}; n) where {R} = zero(R)  # disutility
disutility_given_w(::FixedLabor{R}; w) where {R} = zero(R)  # disutility


# GHH
labor(v::GHH; w) = (w / v.θ)^v.ν  # labor supply
labor_income(v::GHH; w) = w^(1 + v.ν) / (v.θ^v.ν)  # labor income
disutility(v::GHH; n) = v.θ * n^(1 + 1/v.ν) / (1 + 1/v.ν)
disutility_given_w(v::GHH; w) = v.θ * (w / v.θ)^(v.ν + 1) / (1 + 1/v.ν)

