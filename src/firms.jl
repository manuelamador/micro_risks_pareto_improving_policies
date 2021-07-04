################################################
# Technology Methods
###############################################

# These operate in efficiency units

(f::CobbDouglas)(k) = k^(f.α)
_mpn_from_mpk(f::CobbDouglas, mpk) = (1 - f.α) * (f.α / mpk)^(f.α / (1 - f.α))
_k_from_mpk(f::CobbDouglas, mpk) = (f.α / mpk)^(1 / (1 - f.α))
_df_of_k(f::CobbDouglas, k) = f.α * k ^ (f.α - 1)

(f::CES)(k) = (f.α * k^((f.ρ - 1) / f.ρ) + 1 - f.α)^(f.ρ / (f.ρ - 1))
_df_of_k(f::CES, k) = f.α * (f.α + (1 - f.α) * k^((1 - f.ρ)/f.ρ) )^(1/ (f.ρ - 1))
_k_from_mpk(f::CES, mpk) = (((mpk / f.α)^(f.ρ - 1) - f.α) / (1 - f.α))^(f.ρ / (1 - f.ρ))
_mpn_from_mpk(f::CES, mpk) = (1 - f.α)^(f.ρ / (f.ρ - 1)) * (1 - mpk^(1 - f.ρ)* f.α^f.ρ)^(-1 / (f.ρ - 1))


################################################
# General Technology Methods
###############################################

rK_from_r(; t, r) = r + get_δ(t)
mpn_from_factors(t; k, n) = mpn_from_mpk(t, mpk_from_factors(t; k, n))

# TODO: make this technology specific using closed forms.
function golden_rule_k(t, n0; xrange = (0.0, 100.0))
    return find_zero(k -> mpk_from_factors(t; k = k, n = n0) - t.δ, xrange)
end

# ## Specialized methods

# No Markups

output(t::Technology; k, n) = t.f(k / (t.A * n)) * t.A * n
mpk_from_factors(t::Technology; k, n) = _df_of_k(t.f, k / (t.A * n))
mpk_from_after_tax_rK(::Technology, rK) =  rK
rL_from_mpk(t::Technology, mpk) = mpn_from_mpk(t, mpk)
mpn_from_mpk(t::Technology, mpk) = t.A * _mpn_from_mpk(t.f, mpk)
k_from_mpk(t::Technology; mpk, n) = t.A * n * _k_from_mpk(t.f, mpk)

# Markups

output(t::MarkupTechnology; k, n) = (1 - t.m) * t.f(k / (t.A * n)) * t.A * n - t.X
mpk_from_factors(t::MarkupTechnology; k, n) = (1 - t.m) * _df_of_k(t.f, k / (t.A * n))
mpk_from_after_tax_rK(t::MarkupTechnology, rK) = t.μ * rK
rL_from_mpk(t::MarkupTechnology, mpk) = mpn_from_mpk(t, mpk) / t.μ
mpn_from_mpk(t::MarkupTechnology, mpk) = t.A * (1 - t.m) * _mpn_from_mpk(t.f, mpk / (1 - t.m))
k_from_mpk(t::MarkupTechnology; mpk, n) = t.A * n * _k_from_mpk(t.f, mpk / (1 - t.m))
