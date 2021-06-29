################################################
# Technology Methods
###############################################

# These operate in efficiency units

(f::CobbDouglas)(k) = k^(f.α)

function _mpn_from_mpk(f::CobbDouglas, mpk) 
    @unpack α = f
    return (1 - α) * (α / mpk)^(α / (1 - α))
end 

function _k_from_mpk(f::CobbDouglas, mpk) 
    @unpack α = f
    return (α / mpk)^(1 / (1 - α))
end

function _df_of_k(f::CobbDouglas, k)
    @unpack α = f
    return α * k ^ (α - 1)
end 

function (f::CES)(k) 
    @unpack α, ρ = f
    return (α * k^((ρ - 1) / ρ) + 1 - α)^(ρ / (ρ - 1))
end

function _df_of_k(f::CES, k)
    @unpack α, ρ = f
    # return ((α * k^((ρ - 1) / ρ) + 1 - α)^(ρ / (ρ - 1) - 1)) * α * k^((ρ - 1) / ρ - 1)
    return α * (α + (1 - α) * k^((1 - ρ)/ρ) )^(1/ (ρ - 1))  
end 

function _k_from_mpk(f::CES, mpk) 
    @unpack α, ρ = f
    return (((mpk / α)^(ρ - 1) - α) / (1 - α))^(ρ / (1 - ρ))
end 

function _mpn_from_mpk(f::CES, mpk) 
    @unpack α, ρ = f
   return (1 - α)^(ρ / (ρ - 1)) * (1 - mpk^(1 - ρ)* α^ρ)^(-1 / (ρ - 1))
end


################################################
# General Technology Methods
###############################################

get_δ(t::AbstractTechnology) = t.δ
get_y(t, factors) = get_y(t; factors.k, factors.n)
rK_from_r(; t, r) = r + get_δ(t)


# TODO: make this technology specific using closed forms. 
function get_golden_k(t, n0; xrange = (0.0, 100.0)) 
    return find_zero(k -> get_mpk(t; k = k, n = n0)  - t.δ, xrange)
end 

# No Markup 

get_μ(::Technology) = 1
get_Π(t::Technology, factors) = 0
get_mpk(t::Technology; k, n) = _df_of_k(t.f, k / (t.A * n))
get_y(t::Technology; k, n) = t.f(k / (t.A * n)) * t.A * n

mpk_from_after_tax_rK(::Technology, rK) =  rK
rL_from_mpk(t::Technology, mpk) = mpn_from_mpk(t, mpk) 
mpn_from_mpk(t::Technology, mpk) = t.A * _mpn_from_mpk(t.f, mpk)
k_from_mpk(t::Technology; mpk, n) = t.A * n * _k_from_mpk(t.f, mpk) 

# Markups 

get_μ(t::MarkupTechnology) = t.μ
get_mpk(t::MarkupTechnology; k, n) = (1 - t.m) * _df_of_k(t.f, k / (t.A * n))
get_y(t::MarkupTechnology; k, n) = (1 - t.m) * t.f(k / (t.A * n)) * t.A * n - t.X
get_Π(t::MarkupTechnology, factors) = (1 - 1 / t.μ) * get_y(t, factors) - t.X / t.μ

mpk_from_after_tax_rK(t::MarkupTechnology, rK) = t.μ * rK
rL_from_mpk(t::MarkupTechnology, mpk) = mpn_from_mpk(t, mpk) / t.μ

function mpn_from_mpk(t::MarkupTechnology, mpk) 
    @unpack A, m, f = t
    A * (1 - m) * _mpn_from_mpk(f, mpk / (1 - m))
end 

function k_from_mpk(t::MarkupTechnology; mpk, n)
    @unpack A, m, f = t
    A * n * _k_from_mpk(f, mpk / (1 - m)) 
end 



