########################################################
# # Utility types
########################################################

abstract type AbstractUtility end


# Helper structs for utility calculations

struct Log end

struct Power{T1, T2}
    m::T1 # exponent
    inv::T2 # inverse exponent
end

Power(x::Real) = Power(convert(AbstractFloat, x), convert(AbstractFloat, 1/x))


# ## CRRA
struct CRRA{T1, T2} <: AbstractUtility
    ra::T1
    par::T2
end

function CRRA(ra_par)
    ra = ra_par == 1 ? Log() : Power(1 - ra_par)
    return CRRA(ra, ra_par)
end

CRRA(; ra=2.0) = CRRA(ra)

Base.show(io::IO, m::CRRA) = print(io, "CRRA(ra = $(m.par))")

get_ies(u::CRRA) = u.ra
get_ra(u::CRRA) = u.ra


# ## Epstein-Zin

# Epstein-Zin type
struct EZ{T1, T2, T3} <: AbstractUtility
    ies::T1 # IES
    ra::T2 # CRRA
    pars::T3 # Original params
end

Base.show(io::IO, m::EZ) = print(io, "EZ(ies = $(m.pars[1]), ra = $(m.pars[2]))")

# Epstein-Zin constructors
function EZ(ies_par, ra_par)
    ies = ies_par == 1 ? Log() : Power(1 - 1/ies_par)
    ra = ra_par == 1 ? Log() : Power(1 - ra_par)
    return EZ(ies, ra, (ies_par, ra_par))
end

EZ(;ies=1/2, ra=2.0) = EZ(ies, ra)

get_ies(u::EZ) = u.ies
get_ra(u::EZ) = u.ra


# ## Disutility of labor

# Fixed Labor
Base.@kwdef struct FixedLabor{R}
    n::R = 1.0
end

Base.show(io::IO, m::FixedLabor) = print(io, "FixedLabor(n = $(m.n))")


# GHH Preferences
Base.@kwdef struct GHH{R}
    θ::R = 1.0 # weight on labor in utility
    ν::R = 0.2 # Frisch elasticity of labor
end

Base.show(io::IO, m::GHH) = print(io, "GHH(θ = $(m.θ), ν = $(m.ν))")


####################################################################
# # Household type
###################################################################

Base.@kwdef struct Household{R, U<:AbstractUtility, V1, V2, V3, V4, M}
    β::R = 0.95
    u::U = CRRA(2)
    v::V4 = GHH()
    z_grid::V1 = @SArray [0.5, 1.0]
    n::Int64 = length(z_grid)
    P::M = @SMatrix [0.5 0.5; 0.2 0.8];
    grid_points::Int64 = 10_000
    a_max::R = 15.0
    a_min::R = 0.0
    a_grid::V2 = collect(LinRange(a_min, a_max, grid_points))
    Pss::V3 = ergodic(P)
end


function Base.show(io::IO, m::Household)
    @unpack β, u, v,  z_grid, P, grid_points, a_max, a_min = m
    if length(z_grid) > 2
        z_grid_string = "[$(z_grid[1])..$(z_grid[end])]"
        P_string = "[..]"
    else
        z_grid_string = "$z_grid"
        P_string = "$P"
    end
    print(io, "$u, β=$β, v=$v, z_grid=$z_grid_string, P=$P_string, pts=$grid_points, a_max=$a_max, a_min=$a_min")
end

get_u(h::Household) = h.u
get_v(h::Household) = h.v

#####################################################
# # Taxes type
#####################################################

abstract type AbstractTaxes end

struct LinearIncomeTaxes{R} <: AbstractTaxes
    τn::R  # labor tax
    τk::R  # capital tax
    τπ::R  # profit tax
end

LinearIncomeTaxes(;τn, τk, τπ) = LinearIncomeTaxes(τn, τk, τπ)

Base.show(io::IO, τ::LinearIncomeTaxes) = print(io, "LinearIncomeTaxes(τn=$(τ.τn), τk=$(τ.τk), τπ=$(τ.τπ))")


#######################################################
# # Technology types
#######################################################

abstract type AbstractProductionFunction end

Base.@kwdef struct CobbDouglas{R} <: AbstractProductionFunction
    α::R = 0.33 # capital share
end

Base.show(io::IO, f::CobbDouglas) = print(io, "CobbDouglas(α=$(f.α))")

Base.@kwdef struct CES{R} <: AbstractProductionFunction
    α::R = 0.33 # capital share
    ρ::R = 0.5 # elasticity between capital and labor
end

Base.show(io::IO, f::CES) = print(io, "CES(α=$(f.α), ρ=$(f.ρ))")


######################################################

abstract type AbstractTechnology end

get_δ(t::AbstractTechnology) = t.δ


# Technology struct incorporating TFP, and depreciation
Base.@kwdef struct Technology{F <: AbstractProductionFunction, R} <: AbstractTechnology
    f::F = CobbDouglas()
    A::R = 1.0
    δ::R = 0.1
end

function Base.show(io::IO, t::Technology)
    @unpack f, A, δ = t
    print(io, "$f, A=$A, δ=$δ")
end

get_μ(::Technology) = 1


## Technology with markup struct, TFP, and depreciation
#
# This follows a modified description of Ball and Mankiw (2020):
#       Market Power In Neoclassical Growth Models
#       NBER WP 28538
Base.@kwdef struct MarkupTechnology{F <: AbstractProductionFunction, R} <: AbstractTechnology
    f::F = CobbDouglas()
    A::R = 1.0
    δ::R = 0.1
    μ::R = 1.0   # Value added markup = (1 - m)γ /(1 - γm) wher γ is firm's markup
    X::R = 0.0   # Fixed cost -- this is in terms of goods (rather than labor as in BM)
    m::R = 0.0   # This represents the paramer α in Ball and Mankiw (2020)
end

function Base.show(io::IO, t::MarkupTechnology)
    @unpack f, A, δ, μ, m, X = t
    print(io, "$f, A=$A, δ=$δ, μ=$μ, m=$m, X=$X")
end

get_μ(t::MarkupTechnology) = t.μ


#######################################################
# # Economy type
#######################################################

abstract type AbstractEconomy end

Base.@kwdef struct Economy{S, T} <: AbstractEconomy
    h::S = Household()
    t::T = Technology()
end

Base.show(io::IO, e::Economy) = print(io, "Economy: $(e.h), $(e.t)")

get_h(e::Economy) = e.h
get_t(e::Economy) = e.t


#########################################################
# # Stationary Economy type
#########################################################

Base.@kwdef struct StationaryEquilibrium{S1, S2, S3, S4, R}
    e::S1 # economy struct
    # fiscal
    b::R
    transfer::R
    r::R
    k::R
    n::R
    w::R
    y::R
    s::R
    # solution details
    v::S2
    pol::S3
    pdf::S4
end

get_e(s::StationaryEquilibrium) = s.e
get_h(s::StationaryEquilibrium) = get_h(s.e)
get_t(s::StationaryEquilibrium) = get_t(s.e)

Base.show(io::IO, s::StationaryEquilibrium) = print(io, "r=$(s.r), b=$(s.b), k=$(s.k), transfer=$(s.transfer),  y=$(s.y), $(s.e)")
