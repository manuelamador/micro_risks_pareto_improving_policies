#######################################################
# Preferences
#######################################################


abstract type AbstractUtility end
abstract type AbstractSubUtility end

struct Log <: AbstractSubUtility
end 
Base.show(io::IO, ::Log) = print(io, "Log")
get_parameter(::Log) = 1 

struct Power{T} <: AbstractSubUtility   
    a::T
end 
Base.show(io::IO, m::Power) = print(io, "$m.a")
get_parameter(m::Power) = m.a


(::Log)(x)  = log(x) 
inverse(::Log, x) = exp(x) 

(m::Power)(x)  = x^(1 - m.a)
inverse(m::Power, x) = x^(1 / (1 - m.a)) 


# Epstein-Zin type
struct EZ{T1<:AbstractSubUtility, T2<:AbstractSubUtility, T3} <: AbstractUtility
    risk::T1 # RA
    temporal::T2 # IES
    β::T3
end
get_ra(m::EZ) = get_parameter(m.risk)
get_inverse_ies(m::EZ) = get_parameter(m.temporal)  
get_ies(m::EZ) = 1/get_parameter(m.temporal)
get_β(m::EZ) = m.β

Base.show(io::IO, m::EZ) = print(io, "EZ(ra=$(get_ra(m)), ies=$(get_ies(m)), β=$(m.β))")

EZ(;ies=1/2.0, ra=2.0, β) = EZ(ra == 1 ? (@warn("Log aggregator for risk not yet implemented. Code may break"); Log()) : Power(ra), ies == 1 ? Log() : Power(1/ies), β)


# CRRA 

struct CRRA{T1<:AbstractSubUtility, T3} <: AbstractUtility
    risk::T1 # RA
    β::T3
end

get_ra(m::CRRA) = get_power(m.risk)
get_inverse_ies(m::CRRA) = get_power(m.risk)  
get_ies(m::CRRA) = 1/get_power(m.risk)
get_β(m::CRRA) = m.β

Base.show(io::IO, m::CRRA) = print(io, "CRRA(ra=$(get_ra(m)), β=$(m.β))")

CRRA(; ra=2.0, β) = CRRA(ra == 1 ? Log() : Power(ra), β)


# ## Disutility of labor

# Fixed Labor
Base.@kwdef struct FixedLabor{R}
    n::R = 1.0
end

Base.show(io::IO, m::FixedLabor) = print(io, "FixedLabor(n=$(m.n))")


# GHH Preferences
Base.@kwdef struct GHH{R}
    θ::R = 1.0 # weight on labor in utility
    ν::R = 0.2 # Frisch elasticity of labor
end

Base.show(io::IO, m::GHH) = print(io, "GHH(θ=$(m.θ), ν=$(m.ν))")



#######################################################
# Household
#######################################################

Base.@kwdef struct Household{U<:AbstractUtility, V, M1, M2, M3, M4}
    u::U = CRRA(; ra = 2, β = 0.95)
    v::V = GHH()
    z_grid::M1 = [0.5, 1.0]
    P::M2 = [0.5 0.5; 0.2 0.8]
    Pprime::M2 = copy(P')
    a_grid::M3 = grid(; start = 0.0, stop = 15.0, length = 500, scale = :log)  
    Pss::M4 = ergodic(P)
end


function Base.show(io::IO, m::Household)
    (; u, v,  z_grid, P, a_grid) = m
    if length(z_grid) > 2
        z_grid_string = "[$(z_grid[1])..$(z_grid[end])]"
        P_string = "[..]"
    else
        z_grid_string = "$z_grid"
        P_string = "$P"
    end
    print(io, "$u, v=$v, z_grid=$z_grid_string, P=$P_string, a_grid = $(first(a_grid))..$(length(a_grid))pts..$(last(a_grid))")
end


######## HouseholdWorkspace

struct HouseholdWorkspace{H, T1, T2, S}
    h :: H
    v :: T1   
    η :: T1   
    a_pol :: T1  
    pdf :: T1
    lower_index :: S
    lower_weight :: T1

    v_tmp :: T2 
    η_tmp :: T2  
    a_tmp :: T2  
    pdf_tmp :: T2 
end 


function HouseholdWorkspace(h::Household)  
    (; v, η, a_pol, pdf, lower_index, lower_weight, a_tmp) = _generate_base_workspace_matrices(h)

    v_tmp = similar(v)
    η_tmp = similar(v)
    pdf_tmp = similar(v)
     
    return HouseholdWorkspace(h, v, η, a_pol, pdf, lower_index, lower_weight, v_tmp, η_tmp, a_tmp, pdf_tmp)
end 
HouseholdWorkspace(;h, R, T, w) = _initialize_HH_ws!(HouseholdWorkspace(h), R, T, w)


function _generate_base_workspace_matrices(h::Household)
    η = Array{eltype(h.a_grid)}(undef, length(h.a_grid), length(h.z_grid))
    v = similar(η) 
    a_pol = similar(η)
    pdf = similar(η)
    lower_index = similar(η, Int)
    lower_weight = similar(η)
    a_tmp = similar(η)  

    return (; η, v, a_pol, pdf, lower_index, lower_weight, a_tmp)
end 


function _initialize_HH_ws!(ws::HouseholdWorkspace, R, T, w)
    h = ws.h
    u = h.u
    ξ = get_inverse_ies(u) 
    a_min = first(h.a_grid)
    for i in axes(ws.v, 1)
        for s in axes(ws.v, 2)
            z = h.z_grid[s]
            dis = disutility_given_w(h.v; w = w * z)
            lab = labor_income(h.v; w = w * z)
            c = R * h.a_grid[i] + T + lab - a_min
            x = c - dis
            ws.η[i, s] = R^(-1/ξ) * x
            ws.v[i, s] = x
            ws.pdf[i, s] = h.Pss[s] # initializing at erodic weights of z
        end
    end
    rdiv!(ws.pdf, sum(ws.pdf)) # normalizing to distribution

    return ws 
end 


function Base.show(io::IO, ws::HouseholdWorkspace)
    print(io, "Worskpace for $(ws.h)")
end


#######################################################
# Technology
#######################################################


# Technology struct incorporating TFP, and depreciation
Base.@kwdef struct CobbDouglasTechnology{R}
    α::R = 0.33 
    A::R = 1.0
    δ::R = 0.1
    μ::R = 1.0  # = 1.0 means no markip
end

function Base.show(io::IO, t::CobbDouglasTechnology)
    (; α, A, δ) = t
    print(io, "α=$α, A=$A, δ=$δ")
end


####################################################
# Stationary Equilibrium
####################################################

Base.@kwdef struct StationaryEquilibrium{H, T1, R1, W}
    h :: H
    t :: T1  # technology
    ws :: W  # HouseholdWorkspace
    w :: R1
    r :: R1
    T :: R1  # transfer
    n :: R1
    a :: R1
    k :: R1 
    b :: R1 
end 


function Base.show(io::IO, e::StationaryEquilibrium)
    (; h, t, r, w, T, a, k, b, n) = e
    print(io, "Stationary Equilibrium. Household: $h, Technology: $t, r=$r, w=$w, T=$T, a=$a, k=$k, b=$b, n=$n")
end



####################################################
# Jacobian Struct
####################################################

mutable struct JacobianCache{S, I, R, T}
    up :: S
    down :: S
    cap_s :: I
    cap_t:: I 
    ws :: T
    R :: R 
    T :: R 
    ΔR :: R
    ΔT :: R
end 