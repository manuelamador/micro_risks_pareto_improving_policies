# -*- coding: utf-8 -*-
# # Comparison with Auclert et al Jacobian Computation 

import Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using MicroRisks
using DelimitedFiles
using CairoMakie

# +
# Household
h = let     
    
    ar1 = 0.966
    sigmaP = 0.5 * (1 - 0.966^2)^(1/2)  # sqrt(0.0384)/(1.2)
    sigmaIID = 0.0 #.5 # sqrt(0.0522)/(1.2)
    P, z_vals = calibration(5, 1, ar1, sigmaP, sigmaIID)

    # This the income process implied by Auclert et al's code.
    P_Auclert = [
        0.933714431521  0.06459062191600005  0.001675545126000003  1.931791600000005e-05  8.352100000000029e-08;
        0.016147655479000014  0.934552204084  0.048457454874000036  0.0008378560840000014  4.829479000000013e-06;
        0.0002792575210000005  0.03230496991600003  0.934831545126  0.03230496991600003  0.0002792575210000005;
        4.829479000000013e-06  0.0008378560840000015  0.04845745487400005  0.934552204084  0.016147655479000014;
        8.352100000000029e-08  1.931791600000005e-05  0.001675545126000003  0.06459062191600005  0.933714431521
    ]
    
    z_vals_Auclert = [ 
        0.32506854402390545,
        0.5359474228930929,
        0.883627915977697,
        1.4568561402539593,
        2.401949706452274,
    ]

    @assert P ≈ P_Auclert 
    @assert z_vals ≈ z_vals_Auclert

    ies = 1 / 2.0
    crra = 2.0
    β = 0.987475
    u = EZ(ies = ies, ra = crra, β = β)
    v = MicroRisks.FixedLabor()
    Household(u = u, a_grid = grid(stop = 30.0, length = 500, scale = :log),
        v = v, P = P, z_grid = z_vals)
end

# Using the same a_grid function as in the Auclert et al version: 
a_grid = readdlm(joinpath(@__DIR__, "assets_grid.csv"), ',')
h.a_grid .= a_grid


# Technology
t = let
    δ = 0.12
    A = 0.748057
    α = 0.3
    μ = 1.0
    CobbDouglasTechnology(α = α, A = A, δ = δ, μ = μ)
end


@time e_init = stationary_laissez_faire(h, t; r_range = (-0.09, 0.02), 
    verbose = true, hh_problem_kwargs = (; verbose= false, max_iters = 20_000))
# -

@assert is_pol_valid(e_init)

# Reading the jacobian from Auclert et al simulations:
mat = readdlm(joinpath(@__DIR__, "jacobian_SSJ_toolkit.csv"), ',');

# Computing ours
@time jac = jacobian_column([10, 200, 500], e_init, cap_t = 1001, cap_s = 500, ΔR = 1e-4, ΔT = 0.0);


# Comparing in a figure
let 
    fig, ax = lines( jac[1][2:end], linewidth = 5)
    lines!(ax, mat[:, 10], linewidth = 4, linestyle = :dash, color = :black)

    lines!(ax, jac[2][2:end], linewidth = 5)
    lines!(ax, mat[:, 200], linewidth = 4, linestyle = :dash, color=:black)

    lines!(ax, jac[3][2:end], linewidth = 5)
    lines!(ax, mat[:, 500], linewidth = 4, linestyle = :dash, color=:black)
    fig
end 

extrema(mat[1:end, 10] .- jac[1][2:end])
extrema(mat[1:end, 200] .- jac[2][2:end])
extrema(mat[1:end, 500] .- jac[3][2:end])


