# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: Julia 20 threads 1.7.3
#     language: julia
#     name: julia-20-threads-1.7
# ---

# # The present value of asset supply elasticities

# +
import Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using Revise
using MicroRisks
using ProgressMeter
using LaTeXStrings
using Roots
using CairoMakie
using Polyester

ProgressMeter.ijulia_behavior(:clear);
# -

# Household
h = let 
    ar1 = 0.9695
    sigmaP = sqrt(0.0384)/(1.2)
    sigmaIID = sqrt(0.0522)/(1.2)
    P, z_vals = calibration(5, 2 , ar1, sigmaP, sigmaIID)

    ies = 1
    crra = 5.5
    β = 0.993
    u = EZ(ies = ies, ra = crra, β = β)
    v = GHH(θ = 1.0, ν = 0.2)
    Household(u = u, a_grid = grid(; stop = 10.0, length = 500),
        v = v, P = P, z_grid = z_vals)
end

# Technology
t = let
    δ = 0.1
    A = 0.2
    α = 0.3
    μ = 1.4
    CobbDouglasTechnology(α = α, A = A^((1 - α)), δ = δ, μ = μ)
end

@time e_init = stationary_laissez_faire(h, t; r_range = (-0.02, 0.0), verbose = false)

@time jac = jacobian_column(1:50:500, e_init; cap_t = 500, cap_s = 500, ΔR = 1e-4, ΔT = 0.0);

let
    fig = Figure() 
    ax = Axis(fig[1, 1], xlabel = "s", title = "Jacobian R")
    for col in jac
        lines!(ax, col) 
    end
    fig
end

elas = jac ./ e_init.a .* (1 + e_init.r);

let
    fig = Figure() 
    ax = Axis(fig[1, 1], xlabel = "s", title = "Elasticities R")
    for col in elas
        lines!(ax, col) 
    end
    fig
end

full_jac = jacobian(e_init; cap_t = 1_000, ΔR = 1e-4, ΔT = 0.0);
full_elas = full_jac ./ e_init.a .* (e_init.r + 1);

f, ax, = lines(full_elas' * ones(size(full_elas, 1)), title = "Elasticities to a permanent increase in R")
ax.yticks = 0:5:90 
f

@time jac_T = jacobian_column(1:50:500, e_init; cap_t = 500, cap_s = 500, ΔR = 0.0, ΔT = 1e-4);

let
    fig = Figure() 
    ax = Axis(fig[1, 1], xlabel = "s", title = "Jacobian T")
    for col in jac_T
        lines!(ax, col) 
    end
    fig
end

@time pvs = pv_elasticities(1:20:500, e_init; cap_s = 500, cap_t = 2000);


lines(1:20:500, pvs, title = "PV Formula for different s", xlabel = "s")



