# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Julia 8 Threads 1.6.1
#     language: julia
#     name: julia-8-threads-1.6
# ---

# # Case with no capital tax and no markups

# Here we do the simulation of the case where there is no capital tax, and as a result, the level of capital adjusts with the interest rate. This is done for the economy without markups.

# + tags=[]
import Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using Revise
using Aiyagari
using StructArrays
using Plots
using ProgressMeter
using DelimitedFiles
using FastClosures

# + tags=[]
ProgressMeter.ijulia_behavior(:clear);
pgfplotsx();
default(label = "", lw = 2, dpi = 300, left_margin = 0Plots.mm, format=:svg);
# -

_LOAD_GUESSES = true
_SAVE_GUESSES = false
_ITERS = 1

# ## Benchmark calibration

# Dirk and Kurt's calibration (almost)
e = let
    P, z_vals = let
        ar1 = 0.9695
        sigmaP = sqrt(0.0384)/(1.2)
        sigmaIID = sqrt(0.0522)/(1.2)
         calibration(5, 2 , ar1, sigmaP, sigmaIID)
    end
    # Technology
    t = let
        θ = 0.3
        ls = 1 - θ
        δ = 0.1
        #α1, A1 = get_tech_params(1, θ = θ)
        A1 = 1.0
        α1 = 0.3
        Technology(f = CobbDouglas(α = α1), δ = δ)
    end
    # Households
    h = let
        ies = 1.0
        β = 0.99 #* (1 + g)^(1 - 1/ies)
        Household(
            u = EZ(ies = 1.0, ra = 5.5),
            v = GHH(θ = 1.0, ν = 0.2),
            P = P, z_grid = z_vals, β = β, a_max = 100.0)
    end
    Economy(h = h, t = t)
end

r_range = (-0.0172, -0.0171) # narrowing the range
@time laissez_faire = solve_laissez_faire(e;
    r_range = r_range,
    tol =  (value_function = 1e-10, distribution = 1e-13)
)

# ## New stationary outcome

b_target = laissez_faire.y * 0.60

r_range_2 = (-0.0152, -0.0148)  # narrowing the range
@time final_eq = @closure solve_new_stationary_equilibrium_given_k_b(
        laissez_faire; r_range = r_range_2, tol = (value_function = 1e-7, distribution = 1e-8)
    ) do (r)
        # returns (k, b) consistent with r and no capital taxes
        t = get_t(laissez_faire)
        b = b_target
        rK = rK_from_r(;t, r)
        mpk = mpk_from_after_tax_rK(t, rK)
        k = k_from_mpk(t; mpk, laissez_faire.n)
        return (k, b)
    end

# ## Transition

T = 100
H = 50;

# Smooth debt policy
ρB = 0.9
b_list = Array{Float64,1}(undef, T + H)
b_list[1] = 0.0
b_list[2] = laissez_faire.y * 0.05
b_list[T:end] .= b_target
for i in 3:T-1
    b_list[i] = b_list[2] * ρB^(i-2) + (1 - ρB^(i-2)) * b_target
end

r_path = nothing
if _LOAD_GUESSES
    r_path = try
        # load the transfer vector from previous iterations
        readdlm(joinpath(@__DIR__, "..", "output", "tmp_calcs", "tmpINEFF001.txt"))[:,1]
    catch
        nothing
    end
end;

function generate_k_b_no_τk(laissez_faire, final, b_list)
    f = @closure (r, i) -> begin
        if i <= length(b_list)
            b = b_list[i]
        else
            b = final.b
        end
        t = get_t(laissez_faire)
        rK = rK_from_r(;t, r)
        mpk = mpk_from_after_tax_rK(t, rK)
        k = k_from_mpk(t; mpk, laissez_faire.n)
        return (k, b)
    end
    return f
end

my_k_b_fun = generate_k_b_no_τk(laissez_faire, final_eq, b_list)

transition = solve_transition(
    my_k_b_fun,
    laissez_faire,
    final_eq;
    init_r_path = r_path,
    path_length = T + H,
    iterations = _ITERS,
    m = 10,
    beta = -0.01
);

plot(transition.path.r)

plot(transition.path.transfer)

plot(transition.path.k)
plot!(transition.path.b)

_SAVE_GUESSES && open(joinpath(@__DIR__, "..", "output", "tmp_calcs", "tmpINEFF001.txt"), "w") do io
    writedlm(io, transition.path.r)
end

# ## Satistics and Plots

println("INITIAL STEADY STATE")
println("=====================")
summary_statics(laissez_faire)

println("FINAL STEADY STATE")
println("=====================")
summary_statics(final_eq, path = transition.path, laissez_faire = laissez_faire)


f1 = plot(do_plots(transition, laissez_faire)..., size = (800, 400))

savefig(f1, joinpath(@__DIR__, "..", "output", "figures", "transition_inefficient.pdf"))


