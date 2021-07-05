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

# # Dynamic Efficient Case

# + tags=[]
import Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()
using Revise
using Aiyagari
using StructArrays
using ProgressMeter
using DelimitedFiles
using Roots
using Plots
using LaTeXStrings
# -

pgfplotsx();
Threads.nthreads()

ProgressMeter.ijulia_behavior(:clear);
default(label = "", lw = 2, dpi = 300, left_margin = 0Plots.mm, format=:svg);

# +
_LOAD_GUESSES = true # load the initial starting guesses for zeros from disk
_SAVE_GUESSES = false # save zeros to disk
_ITERS = 1

# Setting _LOAD_GUESSES to false will recompute all the calculations from scratch. At that point,
# set _ITERS to a higher number.
# -

# ## The benchmark calibration

#' The parameters for technology and households.
e = let
    #Household:
    ar1 = 0.9695
    sigmaP = sqrt(0.0384)/(1.2)
    sigmaIID = sqrt(0.0522)/(1.2)
    P, z_vals = calibration(5, 2 , ar1, sigmaP, sigmaIID)

    ies = 1.0
    crra = 5.5
    β = 0.993

    hh = Household(u = EZ(ies = ies, ra = crra), grid_points = 5_000,
        v = GHH(θ = 1.0, ν = 0.2), P = P, z_grid = z_vals, β = β, a_max = 10.0)

    #Technology
    δ = 0.1
    μ = 1.4
    α = .3
    #ρ = 0.7
    A_μ = 0.2
    # t_μ = MarkupTechnology(f = CES(α = α, ρ = 0.7), δ = δ, μ = μ, A = A_μ)
    t_μ = MarkupTechnology(f = CobbDouglas(α = α), δ = δ, μ = μ, A = A_μ)

    Economy(h = hh, t = t_μ)
end

# Solve laissez-faire economy
@time laissez_faire = let r_range = (-0.0170, -0.0169)
    solve_laissez_faire(e;
        r_range = r_range,
        tol =  (value_function = 1e-10, distribution = 1e-13)
    )
end

# ## Transitions

T = 100  # period of adjustment
H = 50   # stationary part

# ## Transition with constant k and smooth b

# We target a debt of 59% percent of output

b_target = laissez_faire.y * 0.60
k_target = laissez_faire.k

@time final_eq_1 = solve_new_stationary_equilibrium_given_k_b(
    laissez_faire,
    k_target,
    b_target;
    r_range = (laissez_faire.r, -0.010),
    tol = (value_function = 1e-10, distribution = 1e-13)
)

# Smooth debt policy
ρB = 0.9
b_list = Array{Float64,1}(undef, T + H)
b_list[1] = 0.0
b_list[2] = laissez_faire.y * 0.05
b_list[T:end] .= b_target
for i in 3:T-1
    b_list[i] = b_list[2] * ρB^(i-2) + (1 - ρB^(i-2)) * b_target
end

# k remains constant
k_list = [laissez_faire.k for _ in b_list];

# Computes the transition

r_path = nothing
if _LOAD_GUESSES
    r_path = try
        # load the transfer vector from previous iterations
        readdlm(joinpath(@__DIR__,"..", "output", "tmp_calcs", "tmp_r_001.txt"))
    catch
        nothing
    end
end;

# + tags=[]
@time transition = solve_transition(
    laissez_faire,
    final_eq_1,
    k_list,
    b_list;
    init_r_path = r_path,
    iterations = _ITERS);

# + tags=[]
_SAVE_GUESSES && open(joinpath(@__DIR__, "..", "output", "tmp_calcs", "tmp_r_001.txt"), "w") do io
    writedlm(io, transition.path.r)
end
# -

# ### Plots

plot(transition.path.transfer, fill = 0, alpha = 0.3,)
plot!(transition.path.transfer, color=1)
hline!([final_eq_1.transfer], lw = 1, style = :dash)

plot(transition.path.r, legend=false)
hline!([laissez_faire.r], lw = 1)
hline!([final_eq_1.r], lw = 1)

# ## Transition with increasing k and increasing b

# We go to the golden rule k -- slowly.

k_golden = golden_rule_k(e.t, laissez_faire.n)
ρK = 0.95
k_list_2 = similar(k_list)
for i in eachindex(k_list_2)
    k_list_2[i] = laissez_faire.k * ρK^(i-1) + (1 - ρK^(i-1)) * k_golden
end

# The debt path and its target remains as in the previous case

b_list_2 = b_list
b_target_2 = b_list[end]

# Computing the new stationary equilibrium in this case:

@time final_eq_2 = solve_new_stationary_equilibrium_given_k_b(
    laissez_faire,
    k_list_2[end],
    b_target_2;
    r_range = (laissez_faire.r, -0.010),
    tol = (value_function = 1e-10, distribution = 1e-13)
)

# Computing the transition:

r_path_2 = nothing
if _LOAD_GUESSES
    r_path_2 = try
        # load the transfer vector from previous iterations
        readdlm(joinpath(@__DIR__, "..", "output", "tmp_calcs", "tmp_r_002.txt"))
    catch b
        nothing
    end
end;

# + tags=[]
@time transition_2 = solve_transition(
    laissez_faire,
    final_eq_2,
    k_list_2,
    b_list_2;
    init_r_path = r_path_2,
    beta = -0.02,
    m = 10,
    iterations = _ITERS
);
# -

_SAVE_GUESSES && open(joinpath(@__DIR__, "..", "output", "tmp_calcs", "tmp_r_002.txt"), "w") do io
    writedlm(io, transition_2.path.r)
end

# ### Plots

plot(transition_2.path.transfer, fill = 0, alpha = 0.3)
plot!(transition_2.path.transfer, color = 1)
hline!([0, final_eq_2.transfer], lw = 1)

plot(transition_2.path.r, legend=false)
hline!([laissez_faire.r], lw = 1)
hline!([final_eq_2.r])

# #### Comparison with previous case

f1 = plot(transition_2.path.v[1] ./ laissez_faire.v .- transition.path.v[1] ./ laissez_faire.v,
    legend = :bottom,
    title = "Welfare gain with k increasing wrt to k fixed")
f2 = plot(transition.path.v[1] ./ laissez_faire.v .- 1)
plot(f1, f2, size = (1000, 400))

plot(transition_2.path.transfer, label="new")
plot!(transition.path.transfer, label="old")

plot(transition_2.path.r, label="new")
plot!(transition.path.r, label="old")

plot(
    [ -(1 + r) * b + bprime for (b, bprime, r) in
        zip(transition.path.b, transition.path.b[2:end], transition.path.r)],
    label = "old"
)
plot!(
    [ -(1 + r) * b + bprime for (b, bprime, r) in
        zip(transition_2.path.b, transition_2.path.b[2:end], transition_2.path.r)],
    label = "new"
)

# ## Transition with increasing k and no debt

# No debt:

b_list_3 = [0.0 for _ in k_list]
b_target_3 = b_list_3[end];

# But the same increase in k as before:

k_list_3 = k_list_2;

# Computing the new stationary equilibrium

@time final_eq_3 = solve_new_stationary_equilibrium_given_k_b(
    laissez_faire,
    k_list_3[end],
    b_target_3;
    r_range = (laissez_faire.r, -0.010),
    tol = (value_function = 1e-10, distribution = 1e-13)
)

r_path_3 = nothing
if _LOAD_GUESSES
    r_path_3 = try
        readdlm(joinpath(@__DIR__, "..", "output", "tmp_calcs", "tmp_r_003.txt"))
    catch
        nothing
    end
end;

# + tags=[]
@time transition_3 = solve_transition(
    laissez_faire,
    final_eq_3,
    k_list_3,
    b_list_3;
    init_r_path = r_path_3,
    beta = -0.02,
    m = 10,
    iterations = _ITERS
);
# -

_SAVE_GUESSES && open(joinpath(@__DIR__, "..", "output", "tmp_calcs", "tmp_r_003.txt"), "w") do io
    writedlm(io, transition_3.path.r)
end

# ### Plots

plot(transition_3.path.transfer,
    fill = 0, alpha = 0.3, label = "transfers", legend = :bottom)
plot!(transition_3.path.transfer, color = 1)
hline!([0, final_eq_3.transfer], lw = 1)

# + tags=[]
plot(transition_3.path.r, legend = :bottom, label = "r")
hline!([laissez_faire.r, final_eq_3.r], lw = 1)
# -

# ## Statistics

println("INITIAL STEADY STATE")
println("=====================")
summary_statics(laissez_faire)


println("\nFINAL STEADY STATE -- CONSTANT K AND DEBT")
println("=====================")
summary_statics(final_eq_1,
    laissez_faire = laissez_faire,
    path = transition.path)


println("FINAL STEADY STATE -- GOLDEN K AND DEBT")
println("=====================")
summary_statics(final_eq_2,
    laissez_faire = laissez_faire,
    path = transition_2.path)


# ## Figures for paper


# ### Transition plots

f1 = plot(do_plots(transition, laissez_faire)..., size = (1000, 500))

f2 = plot(do_plots(transition_2, laissez_faire)..., size = (1000, 500))

f3 = plot(do_plots(transition_3, laissez_faire)..., size = (1000, 500))

# ### Steady state segniorage/transfer plot

# We need to increase a_max so that it does not bind as we increase b.

# Increasing the amax so that it doesn't bind
e_2 = let
    #Household:
    ar1 = 0.9695
    sigmaP = sqrt(0.0384)/(1.2)
    sigmaIID = sqrt(0.0522)/(1.2)
    P, z_vals = calibration(5, 2 , ar1, sigmaP, sigmaIID)

    ies = 1.0
    crra = 5.5
    β = 0.993

    hh = Household(u = EZ(ies = ies, ra = crra), grid_points = 5_000,
        v = GHH(θ = 1.0, ν = 0.2), P = P, z_grid = z_vals, β = β, a_max = 50.0)

    #Technology
    δ = 0.1
    μ = 1.4
    α = .3
    #ρ = 0.7
    A_μ = 0.2
    # t_μ = MarkupTechnology(f = CES(α = α, ρ = 0.7), δ = δ, μ = μ, A = A_μ)
    t_μ = MarkupTechnology(f = CobbDouglas(α = α), δ = δ, μ = μ, A = A_μ)

    Economy(h = hh, t = t_μ)
end

# Solve laissez-faire economy
@time laissez_faire_2 = let r_range = (-0.0175, -0.016)
    solve_laissez_faire(e_2;
        r_range = r_range,
        tol =  (value_function = 1e-10, distribution = 1e-13)
    )
end

r_guess = nothing
if _LOAD_GUESSES
    r_guess = try
        # load the transfer vector from previous iterations
        readdlm(joinpath(@__DIR__, "..", "output", "tmp_calcs", "tmp_ss_b_r.txt"))
    catch
        nothing
    end
end;

out = let
    range = [laissez_faire_2.r, 0.0]
    f = (b, i) -> begin
        print(i, " ", b, ", ")
        if !isnothing(r_guess)
            myrange = [r_guess[i, 2] - 1e-5, r_guess[i, 2] + 1e-5]
        else
            myrange = range
        end
        sol = solve_new_stationary_equilibrium_given_k_b(
            laissez_faire_2,
            laissez_faire_2.k,
            b * laissez_faire_2.y;
            r_range = myrange,
            tol = (value_function = 1e-7, distribution = 1e-7),
            verbose = false,
            find_zero_args = (xatol = 1e-5, maxevals = _ITERS)
        )
        range[2] = sol.r  # iterating down -- change top value of range
        sol
    end
    [f(i, b) for (b, i) in enumerate(reverse(0.2:0.2:7.0))]
end
push!(out, laissez_faire_2);

_SAVE_GUESSES && open(joinpath(@__DIR__, "..", "output", "tmp_calcs", "tmp_ss_b_r.txt"), "w") do io
    writedlm(io, [(eq.b, eq.r) for eq in out])
end

f4 = let
    by = [eq.b / eq.y for eq in out]
    rb = [-eq.r * eq.b / eq.y for eq in out]
    deltark = [(eq.r - laissez_faire_2.r) * eq.k / eq.y for eq in out]
    plot(by, rb,  label = L"- r  b / y", legend = :bottom, size = (1000/3, 500/2))
    plot!(by, deltark, fillrange = rb, color = 1, alpha = 0.2)
    plot!(by, deltark, color = 2, label = L" (r - r_0)  k_0 / y")
end

# ## Exporting the figures

savefig(f1, joinpath(@__DIR__, "..", "output", "figures", "transition_efficient_fixed_k.pdf"))

savefig(f2, joinpath(@__DIR__, "..", "output", "figures", "transition_efficient_golden_k.pdf"))

savefig(f3, joinpath(@__DIR__, "..", "output", "figures", "transition_efficient_no_debt.pdf"))

savefig(f4, joinpath(@__DIR__, "..", "output", "figures", "steady_state_transfers.pdf"))


