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

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Aiyagari
using Plots
using PrettyTables
using Roots
using UnPack
using ProgressMeter

ProgressMeter.ijulia_behavior(:clear);
# pgfplotsx();
gr()
default(label = "", lw = 2, dpi = 300, left_margin = 0Plots.mm, format=:png);

# helper function
get_k(r, t, n) = k_from_mpk(t; mpk = mpk_from_after_tax_rK(t, r + t.δ), n)

# # Mark-up Economy

e_μ = let
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
r_range = (-0.0171, -0.016)
@time laissez_faire_μ = solve_laissez_faire(e_μ;
    r_range = r_range,
    tol =  (value_function = 1e-10, distribution = 1e-13)
)

# # Dynamic Efficient -- No change in K

# ## new steady state interest rate and eq

final_eq_μ = let
    r = -0.013930022583173161  # using a good guess
    b_target = laissez_faire_μ.y * 0.60
    k_target = laissez_faire_μ.k
    @time solve_new_stationary_equilibrium_given_k_b(
        laissez_faire_μ,
        k_target,
        b_target;
        r_range = (r - 1e-8, r + 1e-8),
        tol = (value_function = 1e-10, distribution = 1e-13)
    )
end

# # trace out household savings

eq_μ = let
    y_0=laissez_faire_μ.y
    b_range = range(-0.5 * y_0, 3 * y_0, length = 20)
    k_target = laissez_faire_μ.k
    eqs = []
    for b in b_range
        push!(eqs,
            solve_new_stationary_equilibrium_given_k_b(
                laissez_faire_μ,
                k_target,
                b;
                r_range = (-0.05, 0.0),
                tol = (value_function = 1e-6, distribution = 1e-6)
            )
        )
    end
    eqs
end;

let
    t = get_t(laissez_faire_μ)
    n0 = laissez_faire_μ.n
    y0 = laissez_faire_μ.y
    y1 = final_eq_μ.y
    ky_0 = laissez_faire_μ.k / y0
    ky_1 = final_eq_μ.k / y1
    sy_1 = final_eq_μ.s / y1
    r_1 = final_eq_μ.r
    r_0 = laissez_faire_μ.r

    plot([(eq.s / y0, eq.r) for eq in eq_μ ], color = :black)
    plot!([(get_k(r, t, n0)/y0, r) for r in range(-0.03, 0.01, length = 15)], color = :red)
    plot!([(k / y0, mpk_from_factors(t, k = k, n = n0) - t.δ) for k in range(2 * y0, 5 * y0, length = 15)],
        color = :blue, style = :dash)

    ylims!(-0.03, 0.01)
    xlims!(0.0, 5)
    hline!([laissez_faire_μ.r, final_eq_μ.r], color = :black, lw = 0.75)

    plot!([0,ky_0],             # xlims for shade
         [0,0],                   # dummy y coordinates
         fillrange = (r_0,r_1), # apply ylims
         fillalpha = 0.5,         # set transparency
         fillcolor=:red,         # set shade color
         label = "Cost: Δr*K/Y")

    plot!([ky_0,final_eq_μ.s/y0],             # xlims for shade
             [0,0],                   # dummy y coordinates
             fillrange = (r_1,0.0), # apply ylims
             fillalpha = 0.75,         # set transparency
             fillcolor=:gray,         # set shade color
             label = "Revenue: r*dS", legend=false, grid=false)

    hline!([0], color = :black, lw = 2)
    vline!([laissez_faire_μ.k / laissez_faire_μ.y], color = :black, lw = 0.75)
end

savefig(joinpath(@__DIR__, "..", "output", "figures", "CostBenefit_Efficient.pdf"))


# # Crowding in

kGR = golden_rule_k(e_μ.t, laissez_faire_μ.n)

eq_crowdin = let
    y_0 = laissez_faire_μ.y
    eq_crowdin = []
    b_range = range(-0.5 * y_0,0.0, length=5)
    for b in b_range
        push!(eq_crowdin,
            solve_new_stationary_equilibrium_given_k_b(
                laissez_faire_μ,
                laissez_faire_μ.k,
                b;
                r_range = (-0.05, 0.0),
                tol = (value_function = 1e-6, distribution = 1e-6)
            )
        )
    end

    k_range = range(laissez_faire_μ.k,kGR, length = 10)
    for k in k_range
        push!(eq_crowdin,
            solve_new_stationary_equilibrium_given_k_b(
                laissez_faire_μ ,
                k,
                0.0;
                r_range = (-0.05, 0.0),
                tol = (value_function = 1e-6, distribution = 1e-6)
            )
        )
    end

    b_range = range(0.0, 2.5 * y_0, length = 15)
    for b in b_range
        push!(eq_crowdin,
            solve_new_stationary_equilibrium_given_k_b(
                laissez_faire_μ ,
                kGR,
                b;
                r_range = (-0.05, 0.0),
                tol = (value_function = 1e-6, distribution = 1e-6)
            )
        )
    end


    eq_crowdin
end


@time crowdin_final = solve_new_stationary_equilibrium_given_k_b(
    laissez_faire_μ,
    kGR,
    0.6 * laissez_faire_μ.y;
    r_range = (laissez_faire_μ.r, -0.010),
    tol = (value_function = 1e-10, distribution = 1e-13)
)

let
    r_0 = laissez_faire_μ.r
    t = get_t(laissez_faire_μ)
    n0 = laissez_faire_μ.n
    y0 = laissez_faire_μ.y
    y1 = crowdin_final.y
    ky_0 = laissez_faire_μ.k / y0
    ky_1 = crowdin_final.k / y1
    sy_1 = crowdin_final.s / y1
    r_1 = crowdin_final.r

    plot([(eq.s / y0, eq.r) for eq in eq_crowdin], color = :black)
    plot!([(get_k(r, t, n0)/y0, r) for r in range(-0.03, 0.01, length = 15)], color = :red)
    plot!([ (k / y0, mpk_from_factors(t, k = k, n = n0) - t.δ) for k in range(2 * y0, 5 * y0, length = 15)],
        color = :blue, style = :dash)
    ylims!(-0.03, 0.0175)
    xlims!(0.0, 5)
    hline!([laissez_faire_μ.r, crowdin_final.r], color = :black, lw = 0.75)
    vline!([laissez_faire_μ.k / laissez_faire_μ.y, kGR/laissez_faire_μ.y], color = :black, lw = 0.75)

    plot!([0,ky_0],             # xlims for shade
         [0,0],                   # dummy y coordinates
         fillrange = (r_0,r_1), # apply ylims
         fillalpha = 0.5,         # set transparency
         fillcolor=:red,         # set shade color
         label = "Cost: Δr*K/Y")

    plot!([  (k / y0, mpk_from_factors(t, k = k, n = n0) - t.δ) for k in range(laissez_faire_μ.k,kGR, length = 15)],
        fill = (0, 0.5, :blue),lw=0.0)

    plot!([ky_0,crowdin_final.k /y0],             # xlims for shade
        [0,0],                   # dummy y coordinates
        fillrange = (r_1,0.0), # apply ylims
        fillalpha = 0.5,         # set transparency
        fillcolor=:blue,         # set shade color
        label = "Revenue I: (MPK-δ)* ΔK/Y")

    plot!([crowdin_final.k /y0,crowdin_final.s/y0],             # xlims for shade
        [0,0],                   # dummy y coordinates
        fillrange = (r_1,0.0), # apply ylims
        fillalpha = 0.75,         # set transparency
        fillcolor=:gray,         # set shade color
        label = "Revenue II: r*ΔS/Y")

    hline!([0], color = :black, lw = 2,legend=false, grid=false)
end

savefig(joinpath(@__DIR__, "..", "output", "figures", "CostBenefit_Crowdin.pdf"))

# # Competititive Case

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
        A1=1.0
        α1=0.3
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


# Solve laissez-faire economy
r_range = (-0.0172, -0.0171) # narrowing the range
@time laissez_faire = solve_laissez_faire(e;
    r_range = r_range,
    tol =  (value_function = 1e-10, distribution = 1e-13)
)


@time final_eq = let
    tol = (value_function = 1e-7, distribution = 1e-8)
    b_target = laissez_faire.y*0.6
    r_range_2 = (-0.0152, -0.0148)  # narrowing the range
    solve_new_stationary_equilibrium_given_k_b(laissez_faire; r_range = r_range_2, tol) do (r)
        # returns k consistent with r and no capital taxes
        t = get_t(laissez_faire)
        b = b_target
        rK = rK_from_r(;t, r)
        mpk = mpk_from_after_tax_rK(t, rK)
        k = k_from_mpk(t; mpk, laissez_faire.n)
        return (k, b)
    end
end


eq_nomarkup = let  r_range = [-.03,0.0], tol = (value_function = 1e-7, distribution = 1e-8)
    b_range = range(-1.0 * laissez_faire.y, 4.5 * laissez_faire.y, length = 15)
    eq_nomarkup = []
    for b in b_range
        push!(eq_nomarkup,
            solve_new_stationary_equilibrium_given_k_b(laissez_faire; r_range, tol) do  (r)
                # returns k consistent with r and no capital taxes
                t = get_t(laissez_faire)
                rK = rK_from_r(;t, r)
                mpk = mpk_from_after_tax_rK(t, rK)
                k = k_from_mpk(t; mpk, laissez_faire.n)
                return (k, b)
            end
        )
    end
    eq_nomarkup
end


# +
let
    t = get_t(laissez_faire)
    n0 = laissez_faire.n
    y0 = laissez_faire.y
    y1 = final_eq.y
    ky_0 = laissez_faire.k / y0
    ky_1 = final_eq.k / y1
    sy_1 = final_eq.s / y1
    r_1 = final_eq.r

    plot([(eq.s / y0, eq.r) for eq in eq_nomarkup], color = :black)
    plot!([ (k / y0, mpk_from_factors(t, k = k, n = n0) - t.δ) for k in range(2 * y0, 5 * y0, length = 15)], color = :red)

    ylims!(-0.03, 0.01)
    xlims!(0.0, 5)
    hline!([laissez_faire.r, final_eq.r], color = :black, lw = 0.75)

    plot!([0,ky_0],             # xlims for shade
         [0,0],                   # dummy y coordinates
         fillrange = (laissez_faire.r,r_1), # apply ylims
         fillalpha = 0.5,         # set transparency
         fillcolor=:red,         # set shade color
         label = "Cost: Δr*K/Y")

    plot!([ky_0,final_eq.k/y0],             # xlims for shade
         [0,0],                   # dummy y coordinates
         fillrange = (r_1,0.0), # apply ylims
         fillalpha = 0.25,         # set transparency
         fillcolor=:gray,         # set shade color
         label = "Revenue I: r*dK")

    plot!([ky_0,final_eq.s/y0],             # xlims for shade
         [0,0],                   # dummy y coordinates
         fillrange = (r_1,0.0), # apply ylims
         fillalpha = 0.75,         # set transparency
         fillcolor=:gray,         # set shade color
         label = "Revenue II: r*dS")

    hline!([0], color = :black, lw = 2)
    vline!([laissez_faire.k / laissez_faire.y], color = :black, lw = 0.75, legend=false,grid=false)

end
# -

savefig(joinpath(@__DIR__, "..", "output", "figures", "CostBenefit_Nomarkup.pdf"))


