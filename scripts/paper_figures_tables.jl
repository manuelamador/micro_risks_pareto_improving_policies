# -*- coding: utf-8 -*-
# # Micro Risks and Pareto Improving Policies 

import Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using Revise
using MicroRisks
using ProgressMeter
using CairoMakie
using LaTeXStrings
using StatsBase
using Polyester
using Roots

ProgressMeter.ijulia_behavior(:clear);

# ## Benchmark Economy

# Household
h = let 
    # labor supply
    v = GHH(θ = 1.0, ν = 0.2)
 
    # income process
    ar1 = 0.9695
    sigmaP = sqrt(0.0384)/(1 + v.ν)
    sigmaIID = sqrt(0.0522)/(1 + v.ν)
    P, z_vals = calibration(5, 2 , ar1, sigmaP, sigmaIID)

    # consumption preferences
    ies = 1 
    crra = 5.5
    β = 0.993
    u = EZ(ies = ies, ra = crra, β = β)

    Household(u = u, a_grid = grid(stop = 10.0, length = 500, scale = :log),
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

# Initial equilibrium (laissez faire)

@time e_init = stationary_laissez_faire(h, t; r_range = (-0.02, 0.0), verbose = true)


# ## Constant-K Transition

# Transition to a higher debt level

# b = 0.6 y0 
b_target = y(e_init) * 0.60

# Final equilibrium with higher debt and same k

@time e_final = stationary_equilibrium_given_k_b(e_init, e_init.k, b_target; r_range = (-0.02, 0.0), verbose = true)


# Debt policy and capital (constant) along the transition

b_path = let   # Smooth path of increasing debt  
    T = 100  # period of adjustment of debt
    H = 50   # debt level no longer moving
    ρB = 0.9
    b_list = Array{Float64,1}(undef, T + H)
    b_list[1] = 0.0
    b_list[2] = y(e_init) * 0.05
    b_list[T:end] .= b_target
    for i in 3:T-1
        b_list[i] = b_list[2] * ρB^(i-2) + (1 - ρB^(i-2)) * b_target
    end
    b_list
end;

k_path = [e_init.k for _ in b_path];

# Solving the transition 

@time path = solve_transition(e_init, e_final; k_path, b_path, 
                nlsolve_kwargs = (; ftol = 1e-7));

# Implied aggreate savings elasticities: 

elas = b_path ./ (path.r .- e_init.r) ./ (e_init.a) .* (1 + e_init.r)
( elas[2], elas[end] )

# ## Transition Towards Golden Rule

k_golden = golden_rule_k(t; e_init.n)

@time e_final_2 = let 
    b_target = y(e_init) * 0.60
    stationary_equilibrium_given_k_b(e_init, k_golden, b_target; r_range = (-0.02, 0.0), verbose = true)
end 

# Debt policy is same as above. Capital along the transition is:

k_path_2 = let 
    ρK = 0.95
    k_path_2 = similar(k_path)
    for i in eachindex(k_path_2)
        k_path_2[i] = e_init.k * ρK^(i-1) + (1 - ρK^(i-1)) * k_golden
    end
    k_path_2
    end;

@time path_2 = solve_transition(e_init, e_final_2; k_path = k_path_2, b_path, 
    nlsolve_kwargs = (; beta = -.03, ftol = 1e-7));

# ## Transition To Golden Rule Without Debt 

@time e_final_3 = let 
    b_target = 0.0
    stationary_equilibrium_given_k_b(e_init, k_golden, b_target; r_range = (-0.02, 0.0), verbose = true)
end 

b_path_3 = [0.0 for _ in k_path_2];

@time path_3 = solve_transition(e_init, e_final_3; 
    k_path = k_path_2, b_path = b_path_3, 
    nlsolve_kwargs = (; beta = -.03, ftol = 1e-7));

# ## Statistics

statistics = let
    iobuffer = IOBuffer()

    println(iobuffer, "INITIAL STEADY STATE")
    println(iobuffer,"=====================")
    summary_statics(e_init; iobuffer)


    println(iobuffer, "\nFINAL STEADY STATE -- CONSTANT K AND DEBT")
    println(iobuffer, "=====================")
    summary_statics(e_final;
        laissez_faire = e_init,
        path = path, iobuffer)


    println(iobuffer, "\nFINAL STEADY STATE -- GOLDEN K AND DEBT")
    println(iobuffer, "=====================")
    summary_statics(e_final_2;
        laissez_faire = e_init,
        path = path_2, iobuffer)

    String(take!(iobuffer))
end
println("\n", statistics, "\n")

# ##  Plots

function do_plots(path, laissez_faire; legend_pos = :rt, last_plot = false)
    t = laissez_faire.t
    h = laissez_faire.h
    n = laissez_faire.n

    ylist  = map(k -> output(t; k , n), path.k)
    size = (400, 300)
    y0 = y(laissez_faire)


    #################
    #  B

    fontsize_theme = Theme(fontsize = 15)
    set_theme!(fontsize_theme)

    fig = Figure(resolution = (1000, 500)) 


    ax = Axis(fig[1, 1], title = L"$b/y_{0}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, 100))
    ylims!(ax,  (0, 100))
    ax.xticks = 0:25:100
    
    if last_plot  
        lines!(ax, 100 .* path.b ./ y0, color = :blue, linewidth = 4)
    else
        lines!(ax, 100 .* path.b ./ y0, color = :blue)
    end

    ##################
    # K

    ax = Axis(fig[1, 2], title = L"$k/y_{0}$", titlesize = 20)
    xlims!(ax,  (0, 100))
    ylims!(ax, (0, 1.1 * maximum(path.k ./ y0)))
    ax.xticks = 0:25:100
   
    lines!(ax, path.k ./ y0, color = :blue)

    # ##################
    # # TRANSFERS AND B 

    ax = Axis(fig[1, 3], title = LaTeXString("\$T\$ and  \$b^\\prime - R b\$"), ylabel = LaTeXString("%"),  titlesize = 20)
    xlims!(ax,  (0, 100))
    ax.xticks = 0:25:100

    ser = 100 .* path.T ./ y0
    band!(1:length(ser), ser, [0.0 for _ in ser], color = (:blue, 0.2), label = L"$T/y_0$")
    lines!(ser, color = :blue)
    lines!([100 * (bprime - (1 + r) * b) / y0 for (b, bprime, r) in zip(path.b, path.b[2:end], path.r, ylist)], label = L"$\frac{b^\prime - R b}{y_0}$", color = :red, linestyle = :dash)

    axislegend(position = legend_pos)

    # ######################
    # #  R

    ax = Axis(fig[2, 1], title = L"$r$", ylabel = LaTeXString("%"), titlesize = 20) 
    xlims!(ax,  (0, 100))
    ax.xticks = 0:25:100
   
    lines!(ax, path.r .* 100, color = :blue)

    # ######################
    # # CONSUMPTION
   
    ax = Axis(fig[2, 2], title = L"$c/c_0$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, 100))
    ylims!(ax, (-5, 5))
    ax.xticks = 0:25:100

    agg_c = [y + (1 - t.δ) * k - kprime for (y, k, kprime) in zip(ylist, path.k, path.k[2:end])]
    c0 = y(laissez_faire) - t.δ * laissez_faire.k

    band!(ax, 1:length(agg_c), 100 .* (agg_c ./ c0  .- 1), [0.0 for _ in agg_c], color = (:blue, 0.2))
    lines!(ax, 100 .* (agg_c ./ c0  .- 1), color = :blue)
    
    # #########################
    # # STD DEV

    ax = Axis(fig[2, 3],  title = L"$\sigma(\log c) - \sigma(\log c_0)$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, 100))
    ax.xticks = 0:25:100
   
    std_con = map(path) do alloc
        con = reshape(consumption_alloc(h; R = 1 + alloc.r, laissez_faire.w, alloc.T, alloc.a_pol), 1, :)[1, :]
        allpdf = reshape(alloc.pdf, 1, :)[1, :]
        100 * StatsBase.std(log.(con), StatsBase.ProbabilityWeights(allpdf))
    end

    con_init = consumption_alloc(h; R = 1 + laissez_faire.r, laissez_faire.w, laissez_faire.T, laissez_faire.ws.a_pol)
    con_init_ = reshape(con_init,1, :)[1, :]
    allpdf_init = reshape(laissez_faire.ws.pdf, 1, :)[1, :]
    std_con_init = 100 * StatsBase.std(log.(con_init_), StatsBase.ProbabilityWeights(allpdf_init))

    band!(ax, 1:length(std_con), std_con .- std_con_init, [0 for _ in std_con], color = (:blue, 0.2))
    lines!(ax, std_con .- std_con_init, color = :blue)
    
    return fig
end

f1 = do_plots(path, e_init)

f2 = do_plots(path_2, e_init)

f3 = do_plots(path_3, e_init, legend_pos = :rb, last_plot = true)

# ##  Steady State Segniorage Plots

# Increasing the amax so that it doesn't bind
h_2 = let 
    v = GHH(θ = 1.0, ν = 0.2)

    ar1 = 0.9695
    sigmaP = sqrt(0.0384)/(1 + v.ν)
    sigmaIID = sqrt(0.0522)/(1 + v.ν)
    P, z_vals = calibration(5, 2 , ar1, sigmaP, sigmaIID)

    ies = 1 
    crra = 5.5
    β = 0.993
    u = EZ(ies = ies, ra = crra, β = β)
    Household(u = u, a_grid = grid(stop = 50, length = 500, scale = :log), 
        v = v, P = P, z_grid = z_vals)
end

# Solve laissez-faire economy
@time e_init_2 = stationary_laissez_faire(h_2, t; r_range = (-0.02, 0.0), verbose = true)

r_range_ = [e_init_2.r, 0.0]
b_range  =  2.6:-0.2:0.2
out = Array{Any}(undef, length(b_range))
Polyester.disable_polyester_threads() do
    p = Progress(length(b_range))
    @time Threads.@threads for i in eachindex(b_range)
        b = b_range[i]
        sol = stationary_equilibrium_given_k_b(
            e_init_2,
            e_init_2.k,
            b * y(e_init_2);
            r_range = r_range_,
            verbose = false,
            hh_problem_kwargs = (;
                value_tol = 1e-7, 
                policy_tol = 1e-7, 
                pdf_tol = 1e-7)
        )
        out[i] = sol
        next!(p)
    end
end
# out = @showprogress map(f, 2.6:-0.2:0.2)
push!(out, e_init_2);

f4 = let
    b_y = [eq.b / y(eq) for eq in out]
    rb = [-eq.r * eq.b / y(eq) for eq in out]
    deltark = [(eq.r - e_init_2.r) * eq.k / y(eq) for eq in out]
    fig = Figure(resolution = (1002/2, 500))
    ax = Axis(fig[1, 1], xlabel = L"b/y", xlabelsize = 25)
    lines!(ax, b_y, rb,  label=L"$- r  b / y$")
    lines!(ax, b_y, deltark, color = :red, linestyle = :dash, linewidth = 3, label = L"$(r - r_0)  k / y")
    band!(ax, b_y, deltark, rb, color = (:blue, 0.2))
    axislegend(position = :rb)
    fig 
end

# ## Present Value of Elasticities

# ### Partial Equilibrium

cap_s = 50   # time of policy change
cap_t = 1_000  # total period of integration

# The interest rate is fixed.

function pv_PE(;R, ies_range, μ_range, 
    T = 0.0, w = 1.0, 
    cap_s = cap_s, cap_t = cap_t, Δ = 1e-4,
    value_tol = 1e-7, policy_tol = 1e-7, pdf_tol = 1e-7
)
    r = R - 1
    v = GHH(θ = 1.0, ν = 0.2)
    a_grid = grid(; stop = 250.0, length = 200, scale = :log)

    ar1 = 0.9695
    sigmaP = sqrt(0.0384)/(1 + v.ν)
    sigmaIID = sqrt(0.0522)/(1 + v.ν)
    P, z_vals = calibration(5, 2 , ar1, sigmaP, sigmaIID)

    δ = 0.1
    crra = 5.5 
    β = 0.993
    
    lst = similar(ies_range, Any)
    
    Polyester.disable_polyester_threads() do
        p = Progress(length(ies_range))
        Threads.@threads for i in eachindex(ies_range)
            ies = ies_range[i]
            u = EZ(ies = ies, ra = crra, β = β)
            h = Household(u = u, a_grid = a_grid, v = v, P = P, z_grid = z_vals)
            ws = stationary(h; R, T, w, verbose = false, value_tol, policy_tol, pdf_tol)

            !is_pol_valid(ws.a_pol, h) && @warn "R =$R, ies = $ies. Not valid policy!"

            cache = JacobianCache(ws; R, T, w, cap_s, cap_t, ΔR = 1e-4, ΔT = 0.0) # compute cache only once

            pv_lst = map(μ_range) do μ 
                Rk = μ * (r + δ) + 1 - δ
                return (; 
                    pv = Rk > 1 ? pv_elasticities!(cache, cap_s; Rk) : nothing, 
                    r, μ, ies, crra, β, dyn_efficient = Rk > 1)
            end
            lst[i] = pv_lst
            next!(p)
        end 
    end 
    
    return lst
end 

ies_range = range(0.05, 1.5, length = 100)
μ_range = range(1.001, 3.0, length = 75)
lst_PE = pv_PE(;R = 1 + e_init.r, ies_range, μ_range);

fig_pv_PE = let lst = lst_PE
    
    flattened_lst = collect(Base.Iterators.flatten(lst))

    dyn_eff = let 
        tmp = filter(flattened_lst) do x 
            !isnothing(x.pv) && x.pv > 1
        end
        [(x.ies, x.μ) for x in tmp]
    end

    bad = let 
        tmp = filter(flattened_lst) do x 
            !isnothing(x.pv) && x.pv <= 1
        end
        [(x.ies, x.μ) for x in tmp]
    end

    dyn_ine = let 
        tmp = filter(flattened_lst) do x 
            !x.dyn_efficient
        end
        [(x.ies, x.μ) for x in tmp]
    end

    fig = Figure(fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = LaTeXString("IES"), ylabel = LaTeXString("Markup, \$\\mu\$"), 
        xlabelsize = 30, ylabelsize = 30) 
    ax.xticks = [0.05, 0.5, 1., 1.5]
    plot!(ax, [x[1] for x in dyn_eff], [x[2] for x in dyn_eff], color = (:black, 0.8), markersize = 5)
    plot!(ax, [x[1] for x in dyn_ine], [x[2] for x in dyn_ine], color = (:gray, 0.5), markersize = 5)

    scatter!(ax, [1.0], [1.4], color = :white, markersize = 35, strokecolor = :black, strokewidth= 1, marker = :star5)
    fig
end

# ### General Equilibrium

# The equilibrium interest rate is solved for each parameterization.

function pv_GE(
    μ_range;  
    ies = 1, crra = 5.5, β = 0.993,
    r_range = (-0,07, 0.00), verbose = true, cap_s = cap_s, cap_t = cap_t, Δ = 1e-4,
    value_tol = 1e-7, policy_tol = 1e-7, pdf_tol = 1e-7 
)
    h = let 
        v = GHH(θ = 1.0, ν = 0.2)

        ar1 = 0.9695
        sigmaP = sqrt(0.0384)/(1 + v.ν)
        sigmaIID = sqrt(0.0522)/(1 + v.ν)
        P, z_vals = calibration(5, 2 , ar1, sigmaP, sigmaIID)

        u = EZ(ies = ies, ra = crra, β = β)
        Household(u = u, a_grid = grid(stop = 200.0, length = 200), # Extending a_max given w normalization
            v = v, P = P, z_grid = z_vals) 
    end

    δ = 0.1
    α = 0.3

    w = 1.0  # Normalizing the wage
    T = 0.0
    ws = HouseholdWorkspace(; h, R = 1 + r_range[1], T, w)
    n = labor_supply(h; w)

    f = function (r) 
        stationary!(ws; R = 1 + r, T, w, verbose = false, value_tol, policy_tol, pdf_tol)
        aa = asset_supply(h.a_grid, ws.pdf) 
        dis =  aa / n * (1 - α) - α / (r + δ)
        return dis
    end 

    r = find_zero(f, r_range, A42(), atol = 1e-6)
    stationary!(ws; R = 1 + r, T, w, verbose = false)
    !is_pol_valid(ws.a_pol, h) && @warn "Not valid policy!"
    
    cache = JacobianCache(ws; R = 1 + r, T, w, cap_s, cap_t, ΔR = 1e-4, ΔT = 0.0) # compute cache only once

    pv_lst = map(μ_range) do μ 
        Rk = μ * (r + δ) + 1 - δ
        return (; 
            pv = Rk > 1 ? pv_elasticities!(cache, cap_s; Rk) : nothing, 
            r, μ, ies, crra, β, dyn_efficient = Rk > 1)
    end 
    return pv_lst
end 

ϵ_range = range(0.05, 1.5, length = 100)
μ_range = range(1.001, 3.0, length = 75)
lst_GE = similar(ϵ_range, Any)
Polyester.disable_polyester_threads() do
    p = Progress(length(ϵ_range))
    @time Threads.@threads for i in eachindex(collect(reverse(ϵ_range)))
        lst_GE[i] =  pv_GE(μ_range; ies = ϵ_range[i], r_range = (-0.07, 0.00))
        next!(p)
    end
end

fig_pv_GE = let lst = lst_GE
    flattened_lst = collect(Base.Iterators.flatten(lst))

    dyn_eff = let 
        tmp = filter(flattened_lst) do x 
            !isnothing(x.pv) && x.pv > 1
        end
        [(x.ies, x.μ) for x in tmp]
    end

    bad = let 
        tmp = filter(flattened_lst) do x 
            !isnothing(x.pv) && x.pv <= 1
        end
        [(x.ies, x.μ) for x in tmp]
    end

    dyn_ine = let 
        tmp = filter(flattened_lst) do x 
            !x.dyn_efficient
        end
        [(x.ies, x.μ) for x in tmp]
    end

    fig = Figure(fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = LaTeXString("IES"), ylabel = LaTeXString("Markup, \$\\mu\$"), 
        xlabelsize = 30, ylabelsize = 30) 
    ax.xticks = [0.05, 0.5, 1., 1.5]
    plot!(ax, [x[1] for x in dyn_eff], [x[2] for x in dyn_eff], color = (:black, 0.8), markersize = 5)
    plot!(ax, [x[1] for x in dyn_ine], [x[2] for x in dyn_ine], color = (:gray, 0.5), markersize = 5)
    
    scatter!(ax, [1.0], [1.4], color = :white, markersize = 35, strokecolor = :black, strokewidth= 1, marker = :star5)

    fig 
end 

# ## Saving the Figures and Statistics

save(joinpath(@__DIR__, "..", "output", "figures", "transition_efficient_fixed_k.pdf"), f1)
save(joinpath(@__DIR__, "..", "output", "figures", "transition_efficient_golden_k.pdf"), f2)
save(joinpath(@__DIR__, "..", "output", "figures", "transition_efficient_no_debt.pdf"), f3)
save(joinpath(@__DIR__, "..", "output", "figures", "steady_state_transfers.pdf"), f4)

save(joinpath(@__DIR__, "..", "output", "figures", "pv_elasticities_regions_PE.pdf"), fig_pv_PE);
save(joinpath(@__DIR__, "..", "output", "figures", "pv_elasticities_regions_GE.pdf"), fig_pv_GE);

open(joinpath(@__DIR__, "..", "output", "tables", "statistics.txt"), "w") do file
    write(file, statistics)
end
