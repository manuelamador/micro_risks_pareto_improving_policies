using Revise
using Aiyagari
using StructArrays
using ProgressMeter
using DelimitedFiles
using Roots
using Plots
using LaTeXStrings
using PrettyTables

pgfplotsx();
Threads.nthreads()

ProgressMeter.ijulia_behavior(:clear);
default(label = "", lw = 2, dpi = 300, left_margin = 0Plots.mm, format=:svg);

# +
_LOAD_GUESSES = true # load the initial starting guesses for zeros from disk
_SAVE_GUESSES = false # save zeros to disk
_ITERS = 50

#preferences (first element is benchmark)
ies_vals=(1.0, 0.5, 1.5)
crra_vals=(5.5, 2.0, 10.0)
β_vals=(0.993, 0.97, 0.98)

e_vals = NamedTuple{(:ies, :crra, :β, :e),Tuple{Float64,Float64,Float64,Economy}}[]

for ies in ies_vals
    crra=crra_vals[1]
    β=β_vals[1]
            e = let
                #Household:
                ar1 = 0.9695
                sigmaP = sqrt(0.0384)/(1.2)
                sigmaIID = sqrt(0.0522)/(1.2)
                P, z_vals = calibration(5, 2 , ar1, sigmaP, sigmaIID)
            
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
            push!(e_vals, (ies=ies, crra=crra, β = β, e = e))
end


for crra in crra_vals[2:end]
    ies=ies_vals[1]
    β=β_vals[1]
            e = let
                #Household:
                ar1 = 0.9695
                sigmaP = sqrt(0.0384)/(1.2)
                sigmaIID = sqrt(0.0522)/(1.2)
                P, z_vals = calibration(5, 2 , ar1, sigmaP, sigmaIID)
            
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
            push!(e_vals, (ies=ies, crra=crra, β = β, e = e))
end

for β  in β_vals[2:end]    
    ies=ies_vals[1]
    crra=crra_vals[1]
            e = let
                #Household:
                ar1 = 0.9695
                sigmaP = sqrt(0.0384)/(1.2)
                sigmaIID = sqrt(0.0522)/(1.2)
                P, z_vals = calibration(5, 2 , ar1, sigmaP, sigmaIID)
            
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
            push!(e_vals, (ies=ies, crra=crra, β = β, e = e))
end



# ## Transitions

T = 100  # period of adjustment
H = 50   # stationary part


#Iterate over preferences

transitions_vals=[]

Threads.@threads for i=1:length(e_vals)
    println(i)

    e=e_vals[i].e
    
    # Solve laissez-faire economy
    @time laissez_faire = let r_range = (-0.04, -0.00)
        solve_laissez_faire(e;
            r_range = r_range,
            tol =  (value_function = 1e-10, distribution = 1e-13),verbose=false
        )
    end

    # We target a debt of 59% percent of output

    b_target = laissez_faire.y * 0.60
    k_target = laissez_faire.k


    # Solve final SS
    @time final_eq_1 = solve_new_stationary_equilibrium_given_k_b(
        laissez_faire,
        k_target,
        b_target;
        r_range = (laissez_faire.r-.01, laissez_faire.r+.01),
        tol = (value_function = 1e-10, distribution = 1e-13),verbose=false
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

    push!(transitions_vals,(e=e,laissez_faire=laissez_faire,transition=transition))
end

out=NamedTuple{(:ies, :crra, :β,  :y0, :e, :path_r, :path_tr),Tuple{Float64,Float64,Float64,Float64,Economy,Array{Float64, 1},Array{Float64, 1}}}[]

for i=1:length(transitions_vals)
    e=transitions_vals[i].e
    ies=e.h.u.pars[1]
    crra=e.h.u.pars[2]
    β=e.h.β
    y0=transitions_vals[i].laissez_faire.y
    path_r=transitions_vals[i].transition.path.r
    path_tr=transitions_vals[i].transition.path.transfer

    push!(out,(ies=ies,crra=crra,β=β,y0=y0,e=e,path_r=path_r,path_tr=path_tr))
end

benchmark=[y for y in out if y.ies==1.0 && y.crra==5.5 && y.β==0.993]
benchmark=benchmark[1]

#IES
plot(benchmark.path_r,  linecolor=:blue, lw = 1)
plot!([y.path_r for y in out if y.ies==0.5 && y.crra==5.5 && y.β==0.993], linecolor=:black, linestyle=:dash, lw = 1)
plot!([y.path_r for y in out if y.ies==1.5 && y.crra==5.5 && y.β==0.993], linecolor=:red, linestyle=:dot, lw = 1)

savefig(joinpath(@__DIR__, "output", "figures", "IES_interest.pdf"))


plot(benchmark.path_tr./benchmark.y0,  linecolor=:blue, lw = 1)
plot!([y.path_tr./y.y0 for y in out if y.ies==0.5 && y.crra==5.5 && y.β==0.993], linecolor=:black, linestyle=:dash, lw = 1)
plot!([y.path_tr./y.y0 for y in out if y.ies==1.5 && y.crra==5.5 && y.β==0.993], linecolor=:red, linestyle=:dot, lw = 1)

savefig(joinpath(@__DIR__, "output", "figures", "IES_transfers.pdf"))


#Discount
plot(benchmark.path_r,  linecolor=:blue, lw = 1)
plot!([y.path_r for y in out if y.ies==1.0 && y.crra==5.5 && y.β==0.97], linecolor=:black, linestyle=:dash, lw = 1)
plot!([y.path_r for y in out if y.ies==1.0 && y.crra==5.5 && y.β==0.98], linecolor=:red, linestyle=:dot, lw = 1)
savefig(joinpath(@__DIR__, "output", "figures", "Beta_interest.pdf"))

plot(benchmark.path_tr./benchmark.y0,  linecolor=:blue, lw = 1)
plot!([y.path_tr./y.y0 for y in out if y.ies==1.0 && y.crra==5.5 && y.β==0.97], linecolor=:black, linestyle=:dash, lw = 1)
plot!([y.path_tr./y.y0 for y in out if y.ies==1.0 && y.crra==5.5 && y.β==0.98], linecolor=:red, linestyle=:dot, lw = 1)
savefig(joinpath(@__DIR__, "output", "figures", "Beta_transfers.pdf"))

#CRRA
plot(benchmark.path_r,  linecolor=:blue, lw = 1)
plot!([y.path_r for y in out if y.ies==1.0 && y.crra==2.0 && y.β==0.993], linecolor=:black, linestyle=:dash, lw = 1)
plot!([y.path_r for y in out if y.ies==1.0 && y.crra==10.0 && y.β==0.993], linecolor=:red, linestyle=:dot, lw = 1)
savefig(joinpath(@__DIR__, "output", "figures", "CRRA_interest.pdf"))

plot(benchmark.path_tr./benchmark.y0,  linecolor=:blue, lw = 1)
plot!([y.path_tr./y.y0 for y in out if y.ies==1.0 && y.crra==2.0 && y.β==0.993], linecolor=:black, linestyle=:dash, lw = 1)
plot!([y.path_tr./y.y0 for y in out if y.ies==1.0 && y.crra==10.0 && y.β==0.993], linecolor=:red, linestyle=:dot, lw = 1)
savefig(joinpath(@__DIR__, "output", "figures", "CRRA_transfers.pdf"))


#numbers

rates=[]
push!(rates,(Elasticity="Benchmark", Initial=benchmark.path_r[1],Final=benchmark.path_r[end],Max=maximum(benchmark.path_r)))
push!(rates,
    (Elasticity="Inelastic", Initial=[y.path_r for y in out if y.ies==0.5 && y.crra==5.5 && y.β==0.993][1][1],
    Final=[y.path_r for y in out if y.ies==0.5 && y.crra==5.5 && y.β==0.993][1][end], 
    Max=maximum([y.path_r for y in out if y.ies==0.5 && y.crra==5.5 && y.β==0.993][1])
    )
)
push!(rates,
    (Elasticity="Elastic", Initial=[y.path_r for y in out if y.ies==1.5 && y.crra==5.5 && y.β==0.993][1][1],
    Final=[y.path_r for y in out if y.ies==1.5 && y.crra==5.5 && y.β==0.993][1][end], 
    Max=maximum([y.path_r for y in out if y.ies==1.5 && y.crra==5.5 && y.β==0.993][1])
    )
)


transfers=[]
push!(transfers,(Elasticity="Benchmark", y0=benchmark.y0, Initial=benchmark.path_tr[1],Final=benchmark.path_tr[end],Min=minimum(benchmark.path_tr)))
push!(transfers,
    (Elasticity="Inelastic", y0=[y.y0 for y in out if y.ies==0.5 && y.crra==5.5 && y.β==0.993], Initial=[y.path_tr./y.y0 for y in out if y.ies==0.5 && y.crra==5.5 && y.β==0.993][1][1],
    Final=[y.path_tr./y.y0 for y in out if y.ies==0.5 && y.crra==5.5 && y.β==0.993][1][end], 
    Min=minimum([y.path_tr./y.y0 for y in out if y.ies==0.5 && y.crra==5.5 && y.β==0.993][1])
    )
)
push!(transfers,
    (Elasticity="Elastic", y0=[y.y0 for y in out if y.ies==1.5 && y.crra==5.5 && y.β==0.993], Initial=[y.path_tr./y.y0 for y in out if y.ies==1.5 && y.crra==5.5 && y.β==0.993][1][1],
    Final=[y.path_tr./y.y0 for y in out if y.ies==1.5 && y.crra==5.5 && y.β==0.993][1][end], 
    Min=minimum([y.path_tr./y.y0 for y in out if y.ies==1.5 && y.crra==5.5 && y.β==0.993][1])
    )
)



pretty_table([y for y in rates])

pretty_table([y for y in transfers])

