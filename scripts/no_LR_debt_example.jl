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
_ITERS = 50

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

T1 = 100  # period of increase in debt
T2= 200 #period of decrease in debt
T3= 0 # period of increasing k
H = 50   # stationary part


# ## Transition with increasing k and hump-shaped b

# # Smooth debt policy
b_target = laissez_faire.y * 0.60
ρB = 0.9
b_list_3 = Array{Float64,1}(undef, T1+T2 + H)
b_list_3[1] = 0.0
b_list_3[2] = laissez_faire.y * 0.05
for i in 3:T1-1
    b_list_3[i] = b_list_3[2] * ρB^(i-2) + (1 - ρB^(i-2)) * b_target
end
b_list_3[T1]=b_target
for i in T1+1:T1+T2-1
    b_list_3[i] =  b_target*(1-(i-T1)/T2)
end
b_list_3[T1+T2:end] .= 0.0

#Linear debt policy
# b_target = laissez_faire.y * 0.60
# b_list_3 = Array{Float64,1}(undef, T1+T2 + H)
# b_list_3[1] = 0.0
# for i in 2:T1
#     b_list_3[i] =  0.0+b_target*i/T1
# end

# for i in T1+1:T1+T2-1
#     b_list_3[i] =  b_target*(1-(i-T1)/T2)
# end
# b_list_3[T1+T2:end] .= 0.0


# We go to the golden rule k .

k_golden = golden_rule_k(e.t, laissez_faire.n)
ρK = 0.95
k_list_3 = similar(b_list_3)
for i in eachindex(k_list_3)
    k_list_3[i] = laissez_faire.k * ρK^(i-1) + (1 - ρK^(i-1)) * k_golden
end


# Computing the new stationary equilibrium in this case:

@time final_eq_3 = solve_new_stationary_equilibrium_given_k_b(
    laissez_faire,
    k_list_3[end],
    0.0;
    r_range = (laissez_faire.r, -0.010),
    tol = (value_function = 1e-10, distribution = 1e-13)
)

# Computing the transition:

r_path_3 = nothing
if _LOAD_GUESSES
    r_path_3 = try
        # load the transfer vector from previous iterations
        readdlm(joinpath(@__DIR__, "..", "output", "tmp_calcs", "tmp_r_003.txt"))
    catch b
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

# _SAVE_GUESSES && open(joinpath(@__DIR__, "..", "output", "tmp_calcs", "tmp_r_002.txt"), "w") do io
#     writedlm(io, transition_2.path.r)
# end

# ### Plots

#plot(transition_3.path.transfer, fill = 0, alpha = 0.3)
plot(transition_3.path.transfer, color = 1)
hline!([0, final_eq_3.transfer], lw = 1)
#vline!([T1,T1+T2,T3], lw = 1)

savefig(joinpath(@__DIR__, "output", "figures", "transfers_no_LR_debt.pdf"))


plot(transition_3.path.r, legend=false)
hline!([laissez_faire.r], lw = 1)
hline!([final_eq_3.r])

savefig(joinpath(@__DIR__, "output", "figures", "interest_no_LR_debt.pdf"))
