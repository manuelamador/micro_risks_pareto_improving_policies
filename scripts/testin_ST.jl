


using MicroRisks 
using CairoMakie

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

# # Transition to a higher debt level

# # b = 0.6 y0 
# b_target = y(e_init) * 0.60

# # Final equilibrium with higher debt and same k

# @time e_final = stationary_equilibrium_given_k_b(e_init, e_init.k, b_target; r_range = (-0.02, 0.0), verbose = true)


# # Debt policy and capital (constant) along the transition

# b_path = let   # Smooth path of increasing debt  
#     T = 100  # period of adjustment of debt
#     H = 50   # debt level no longer moving
#     ρB = 0.9
#     b_list = Array{Float64,1}(undef, T + H)
#     b_list[1] = 0.0
#     b_list[2] = y(e_init) * 0.05
#     b_list[T:end] .= b_target
#     for i in 3:T-1
#         b_list[i] = b_list[2] * ρB^(i-2) + (1 - ρB^(i-2)) * b_target
#     end
#     b_list
# end;

# k_path = [e_init.k for _ in b_path];

# # Solving the transition 

# @time path = MicroRisks.solve_RPI_transition(e_init, e_final; k_path, b_path);


# Some simple simulations 

tree = MicroRisks.initialize_tree(;h, t, Avec = [1.02, 0.98], ρ = 0.9, init_pdf = nothing)
out = MicroRisks.laissez_faire_transition(; tree)

fig = Figure()

axis1 = fig[1, 1]
axis2 = fig[1, 2]
axis3 = fig[2, 2]

lines(axis1, [node.t.A for node in out.paths[1]])
lines!(axis1,  [node.t.A for node in out.paths[2]])

lines(axis2,  [node.k for node in out.paths[1]])
lines!(axis2,  [node.k for node in out.paths[2]])
hlines!(axis2, [e_init.k])

lines(axis3, [node.rk - node.t.δ for node in out.paths[1]])
lines!(axis3,  [node.rk - node.t.δ for node in out.paths[2]])

fig
