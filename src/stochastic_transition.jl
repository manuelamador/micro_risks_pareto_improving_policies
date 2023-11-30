# This file contains the code for the stochastic transition model. In this case, there is an initial period, t = 0, where households decide on capital and debt holdings. In period t= 1, an aggregate productivity shock (with a deterministic path of productivity for all subsequent t) is realized and the economy transitions to a new stationary equilibrium. 


Base.@kwdef struct TransitionTree{T1, T2, T3, T4, T5}
    h::T1
    first::T2 
    paths::T3 
    prob::T4
    leaves::T5
end 


function initialize_tree(; h, t, path_lengths = [100, 100], ρ = 1.0, Avec = [1.05, 0.95], prob = [0.5, 0.5], r_range = nothing, init_pdf = nothing)

    isnothing(r_range) && (r_range = [nextfloat(- t.δ), 1 / h.u.β] )  
    isnothing(init_pdf) && (init_pdf = stationary_laissez_faire(h, t; r_range).ws.pdf) 

    paths = [
        [_transition_item(
                h, CobbDouglasTechnology(; 
                    t.α, t.δ, t.μ, 
                    A = (t.A * (a1 - 1) * ρ^(i-1) + t.A)
                )
            ) for i in 1:path_length] 
        for (path_length, a1) in zip(path_lengths, Avec)
    ]
    first = _transition_item(h, t)
    first.pdf .= init_pdf

    # compute equilibria at the end points of the paths
    leaves = [
        stationary_laissez_faire(h, path[end].t; r_range) for path in paths
    ] 

    return TransitionTree(; h, first, paths, prob, leaves)
end


function laissez_faire_transition(; tree::TransitionTree, verbose = true, nlsolve_kwargs = nothing)

    h = tree.h 

    nlsolve_baseline_kwargs = (
        ftol = _ZERO_FTOL,
        show_trace = verbose,
        method = :anderson,
        m = 10,
        iterations = _ZERO_MAX_ITERS,
        beta = -0.1  # adjustment parameter for the fixed point algorithm
    )
    
    nlsolve_kwargs_merged = isnothing(nlsolve_kwargs) ?             
        nlsolve_baseline_kwargs :
        merge(nlsolve_baseline_kwargs, nlsolve_kwargs)
    

    K0 = aggregate_assets(h.a_grid, tree.first.pdf)  # initial distribution of wealth -- get the initial capital stock 

    # Creating a vector of capital guesses.  
    # The first element is K1 (the capital choice at t = 0)
    # The next elements stack the choices for the paths one after the other. \
    # We start with the guess equal to the initial capital level. 
    K_guess = [K0]
    residuals = [zero(K0)]
    for path in tree.paths
        for _ in path[2:end] 
            push!(K_guess, K0)
            push!(residuals, zero(K0)) 
        end 
    end

    normalized_N = aggregate_labor(h; w = 1)
    assign_given_K_LF!(; node = tree.first, K = K0, normalized_N) 
    laissez_faire_transition_residuals!(residuals, K_guess, tree; h)


    f! = (F, x) -> laissez_faire_transition_residuals!(F, x, tree; h)
    k_sol = nlsolve(f!, K_guess; nlsolve_kwargs_merged...)
    laissez_faire_transition_residuals!(similar(k_sol.zero), k_sol.zero, tree; h)
    
    return tree 
end


get_labor_given_K_LF(; h, t, K, normalized_N) = get_labor_given_K_LF(h.v, t, K, normalized_N)
get_labor_given_K_LF(v::GHH, t::CobbDouglasTechnology, K, normalized_N) = (normalized_N * (1 - t.α)^v.ν)^(1/(1 + t.α * (v.ν))) * (t.A / t.μ * K^t.α) ^(v.ν / (1 + t.α * v.ν))
get_labor_given_K_LF(::FixedLabor, _, _, normalized_N) = normalized_N


function assign_given_K_LF!(; node, K, normalized_N)
    node.k = K # the initial capital in the path is K1 across all paths
    node.a = K
    node.n = get_labor_given_K_LF(; node.h, node.t, K, normalized_N) # set the labor supply
    node.rk = mpk(node.t; node.k, node.n) / node.t.μ  # return to capital 
    node.w = mpl_from_mpk(node.t; mpk = node.rk * node.t.μ) / node.t.μ
    node.T = 0.0 
    node.b = 0.0 
end 


function laissez_faire_transition_residuals!(residuals, K_guess, tree::TransitionTree; h, normalized_N = aggregate_labor(h; w = 1))

    leaves = tree.leaves

    guess_i = 2
    for path in tree.paths
        assign_given_K_LF!(; node = path[1], K = K_guess[1], normalized_N)
        for i in firstindex(path)+1:lastindex(path)
            assign_given_K_LF!(; node = path[i], K = K_guess[guess_i], normalized_N)
            guess_i += 1
        end 
    end

    # BACKWARD ITERATION 

    for i_agg in eachindex(tree.paths) 
        (path, leave) = (tree.paths[i_agg], leaves[i_agg]) 
        for i in reverse(eachindex(path))
            node = path[i]
            next_values = i == lastindex(path) ? leave.ws : path[i+1]
            R = 1 + node.rk - node.t.δ
            backwards_once!(h, node; next_values = next_values, R, node.T, node.w, node.a_tmp)
            asset_policy_given_η!(node.a_pol; h, node.η, R, node.T, node.w)
            interpolate_asset_policy!(node.lower_index, node.lower_weight; h.a_grid, node.a_pol)    
        end 
    end 

    # the time 0 policies. 
    let node = tree.first, R = 1 + node.rk - node.t.δ
        backwards_once!(h, node; next_values = tree, R , node.T, node.w, node.a_tmp)
        asset_policy_given_η!(node.a_pol; h, node.η, R, node.T, node.w)
        interpolate_asset_policy!(node.lower_index, node.lower_weight; h.a_grid, node.a_pol)    
    end 

    # FORWARD ITERATION 

    # the time 0 policy
    forward_pdf!(tree.paths[1][1].pdf; h, tree.first.pdf, tree.first.lower_index, tree.first.lower_weight)
    a1 = aggregate_assets(h.a_grid, tree.paths[1][1].pdf)
    tree.paths[1][1].a = a1
    for i_agg in 2:length(tree.paths)
        tree.paths[i_agg][1].pdf .= tree.paths[1][1].pdf
        tree.paths[i_agg][1].a = a1
    end 

    residuals[1] = a1 - tree.paths[1][1].k

    # forward iteration of the distribution and getting the aggregates

    residual_index = 2
    for i_agg in eachindex(tree.paths)
        path = tree.paths[i_agg]
        for i in firstindex(path):length(path)-1
            forward_pdf!(path[i+1].pdf; h, path[i].pdf, path[i].lower_index, path[i].lower_weight)
            a = aggregate_assets(h.a_grid, path[i+1].pdf)
            path[i+1].a = a
            residuals[residual_index] =  a - path[i+1].k
            residual_index += 1
        end
    end 

    residuals
end 



# Specializing the backward_euler functions for t = 0 in the stochastic case. 

function backwards_euler_x(u::CRRA, P, i, s, tree::TransitionTree)
    ξ = get_inverse_ies(u)
    β = get_β(u)
    Evx = 0.0
    for (p, next_values) in zip(tree.prob, tree.paths)
        # ^ we need to iterate across aggregate realizations tomorrow  
        for j in axes(P, 2)
            Evx += p * P[s, j] * (next_values[1].η[i, j])^(-ξ)
        end 
    end 
    x = (β * Evx)^(-1/ξ)    
    return x
end


function backwards_euler_x(u::EZ, P, i, s, tree::TransitionTree)
    ξ = get_inverse_ies(u)
    γ = get_ra(u)
    β = get_β(u)
    Ev = Evx = 0.0
    for (p, next_values) in zip(tree.prob, tree.paths)
        # ^ we need to iterate across aggregate realizations tomorrow 
        for j in axes(P, 2)
            Evx += p * P[s, j] * (next_values[1].v[i, j])^((ξ - γ)) * (next_values[1].η[i, j])^(-ξ)
            Ev += p * P[s, j] * u.risk(next_values[1].v[i, j])
        end
    end 
    x = (β * inverse(u.risk, Ev^((γ - ξ))) * Evx)^(-1/ξ)
    return x
end


function _backward_once_helper!(u::EZ, current_values, tree::TransitionTree, h, a_prime, P, x, i, s
)
    β = u.β

    Ev = zero(a_prime)
    for (p, next_values) in zip(tree.prob, tree.paths)
        # ^ we need to iterate across aggregate realizations tomorrow 
        for s2 in axes(next_values[1].v, 2)
            v_prime = interp1D(a_prime, h.a_grid, view(next_values[1].v, :, s2))
            Ev += p * P[s, s2] * u.risk(v_prime)
        end
    end 
    current_values.v[i, s] =
        inverse(u.temporal, (1 - β) * u.temporal(x) + β * u.temporal(inverse(u.risk, Ev)))
end 
