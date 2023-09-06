#############################################################################
# DETERMINISTIC TRANSITION PATHS
#############################################################################


"""
   _transition_item(h::Household)
   
Creates a workspace for the transition path at one date.
"""
function _transition_item(h::Household)  
    (; v, η, a_pol, pdf, lower_index, lower_weight, a_tmp) = _generate_base_workspace_matrices(h)

    w = r = T = a = k = b = zero(eltype(h.a_grid)) 
    return (; h, v, η, a_pol, pdf, lower_index, lower_weight, a_tmp, r, w, T, a, k, b)
end 


"""
    _initialize_path(e_init, e_final; k_path)

Initializes the transition path from `e_init` to `e_final` with the given `k_path`. The length of the transition path 
is the length of `k_path`.
"""
function _initialize_path(e_init, e_final; k_path)
    h = e_init.h

    path = StructArray(_transition_item(h) for _ in k_path)
    path[1].pdf .= e_init.ws.pdf  # initial distribution
 
    path[end].v .= e_final.ws.v 
    path[end].η .= e_final.ws.η
    path[end].a_pol .= e_final.ws.a_pol
    path[end].lower_index .= e_final.ws.lower_index
    path[end].lower_weight .= e_final.ws.lower_weight

    path.r[end] = e_final.r
    path.b[end] = e_final.b
    path.k[end] = e_final.k
    path.a[end] = e_final.a

    path.w .= e_init.w  # the wage is constant along the path. 

    return path 
end 


"""
    solve_RPI_transition_robust(e_init, e_final; k_path, b_path, r_path, verbose, nlsolve_kwargs)

Solves for the RPI transition (constant wage/profits) path from a stationary equilibrium `e_init` to a sationary equilibrium `e_final` with the given `k_path`, `b_path`. The function uses a robust but potentially slow algorithm to solve for the transition path.

# Arguments
- `e_init`: Initial equilibrium.
- `e_final`: Final equilibrium.
- `k_path`: Path of capital levels.
- `b_path`: Path of government debt levels.
- `r_path`: Path of interest rates. If `nothing`, the function assumes that the interest rate is constant at the final equilibrium value.
- `verbose`: Whether to show progress.
- `nlsolve_kwargs`: Keyword arguments to pass to the `nlsolve` solver.
"""
function solve_RPI_transition_robust(
    e_init, e_final; k_path, b_path, r_path = nothing, verbose = true,
    nlsolve_kwargs = nothing
)
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

    path = _initialize_path(e_init, e_final; k_path)

    if isnothing(r_path) 
        r_path_rest = fill(e_final.r, length(k_path) - 1)
    else 
        r_path_rest = r_path[2:length(k_path)]
    end 

    f! = (F, x) -> _RPI_transition_residuals(F, x; path, e_init, e_final, k_path, b_path)
    r_sol = nlsolve(f!, r_path_rest; nlsolve_kwargs_merged...)
    _RPI_transition_residuals(similar(r_sol.zero), r_sol.zero; path, e_init, e_final, k_path, b_path)
    
    return path 
end


"""
    solve_RPI_transition(e_init, e_final; k_path, b_path, r_path, verbose, max_iters, tol)

Solves for the RPI transition (constant wage/profits) path from `e_init` to `e_final` with the given `k_path`, `b_path`. It uses the Jacobian to accelerate the root finding algorithm. 

# Arguments
- `e_init`: Initial equilibrium.
- `e_final`: Final equilibrium.
- `k_path`: Path of capital levels.
- `b_path`: Path of government debt levels.
- `r_path`: Path of interest rates. If `nothing`, the function assumes that the interest rate is constant at the final equilibrium value.
- `verbose`: Whether to show progress.
- `max_iters`: Maximum number of iterations.
- `tol`: Tolerance for convergence.
"""
function solve_RPI_transition(
    e_init, e_final; k_path, b_path, r_path = nothing, verbose = true, 
    max_iters = _ZERO_MAX_ITERS, tol = _ZERO_FTOL
)    
    path = _initialize_path(e_init, e_final; k_path)

    if isnothing(r_path) 
        r_path_rest = fill(e_final.r, length(k_path) - 1)
    else 
        r_path_rest = r_path[2:length(k_path)]
    end 

    jac_R = jacobian(e_final; cap_t = length(r_path_rest) + 1, ΔR = 1e-4, ΔT = 0.0)
    jac_T = jacobian(e_final; cap_t = length(r_path_rest) + 1, ΔT = 1e-4, ΔR = 0.0)
    jac = similar(jac_T)
    for t in axes(jac, 2)
        for s in axes(jac, 1)
            jac[t, s] = jac_R[t, s] - jac_T[t, s] * (b_path[end] + k_path[end])
        end
    end 
    jac = jac[2:end, 2:end]

    cond(jac) > 1e+8 && @warn("Jacobian seems ill-conditioned. If this doesn't converge, try solve_RPI_transition_robust.") 
    Afact = factorize(jac)

    residuals = fill(zero(e_final.a), length(r_path_rest))
    f! = (F, x) -> _RPI_transition_residuals(F, x; path, e_init, e_final, k_path, b_path)

    iter = 0
    p = ProgressUnknown()
    while true 
        f!(residuals, r_path_rest)
        dis = maximum((abs(x) for x in residuals))
        verbose && next!(p,  showvalues = [(:error, dis)])
        dis < tol && (finish!(p, showvalues = [(:error, dis)]); break)
        iter += 1
        iter >= max_iters && (@warn("Did not converge after $max_iters iterations"); break)
        r_path_rest .= r_path_rest .-  Afact \ residuals
    end  

    return path 
end


"""
    _RPI_transition_residuals(residuals, r_path_rest; path, e_init, e_final, k_path, b_path)

Calculates the residuals of the transition path given the interest rate path `r_path_rest`. 
It modifies `residuals` in place with the excess asset supply at each date.
It also modifies the `path` workspace in place when doing the calculations.
"""
function _RPI_transition_residuals(residuals, r_path_rest; path, e_init, e_final, k_path, b_path)
    (; h, t, w) = e_init
    n0, r0, k0 = e_init.n, e_init.r, e_init.k

    r_min = - t.δ + eps()  # minimum interest rate
    r_max = 1 / h.u.β  # a maximum -- could be relaxed
    min_T = nextfloat(minimum_feasible_transfer(h, w)) # minimum transfer
    
    # backward iteration solving for the optimal policies
    for i in reverse(eachindex(path, b_path, k_path))
        b, k = b_path[i], k_path[i]
        bprime = i == lastindex(b_path) ? e_final.b : b_path[i + 1] 
        r = (i == firstindex(path)) ? r0 : r_path_rest[i - 1]   # start from the initial interest rate r0
        r = clamp(r, r_min, r_max)
        R = 1 + r

        T = get_T(t; b, bprime, r, k, r0, k0, n0)
        T = max(T, min_T) 
        
        path.T[i] = T
        path.r[i] = r
        next_values = i == lastindex(path) ? e_final.ws : path[i+1]
        backwards_once!(h, path[i]; next_values = next_values, R, T, w, path[i].a_tmp)
        asset_policy_given_η!(path[i].a_pol; h, path[i].η, R, T, w)
        interpolate_asset_policy!(path[i].lower_index, path[i].lower_weight; h.a_grid, path[i].a_pol)
    end

    # forward iteration of the distribution and getting the aggregates 
    for i in eachindex(path) 
        a = aggregate_assets(h.a_grid, path[i].pdf)
        path.a[i] = a
        path.b[i] = b_path[i]
        path.k[i] = k_path[i]
        if i > firstindex(path)
            residuals[i - 1] =  a - k_path[i] - b_path[i]
        end 
        (i == lastindex(path)) && break
        forward_pdf!(path[i+1].pdf; h, path[i].pdf, path[i].lower_index, path[i].lower_weight)
    end
end 

