# This file contains the code for solving for the stationary equilibria of the economy, both under laissez-faire and in the presence of taxes and debt. 



"""
    stationary_equilibrium(h, t; r_range[, verbose, hh_problem_kwargs, solver, find_zero_kwargs])

Solve for the stationary equilibrium of the economy with the given household problem and technology with no taxes or debt (laissez_faire). 

# Arguments
- `h`: Household struct
- `t`: Technology struct
- `r_range`: Range of interest rates to search over 
- `hh_problem_kwargs`: Keyword arguments to pass to the household problem solver.
- `find_zero_kwargs`: Keyword arguments to pass to the `find_zero` solver.
- `solver`: Solver to use for finding the interest rate.
- `verbose`: Whether to show progress of the solver.
"""
function stationary_laissez_faire(h, t; r_range, 
    verbose = true, 
    hh_problem_kwargs = nothing,
    find_zero_kwargs = (;atol = _ZERO_FTOL),
    solver = A42()
)

    hh_problem_options_baseline = (; verbose = false) 
    stationary_kwargs = isnothing(hh_problem_kwargs) ?
            hh_problem_options_baseline : 
            merge(hh_problem_options_baseline, hh_problem_kwargs)
    
    r0 = r_range[1]
    w = mpl_from_mpk(t; mpk = r0 + t.δ)
    pars = (R = 1 + r0, T = zero(r0), w = w)
    ws = HouseholdWorkspace(; h, pars...)

    p = ProgressUnknown()
    f = (x) -> _residual_stationary_laissez_faire!(ws, x, h, t, verbose, p, stationary_kwargs).residual 
    r = find_zero(f, r_range, solver; find_zero_kwargs...)
    # get the rest of the equilibrium objects
    (; w, n, k, a) = _residual_stationary_laissez_faire!(ws, r, h, t, false, p, stationary_kwargs)
    verbose && finish!(p, showvalues = [(:r, r), (:error, a - k)])

    e = StationaryEquilibrium(;
        h, t, ws, w, n, k, a, r, 
        T = zero(r), b = zero(r)
    )
    return e
end


"""
    _residual_stationary_laissez_faire!(ws, r, h, t; verbose, p, stationary_kwargs)

Get the residual of the stationary equilibrium given the interest rate `r` and the household and technology structs `h` and `t`.
It modifies the `ws` workspace in place when doing the calculations. 
"""
function _residual_stationary_laissez_faire!(ws, r, h, t, verbose, p, stationary_kwargs)
    w = mpl_from_mpk(t; mpk = r + t.δ)
    pars = (R = 1 + r, T = zero(r), w)
    stationary!(ws; stationary_kwargs..., pars...)
    a = aggregate_assets(h.a_grid, ws.pdf)
    n = aggregate_labor(h; w = w)
    k = k_from_mpk_n(t; mpk = r + t.δ, n)
    residual = a - k
    verbose && next!(p, showvalues = [(:r, r), (:error, residual)])
    return (; residual, w, a, n, k)
end


"""
    stationary_equilibrium_given_k_b(e_init, k, b; r_range[, verbose, hh_problem_kwargs, solver, find_zero_kwargs])

Solves for the stationary equilibrium of the economy with the given household problem and technology with a given level of capital and debt. It starts from the equilibrium described in `e_init` and solves for the interest rate `r` that solves the household problem, making sure that the government budget constraint holds with equality. In this code, it is assumed that the policy maintains the wage at the level of the initial equilibrium. 

# Arguments 
- `e_init`: Initial laissez-faire equilibrium.
- `k`: Capital level to solve for.
- `b`: Debt level to solve for.
- `r_range`: Range of interest rates to search over.
- `hh_problem_kwargs`: Keyword arguments to pass to the household problem solver.
- `find_zero_kwargs`: Keyword arguments to pass to the `find_zero` solver.
- `solver`: Solver to use for finding the interest rate.
- `verbose`: Whether to show progress of the solver.
"""
function stationary_equilibrium_given_k_b(e_init, k, b; r_range, 
    verbose = true,   
    hh_problem_kwargs = nothing,
    solver = A42(), 
    find_zero_kwargs = (; atol = _ZERO_FTOL)
)
    hh_problem_options_baseline = (; 
        verbose = false
    ) 

    stationary_kwargs = isnothing(hh_problem_kwargs) ?
        hh_problem_options_baseline :
        merge(hh_problem_options_baseline, hh_problem_kwargs)

    (;h, t, w) = e_init
    
    k0, n0, r0, w0 = e_init.k, e_init.n, e_init.r, e_init.w
    r = r_range[1]
    T = get_T(t; b, bprime = b, r, k, r0, k0, n0)
    pars = (R = 1 + r, T = T, w = w0)

    ws = HouseholdWorkspace(; h, pars...)
    min_T = minimum_feasible_transfer(h, w0)

    p = ProgressUnknown()
    f = (x) -> _k_b_residual!(ws, x, e_init, k, b, min_T, verbose, p, stationary_kwargs).residual
    r = find_zero(f, r_range, solver; find_zero_kwargs...)
    (; a, T) = _k_b_residual!(ws, r, e_init, k, b, min_T, false, p, stationary_kwargs)
    verbose && finish!(p, showvalues = [(:r, r), (:error, a - k - b)])

    e = StationaryEquilibrium(; h, t, ws, w = w0, n = n0, k, a, r, b, T)
    return e
end


"""
    _k_b_residual!(ws, r, e_init, k, b, min_T, verbose, p, stationary_kwargs)

Get the residual of the stationary equilibrium given `k` and `b`, for the interest rate `r` and the household and technology structs `h` and `t`.
Modifies the `ws` workspace in place when doing the calculations.
"""
function _k_b_residual!(ws, r, e_init, k, b, min_T, verbose, p, stationary_kwargs)
    (; h, t) = e_init
    k0, n0, r0, w0 = e_init.k, e_init.n, e_init.r, e_init.w 
    T = get_T(t; b, bprime = b, r, k, r0, k0, n0)
    (T <= min_T) && (return (; residual = 1e10))   # TODO: Deal with the potential break of feasibility. 
    
    stationary!(ws; R = 1 + r, T = T, w = w0, stationary_kwargs...)
    a = aggregate_assets(h.a_grid, ws.pdf)
    residual = a - k - b
    verbose && next!(p, showvalues = [(:r, r), (:error, residual)])
    return (; residual, a, T)
end

