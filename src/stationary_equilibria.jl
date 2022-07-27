

function stationary_laissez_faire(h, t; r_range, 
    verbose = true, 
    hh_problem_kwargs = nothing,
    find_zero_pars = (;atol = _ZERO_FTOL)
)

    hh_problem_options_baseline = (; 
        value_tol = _TOL, 
        policy_tol = _TOL, 
        pdf_tol = _TOL,
        max_iters = _MAX_ITERS,
        verbose = false
    ) 

    stationary_kwargs = isnothing(hh_problem_kwargs) ?
            hh_problem_options_baseline : 
            merge(hh_problem_options_baseline, hh_problem_kwargs)
    
    r0 = r_range[1]
    w = w_from_r(t; r = r0)
    pars = (R = 1 + r0, T = zero(r0), w = w)

    ws = HouseholdWorkspace(; h, pars...)

    f = (r) -> _laissez_faire!(ws, r, h, t, verbose, stationary_kwargs).excess 
    r = find_zero(f, r_range, A42(); find_zero_pars...)
    (; w, n, k, a) = _laissez_faire!(ws, r, h, t, false, stationary_kwargs)

    e = StationaryEquilibrium(;
        h, t, ws, w, n, k, a, r, 
        T = zero(r), b = zero(r)
    )
    return e
end


function _laissez_faire!(ws, r, h, t, verbose, stationary_kwargs)
    w = w_from_r(t; r)
    pars = (R = 1 + r, T = zero(r), w)
    stationary!(ws; stationary_kwargs..., pars...)
    a = asset_supply(h.a_grid, ws.pdf)
    n = labor_supply(h; w = w)
    k = k_from_r_n(t; r, n)
    excess = a - k
    verbose && println("error: $excess r:$r w:$w n:$n k:$k a:$a")
    return (; excess, w, a, n, k)
end


function stationary_equilibrium_given_k_b(e_init, k, b; r_range, 
    verbose = true,   
    hh_problem_kwargs = nothing,
    find_zero_pars = (; atol = _ZERO_FTOL)
)

    hh_problem_options_baseline = (; 
        value_tol = _TOL, 
        policy_tol = _TOL, 
        pdf_tol = _TOL,
        max_iters = _MAX_ITERS,
        verbose = false
    ) 

    stationary_kwargs = isnothing(hh_problem_kwargs) ?
        hh_problem_options_baseline :
        merge(hh_problem_options_baseline, hh_problem_kwargs)

    (;h, t, w) = e_init
    
    k0, n0, r0, w0 = e_init.k, e_init.n, e_init.r, e_init.w
    r = r_range[1]
    T = get_T(t; b, bprime = b, k, r, r0, k0, n0)
    pars = (R = 1 + r, T = T, w = w0)

    ws = HouseholdWorkspace(; h, pars...)
    min_T = minimum_feasible_transfer(h, w0)

    f = (r) -> _k_b_excess!(ws, r, e_init, k, b, min_T, verbose, stationary_kwargs).excess
    r = find_zero(f, r_range, A42(); find_zero_pars...)
    (; a, T) = _k_b_excess!(ws, r, e_init, k, b, min_T, false, stationary_kwargs)

    e = StationaryEquilibrium(; h, t, ws, w = w0, n = n0, k, a, r, b, T)
    return e
end


function _k_b_excess!(ws, r, e_init, k, b, min_T, verbose, stationary_kwargs)
    (; h, t) = e_init
    k0, n0, r0, w0 = e_init.k, e_init.n, e_init.r, e_init.w
    T = get_T(t; b, bprime = b, k, r, r0, k0, n0)
    (T <= min_T) && (return (; excess = 1e10))
    pars = (R = 1 + r, T = T, w = w0)
    
    stationary!(ws; pars..., stationary_kwargs...)
    a = asset_supply(h.a_grid, ws.pdf)
    excess = a - k - b
    verbose && println("error: $excess r:$r T:$T a:$a k:$k b:$b")
    return (; excess, a, T)
end

