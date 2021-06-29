#################################################################
# Methods to solve for the General Equilibrium
#################################################################


# Solves for the laissez_faire allocation taking as an input an economy struct.
# r_range defines the range of interest rates to look for an equilibrium. and
# find_zero_args controls the parameters of the root finding. 
function solve_laissez_faire(
    e; 
    tol = (value_function = _TOL_VALUE, distribution = _TOL_PDF),
    r_range = (-0.02, 0.02),
    verbose = true,
    find_zero_args = (method = FalsePosition(3), atol = 1e-10)
)
    t = get_t(e)
    h = get_h(e)

    tmp_arrays_1 = prealloc_and_initialize_ss_PE(
        h; 
        r = r_range[1], 
        w = rL_from_mpk(t, rK_from_r(; t, r = r_range[1])), 
        transfer = 0.0)
    tmp_arrays_2 = prealloc_ss_pdf(h)

    p = ProgressUnknown()
    f = @closure (r) -> begin
        local mpk = mpk_from_after_tax_rK(t, rK_from_r(; t, r))
        local w = rL_from_mpk(t, mpk)
        local sol = solve_stationary_household!(tmp_arrays_1, h, r, w; 
            transfer = 0.0, 
            tol = tol.value_function)
        local pdf = stationary_pdf!(tmp_arrays_2, h, sol; 
            tol = tol.distribution)
        local s = asset_supply(h, pdf)
        local n = labor_supply(h, w)  
        local k = k_from_mpk(t; mpk, n)
        excess = s - k 
        verbose && ProgressMeter.next!(p, 
                        showvalues=[
                            (:r, r),
                            (:excess, excess)
                        ])
        return excess
    end

    # Look for an equilibrium interest rate. 
    r = find_zero(f, r_range, find_zero_args.method, 
        atol = find_zero_args.atol)

    # Construct the rest of the equilibrium given r
    mpk = mpk_from_after_tax_rK(t, rK_from_r(; t, r))
    w = rL_from_mpk(t, mpk)  # Laissez-faire: no taxes, w = rL 
    n = labor_supply(h, w)   # Aggregate labor supply
    k = k_from_mpk(t; mpk, n)
    sol = solve_stationary_household!(tmp_arrays_1, h, r, w; 
        transfer = 0.0, 
        tol = tol.value_function)
    pdf = stationary_pdf!(tmp_arrays_2, h, sol; 
        tol = tol.distribution)
    
    return StationaryEquilibrium(;
        e, r, k, n, w, sol.v, sol.pol, pdf,
        s = asset_supply(h, pdf), 
        b = zero(r), 
        transfer = zero(r),
        y = get_y(t; k, n)
    )
end

# Solves for a new stationary equilibrium given a level of k and a level of b. 
function solve_new_stationary_equilibrium_given_k_b(
    laissez_faire, k::Real, b::Real;
    tol = (value_function = _TOL_VALUE, distribution = _TOL_PDF),
    r_range = (-0.05, 0.0),
    verbose = true,
    find_zero_method = Bisection(), 
    find_zero_args = (atol = 1e-7, xatol = 1e-7)
)
    solve_new_stationary_equilibrium_given_k_b(
        laissez_faire; 
        k_b_fun = (r) -> (k, b), 
        tol, r_range, verbose, find_zero_method, find_zero_args
    )
end 


# Solves for a new stationary equilibrium given a laisses faire economy and a
# k_b_fun describing the level of k and b given an interest rate r
function solve_new_stationary_equilibrium_given_k_b(
    laissez_faire;
    k_b_fun,
    tol = (value_function = _TOL_VALUE, distribution = _TOL_PDF),
    r_range = (-0.05, 0.0),
    verbose = true,
    find_zero_method = Bisection(), 
    find_zero_args = (atol = 1e-7, xatol = 1e-7)
)
    e = get_e(laissez_faire)
    h = get_h(laissez_faire)
    t = get_t(laissez_faire)

    n = laissez_faire.n
    w = laissez_faire.w
    x = laissez_faire.y - (t.δ + laissez_faire.r) * laissez_faire.k

    arrays1 = prealloc_ss_PE(h)
    arrays1.v_0 .= laissez_faire.v
    arrays2 = prealloc_ss_pdf(h)

    min_transfer =  minimum_feasible_transfer(h, w)

    p = ProgressUnknown()
    f! = @closure (r) -> begin
        local (k, b) = k_b_fun(r)  # get the k,b given r 
        local s = k + b
        local c = get_y(t; k, n) - t.δ * k # new stationary aggregate consumption level
        local tr = (c -  x - s * r)

        if tr <= min_transfer
            # transfer is too low 
            distance = - Inf 
        else 
            local sol = solve_stationary_household!(arrays1, h, r, w; 
                transfer = tr, 
                tol = tol.value_function)
            local pdf = stationary_pdf!(arrays2, h, sol; tol = tol.distribution)
            # excess supply 
            distance =  asset_supply(h, pdf) - s 
        end 
        verbose && ProgressMeter.next!(p, 
                        showvalues=[
                            (:tr, tr),
                            (:excess, distance),
                            (:r, r)
                        ])
        return distance 
    end

    r = find_zero(f!, r_range, find_zero_method; find_zero_args...)
    k, b = k_b_fun(r)
    s = k + b
    c = get_y(t; k, n) - t.δ * k
    tr = (c - x - r * s)
    sol = solve_stationary_household!(arrays1, h, r, w; 
        transfer = tr, 
        tol = tol.value_function)
    pdf = stationary_pdf!(arrays2, h, sol; tol = tol.distribution)

    out = StationaryEquilibrium(;
        e, r, b, k, n, w, sol.v, sol.pol, pdf,
        s = asset_supply(h, pdf), 
        transfer = tr,
        y = get_y(t; k, n)
    )

    is_valid(out) || @warn "Equilibrium is not valid"
    return out
end 

