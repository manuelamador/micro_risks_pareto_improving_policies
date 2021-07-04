#################################################################
# Methods to solve for the Stationary General Equilibrium
#################################################################


# Solves for the laissez_faire allocation taking as an input an economy struct.
# r_range defines the range of interest rates to look for an equilibrium. and
# find_zero_args controls the parameters of the root finding.
function solve_laissez_faire(
    e;
    tol = (value_function = _TOL_VALUE, distribution = _TOL_PDF),
    r_range = (-0.02, 0.02),
    verbose = true,
    find_zero_method = Bisection(), # Also try FalsePosition(3) - faster
                                    # but less reliable.
    find_zero_args = (atol = 1e-10, xatol = 1e-8)
)
    t = get_t(e)
    h = get_h(e)

    w0 = rL_from_mpk(t, rK_from_r(; t, r = r_range[1]))
    tmp_arrays_1 = initialize_stationary_household(h; r = r_range[1], w = w0, transfer = 0.0)
    tmp_arrays_2 = prealloc_pdfs(h)

    p = ProgressUnknown()

    # Look for an equilibrium interest rate.
    r = find_zero(r_range, find_zero_method; find_zero_args...) do (r)
        _laissez_faire_excess!((tmp_arrays_1, tmp_arrays_2), r, e, tol, verbose, p)
    end

    # Construct the rest of the equilibrium given r
    (k, w, n) = _laissez_faire_alloc_from_r(e, r)
    transfer = zero(w0)
    sol = solve_stationary_household!(tmp_arrays_1, h, r, w; transfer, tol = tol.value_function)
    pdf = stationary_pdf!(tmp_arrays_2, h, sol; tol = tol.distribution)

    return StationaryEquilibrium(;
        e, r, k, n, w, sol.v, sol.pol, pdf,
        s = asset_supply(h, pdf),
        b = zero(r),
        transfer = zero(r),
        y = output(t; k, n)
    )
end


function _laissez_faire_alloc_from_r(e, r)
    t = get_t(e)
    h = get_h(e)

    rK = rK_from_r(; t, r)
    mpk = mpk_from_after_tax_rK(t, rK)
    w = rL_from_mpk(t, mpk)  # Laissez-faire: no taxes, w = rL
    n = labor_supply(h, w)   # Aggregate labor supply
    k = k_from_mpk(t; mpk, n)
    return (k, w, n)
end


# Helper function to compute the excess asset supply in the laissez_faire calculations
function _laissez_faire_excess!(tmp_arrays, r, e, tol, verbose, p)
    h = get_h(e)
    (tmp_arrays_1, tmp_arrays_2) = tmp_arrays

    (k, w, _) = _laissez_faire_alloc_from_r(e, r)

    transfer = zero(w)
    sol = solve_stationary_household!(tmp_arrays_1, h, r, w; transfer, tol = tol.value_function)
    pdf = stationary_pdf!(tmp_arrays_2, h, sol; tol = tol.distribution)
    s = asset_supply(h, pdf)

    excess = s - k
    verbose && ProgressMeter.next!(p, showvalues=[(:r, r), (:excess, excess)])
    return excess
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
        (r) -> (k, b), laissez_faire;
        tol, r_range, verbose, find_zero_method, find_zero_args)
end


# Solves for a new stationary equilibrium given a laisses faire economy and a
# k_b_fun describing the level of k and b given an interest rate r
function solve_new_stationary_equilibrium_given_k_b(
    k_b_fun,
    laissez_faire;
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

    arrays1 = prealloc_stationary_household(h)
    arrays1.v_0 .= laissez_faire.v
    arrays2 = prealloc_pdfs(h)

    tr_min =  minimum_feasible_transfer(h, w)

    p = ProgressUnknown()
    r = find_zero(r_range, find_zero_method; find_zero_args...) do (r)
        _new_equilibrium_excess!(tmp_arrays, r, k_b_fun, laissez_faire, tr_min, tol, verbose, p)
    end

    k, b = k_b_fun(r)
    c = output(t; k, n) - t.δ * k
    x = laissez_faire.y - (t.δ + laissez_faire.r) * laissez_faire.k
    tr = (c - x - r * (k + b))
    sol = solve_stationary_household!(arrays1, h, r, w; transfer = tr, tol = tol.value_function)
    pdf = stationary_pdf!(arrays2, h, sol; tol = tol.distribution)
    s = asset_supply(h, pdf)
    y = output(t; k, n)

    out = StationaryEquilibrium(;e, r, b, k, n, w, sol.v, sol.pol, pdf, s, transfer = tr, y)
    is_valid(out) || @warn "Equilibrium is not valid"
    return out
end


# Helper function to compute the excess asset supply in the new equilibrium calculations
function _new_equilibrium_excess!(tmp_arrays, r, k_b_fun, laissez_faire, tr_min, tol, verbose, p)
    t = get_t(laissez_faire)
    h = get_h(laissez_faire)

    n = laissez_faire.n
    w = laissez_faire.w
    x = laissez_faire.y - (t.δ + laissez_faire.r) * laissez_faire.k

    tmp_arrays_1, tmp_arrays_2 = tmp_arrays
    (k, b) = k_b_fun(r)  # get k, b given r
    s = k + b
    c = output(t; k, n) - t.δ * k # new stationary aggregate consumption level
    tr = (c -  x - s * r)

    if tr <= tr_min
        # transfer is too low
        distance = - Inf
    else
        sol = solve_stationary_household!(tmp_arrays_1, h, r, w; transfer = tr, tol = tol.value_function)
        pdf = stationary_pdf!(tmp_arrays_2, h, sol; tol = tol.distribution)
        # excess supply
        distance =  asset_supply(h, pdf) - s
    end

    verbose && ProgressMeter.next!(p, showvalues=[(:tr, tr), (:excess, distance), (:r, r)])
    return distance
end
