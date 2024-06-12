#######################################################
#              AGGREGATE RISK PATHS                   #
#######################################################

"""
   _transition_item_AR(h::Household)
   
Creates a workspace for the transition path with aggregate shocks at one date.
"""
function _transition_item_AR(h::Household)  
    (; v, η, a_pol, pdf, lower_index, lower_weight, a_tmp) = _generate_base_workspace_matrices(h)

    r = T = a = k = b = w = n = zero(eltype(h.a_grid)) 
    return (; h, v, η, a_pol, pdf, lower_index, lower_weight, a_tmp, r, T, a, k, b, w, n)
end 


"""
    _initialize_path_AR(e_init, e_final; b_path)

Initializes the transition path from `e_init` to `e_final` with the given `b_path`. The length of the transition path 
is the length of `b_path`. Here, there is a path for each possible realization of the aggregate productivity shock
"""
function _initialize_path_AR(e_init, e_final; b_path)
    h = e_init.h


    lst = StructArray(_transition_item_AR(h) for _ in b_path)
    for i = 1:size(b_path,2)
        lst[1,i].pdf .= e_init.ws.pdf

        lst[end,i].v .= e_final.ws.v 
        lst[end,i].η .= e_final.ws.η
        lst[end,i].a_pol .= e_final.ws.a_pol
        lst[end,i].lower_index .= e_final.ws.lower_index
        lst[end,i].lower_weight .= e_final.ws.lower_weight

        lst.r[end,i] = e_final.r
        lst.b[end,i] = e_final.b
        lst.k[end,i] = e_final.k
        lst.a[end,i] = e_final.a
        lst.w[end,i] = e_final.w
        lst.n[end,i] = e_final.n
    end

    return lst 
end

"""
    solve_laissez_faire_transition_AR(e_init, e_final; b_path, A_path, probH, r_path, verbose, nlsolve_kwargs)

Solves for the RPI transition (constant wage/profits) path from a stationary equilibrium `e_init` to a sationary equilibrium `e_final` with the given `k_path`, `b_path`. The function uses a robust but potentially slow algorithm to solve for the transition path.

# Arguments
- `e_init`: Initial equilibrium.
- `e_final`: Final equilibrium.
- `b_path`: Paths of government debt levels, one for each possible realization of the aggregate productivity shock. This is a matrix of zeros in the laissez-faire case.
- `A_path`: Paths of aggregate productivity, one for each possible realization of the aggregate productivity shock.
- `probH`: Probability of the high aggregate productivity shock
- `r_path`: Path of interest rates. If `nothing`, the function assumes that the interest rate is constant at the final equilibrium value.
- `verbose`: Whether to show progress.
- `nlsolve_kwargs`: Keyword arguments to pass to the `nlsolve` solver.
"""

function solve_laissez_faire_transition_AR(
    e_init, e_final; b_path, A_path, probH, r_path = nothing, verbose = true,
    nlsolve_kwargs = nothing)
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

    path = _initialize_path_AR(e_init, e_final; b_path)

    if isnothing(r_path) 
        r_path_rest = fill(e_final.r, (size(b_path,1)-1,size(b_path,2)))
    else 
        r_path_rest = r_path[2:end,:]
    end 
    r_path_rest[end,:] .= e_final.r

    f! = (F, x) -> _laissez_faire_transition_residuals_AR(F, x; path, e_init, e_final, b_path, A_path, probH)
    r_sol = nlsolve(f!, r_path_rest; nlsolve_kwargs_merged...)
    _laissez_faire_transition_residuals_AR(similar(r_sol.zero), r_sol.zero; path, e_init, e_final, b_path, A_path, probH)
   
    return path 
end

"""
    _laissez_faire_transition_residuals_AR(residuals, r_path_rest; path, e_init, e_final, b_path, Apth, probH)

Calculates the residuals of the transition path given the interest rate path `r_path_rest`. 
It modifies `residuals` in place with the excess asset supply at each date.
It also modifies the `path` workspace in place when doing the calculations.
"""

function _laissez_faire_transition_residuals_AR(residuals, r_path_rest; path, e_init, e_final, b_path, A_path, probH)
    (; h, t) = e_init

    r_min = - t.δ + eps()  # minimum interest rate
    r_max = 1 / h.u.β  # a maximum -- could be relaxed
    
    for j in eachindex(b_path[1,:])
        for i in reverse(eachindex(path[:,1]))
            if i > firstindex(path[:,j])
                r = r_path_rest[i - 1,j]
                r = clamp(r, r_min, r_max)
                R = 1 + r
                T = 0.0            
                
                path.T[i,j] = T
                path.r[i,j] = r
                next_values = i == lastindex(path[:,j]) ? e_final.ws : path[i+1,j]
                t.A = A_path[i,j]
                w = mpl_from_mpk(t; mpk = r + t.δ)
                path.w[i,j] = w
                backwards_once!(h, path[i,j]; next_values = next_values, R, T, w, path[i,j].a_tmp) # Here, when going backwards, it needs to take into account that A changes
                asset_policy_given_η!(path[i,j].a_pol; h, path[i,j].η, R, T, w)
                interpolate_asset_policy!(path[i,j].lower_index, path[i,j].lower_weight; h.a_grid, path[i,j].a_pol)
                n = aggregate_labor(h; w = w)
                path.n[i,j] = n
                path.k[i,j] = k_from_mpk_n(t; mpk = r + t.δ, n)
            end
        end
    end
    # Period 1
    r = e_init.r
    R = 1 + r
    path.r[1,:] .= r
    T = 0.0
    for j in eachindex(b_path[1,:])
        t.A = A_path[1,j]
        w = mpl_from_mpk(t; mpk = r + t.δ)
        path.T[1,j] = T
        path.w[1,j] = w
        backwards_once_t0!(h, path[1,j]; next_values_H = path[2,1], next_values_L = path[2,2], R, T, w, path[1,j].a_tmp, probH)
        asset_policy_given_η!(path[1,j].a_pol; h, path[1,j].η, R, T, w)
        interpolate_asset_policy!(path[1,j].lower_index, path[1,j].lower_weight; h.a_grid, path[1,j].a_pol)
        n = aggregate_labor(h; w = w)
        path.n[1,j] = n
        path.k[1,j] = k_from_mpk_n(t; mpk = r + t.δ, n)
    end

    for j in eachindex(b_path[1,:])       
        for i in eachindex(path[:,j]) 
            a = aggregate_assets(h.a_grid, path[i,j].pdf)
            path.a[i,j] = a
            path.b[i,j] = b_path[i,j]
            if i > firstindex(path[:,j])
                residuals[i - 1, j] =  a - path.k[i,j] - path.b[i,j]
            end
            (i == lastindex(path[:,j])) && break
            forward_pdf!(path[i+1,j].pdf; h, path[i,j].pdf, path[i,j].lower_index, path[i,j].lower_weight)
        end
    end
end

"""
    solve_RPI_transition_AR(e_init, e_final; k_path, b_path, r_path_lf, w_path_lf, probH, r_path, verbose, nlsolve_kwargs)

Solves for the RPI transition (laissez_faire economy with agg. shocks) path from `e_init` to `e_final` with the given `k_path`, `b_path` for each possible realization of the aggregate productivity shock.

# Arguments
- `e_init`: Initial equilibrium.
- `e_final`: Final equilibrium.
- `k_path`: Path of capital levels, one for each possible realization of the aggregate productivity shock.
- `b_path`: Path of government debt levels, one for possible realization of the aggregate productivity shock.
- `r_path_lf`: Path of interest rates in the laissez_faire economy. 
- `w_path_lf`: Path of wages in the laissez_faire economy.
- `probH`: Probability of the high aggregate productivity shock
- `r_path`: Path of interest rates. If `nothing`, the function assumes that the interest rate is constant at the final equilibrium value.
- `verbose`: Whether to show progress.
- `nlsolve_kwargs`: Keyword arguments to pass to the `nlsolve` solver.
"""

function solve_RPI_transition_AR(e_init, e_final; k_path, b_path, r_path_lf, w_path_lf, probH, r_path = nothing, verbose = true,
    nlsolve_kwargs = nothing)
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

    path = _initialize_path_AR(e_init, e_final; b_path)

    if isnothing(r_path) 
        r_path_rest = fill(e_final.r, (size(k_path,1)-1,size(k_path,2)))
    else 
        r_path_rest = r_path[2:end,:]
    end 
    r_path_rest[end,:] .= e_final.r

    f! = (F, x) -> _RPI_transition_AR_residuals(F, x; path, e_init, e_final, k_path, b_path, r_path_lf, w_path_lf, probH)
    r_sol = nlsolve(f!, r_path_rest; nlsolve_kwargs_merged...)
    _RPI_transition_AR_residuals(similar(r_sol.zero), r_sol.zero; path, e_init, e_final, k_path, b_path, r_path_lf, w_path_lf, probH)
    
    return path 
end

"""
    _RPI_transition_AR_residuals(residuals, r_path_rest; path, e_init, e_final, k_path, b_path, r_path_lf, w_path_lf, probH)

Calculates the residuals of the transition path given the interest rate path `r_path_rest`. 
It modifies `residuals` in place with the excess asset supply at each date and for each possible realization of the aggregate productivity shock.
It also modifies the `path` workspace in place when doing the calculations.
"""

function _RPI_transition_AR_residuals(residuals, r_path_rest; path, e_init, e_final, k_path, b_path, r_path_lf, w_path_lf, probH)
    (; h, t) = e_init
    
    r_min = - t.δ + eps()  # minimum interest rate
    r_max = 1 / h.u.β  # a maximum -- could be relaxed
    
    for j in eachindex(k_path[1,:])
        for i in reverse(eachindex(k_path[:,1],b_path[:,1],path[:,1]))
            if i > firstindex(path[:,j])
                b, k = b_path[i,j], k_path[i,j]
                bprime = i == lastindex(b_path[:,j]) ? e_final.b : b_path[i + 1, j] 
                r = r_path_rest[i - 1, j]
                r = clamp(r, r_min, r_max)
                R = 1 + r

                r_lf = r_path_lf[i,j] # Note here that the index is i, instead of i-1, becuase r_path_lf has size 150x1
                T = get_T_AR(; b, bprime, r = r, r_lf, k)

                w = w_path_lf[i,j] # Get w from the guess of r
                min_T = nextfloat(minimum_feasible_transfer(h, w)) # minimum transfer
                T = max(T, min_T)
            
                path.T[i,j] = T
                path.r[i,j] = r
                path.w[i,j] = w
                next_values = i == lastindex(path[:,j]) ? e_final.ws : path[i+1,j]
                backwards_once!(h, path[i,j]; next_values = next_values, R, T, w, path[i,j].a_tmp)
                asset_policy_given_η!(path[i,j].a_pol; h, path[i,j].η, R, T, w)
                interpolate_asset_policy!(path[i,j].lower_index, path[i,j].lower_weight; h.a_grid, path[i,j].a_pol)
            end
        end
    end
    # Period 1
    r = e_init.r
    R = 1 + r
    path.r[1,:] .= r
    for j in eachindex(k_path[1,:])
        b, k = b_path[1,j], k_path[1,j]
        bprime =  b_path[2,j]
        w = w_path_lf[1,j]
        min_T = nextfloat(minimum_feasible_transfer(h, w)) # minimum transfer

        r_lf = r_path_lf[1,j]
        T = get_T_AR(; b, bprime, r = r, r_lf, k) 
        T = max(T, min_T)
        path.T[1,j] = T
        path.w[1,j] = w
        backwards_once_t0!(h, path[1,j]; next_values_H = path[2,1], next_values_L = path[2,2], R, T, w, path[1,j].a_tmp, probH)
        asset_policy_given_η!(path[1,j].a_pol; h, path[1,j].η, R, T, w)
        interpolate_asset_policy!(path[1,j].lower_index, path[1,j].lower_weight; h.a_grid, path[1,j].a_pol)
    end

    for j in eachindex(k_path[1,:])    
        for i in eachindex(path[:,j]) 
            a = aggregate_assets(h.a_grid, path[i,j].pdf)
            path.a[i,j] = a
            path.b[i,j] = b_path[i,j]
            path.k[i,j] = k_path[i,j]
            if i > firstindex(path[:,j])
                residuals[i - 1, j] =  a - k_path[i,j] - b_path[i,j]
            end
            (i == lastindex(path[:,j])) && break
            forward_pdf!(path[i+1,j].pdf; h, path[i,j].pdf, path[i,j].lower_index, path[i,j].lower_weight)
        end
    end
end 

################################################
# Portfolio problem
###############################################

"""
    solve_RPI_transition_portfolio_AR(e_init, e_final; k_path, b_path, r_path_lf, w_path_lf, probH, r_path, verbose, nlsolve_kwargs)

Solves for the RPI transition (laissez_faire economy with agg. shocks) path from `e_init` to `e_final` with the given `k_path`, `b_path` for each possible realization of the aggregate productivity shock.
Households now solve a portfolio problem, choosing how many units of bonds and capital to hold.

# Arguments
- `e_init`: Initial equilibrium.
- `e_final`: Final equilibrium.
- `k_path`: Path of capital levels, one for each possible realization of the aggregate productivity shock.
- `b_path`: Path of government debt levels, one for possible realization of the aggregate productivity shock.
- `r_path_lf`: Path of interest rates in the laissez_faire economy. 
- `w_path_lf`: Path of wages in the laissez_faire economy.
- `probH`: Probability of the high aggregate productivity shock
- `r_path`: Path of interest rates. If `nothing`, the function assumes that the interest rate is constant at the final equilibrium value.
- `verbose`: Whether to show progress.
- `nlsolve_kwargs`: Keyword arguments to pass to the `nlsolve` solver.
"""
function solve_RPI_transition_portfolio_AR(
    e_init, e_final; k_path, b_path, r_path_lf, w_path_lf, probH, r_path = nothing, verbose = true,
    nlsolve_kwargs = nothing)
    nlsolve_baseline_kwargs = (
        ftol = 1e-06,
        show_trace = verbose,
        method = :anderson, 
        m = 10, 
        iterations = 500, #250
        beta = -0.001 #-0.01  # adjustment parameter for the fixed point algorithm
    )
    
    nlsolve_kwargs_merged = isnothing(nlsolve_kwargs) ?             
        nlsolve_baseline_kwargs : nlsolve_kwargs

    path = _initialize_path_AR(e_init, e_final; b_path)

    port = _generate_base_workspace_matrices_T12(e_init.h)

    if isnothing(r_path) 
        r_path_rest = fill(e_final.r, (size(b_path,1),size(b_path,2)))
    else 
        r_path_rest = r_path
    end 
    r_path_rest[end,:] .= e_final.r

    f! = (F, x) -> _RPI_transition_portfolio_AR_residuals(F, x; path, port, e_init, e_final, k_path, b_path, r_path_lf, w_path_lf, probH)
    r_sol = nlsolve(f!, r_path_rest; nlsolve_kwargs_merged...)
    _RPI_transition_portfolio_AR_residuals(similar(r_sol.zero), r_sol.zero; path, port, e_init, e_final, k_path, b_path, r_path_lf, w_path_lf, probH)

    return path, port
end

"""
    _RPI_transition_portfolio_AR_residuals(residuals, r_path_rest; path, e_init, e_final, k_path, b_path, r_path_lf, w_path_lf, probH)

Calculates the residuals of the transition path given the interest rate path `r_path_rest`. 
It modifies `residuals` in place with the excess asset supply at each date and for each possible realization of the aggregate productivity shock.
It also modifies the `path` workspace in place when doing the calculations.
"""

function _RPI_transition_portfolio_AR_residuals(residuals, r_path_rest; path, port, e_init, e_final, k_path, b_path, r_path_lf, w_path_lf, probH)
    (; h, t) = e_init
    (; u) = h
    r_min = - t.δ + eps()  # minimum interest rate
    r_max = 1 / h.u.β  # a maximum -- could be relaxed
    P′ = h.Pprime
    ξ = get_inverse_ies(u)
    γ = get_ra(u)
    a_min = first(h.a_grid)
    a_max = last(h.a_grid)
    β = u.β

    # Transition from period T up to period 3, inclusive
    for j in eachindex(k_path[1,:])
        for i in reverse(eachindex(k_path[:,1],b_path[:,1],path[:,1]))
            if i > firstindex(path[:,j])
                b, k = b_path[i,j], k_path[i,j]
                bprime = i == lastindex(b_path[:,j]) ? e_final.b : b_path[i + 1, j] 
                r = r_path_rest[i, j] # Note here that r_path_rest now has dimension 150x2. 
                r = clamp(r, r_min, r_max)
                R = 1 + r

                r_lf = r_path_lf[i,j] # Note here that the index is i, instead of i-1, becuase r_path_lf has size 150x1
                w = w_path_lf[i,j] # Get w from the guess of r
                T = get_T_AR(; b, bprime, r = r, r_lf, k)
                min_T = nextfloat(minimum_feasible_transfer(h, w)) # minimum transfer
                T = max(T, min_T)
            
                path.T[i,j] = T
                path.r[i,j] = r
                path.w[i,j] = w
                next_values = i == lastindex(path[:,j]) ? e_final.ws : path[i+1,j]
                backwards_once!(h, path[i,j]; next_values = next_values, R, T, w, path[i,j].a_tmp) # Here, when going backwards, it needs to take into account that A changes
                asset_policy_given_η!(path[i,j].a_pol; h, path[i,j].η, R, T, w)
                interpolate_asset_policy!(path[i,j].lower_index, path[i,j].lower_weight; h.a_grid, path[i,j].a_pol)
                (i == 3) && break
            end
        end
    end
    ## Period 2
    r_b = r_path_rest[1, 2] # Note here that we are storing the return on bonds in period 2 (i.e, the "risk-free rate") in r_path_rest[1, 2], even though "1" refers to the first period
    r_k_h = r_path_rest[2, 1]
    r_k_l = r_path_rest[2, 2]
    R_b = 1 + r_b
    R_k_h = 1 + r_k_h
    R_k_l = 1 + r_k_l
    path.r[2,1] = r_k_h
    path.r[2,2] = r_k_l
    for j in eachindex(k_path[1,:])
        b, k = b_path[2,j], k_path[2,j]
        bprime = b_path[3, j] 

        r_lf = r_path_lf[2,j]
        T = j == 1 ? get_T_portfolio_AR(; b, bprime, r_b = r_b, r_k = r_k_h, r_lf, k) : get_T_portfolio_AR(; b, bprime, r_b = r_b, r_k = r_k_l, r_lf, k)

        w = w_path_lf[2,j] 
        min_T = nextfloat(minimum_feasible_transfer(h, w)) # minimum transfer
        T = max(T, min_T)
        path.T[2,j] = T
        path.w[2,j] = w
        next_values = path[3,j]

        R_a_tmp = similar(path[2,j].η)
        x_tmp = similar(path[2,j].η)
        v_tmp = similar(path[2,j].η)
        u_c_tmp = similar(path[2,j].η)
        
        for s in eachindex(h.z_grid)
            z = h.z_grid[s]
            add =  labor_income(h.v; w = w * z)  + T - disutility_given_w(h.v; w = w * z)

            for i in eachindex(h.a_grid)
                x = _backwards_euler_x_helper(h.u, P′, i, s, next_values)
                x_tmp[i,s] = x
                a_next = h.a_grid[i]
                R_a_tmp[i, s]  = (x + a_next - add) # R[i,s]*a_tmp[i,s] (this is cash-in-hand in period 2: cih(a_3,z_2))

                # Solve for the value for each z and a_next [i,s]
                Ev = zero(a_next)
                for s2 in axes(next_values.v, 2)
                    v_prime = next_values.v[i,s2] #v_prime = interp1D(a_prime, h.a_grid, view(next_values.v, :, s2))
                    Ev += P′[s2, s] * u.risk(v_prime)
                end
                v_tmp[i, s] =
                    inverse(u.temporal, (1 - β) * u.temporal(x) + β * u.temporal(inverse(u.risk, Ev))) #v(a_3, z_2)

                u_c_tmp[i, s] = (v_tmp[i, s])^((ξ - γ)) * (x_tmp[i,s])^(-ξ) # u'(c_2) where [i,s] refers to (a_3, z_2). 
            end
            x_itp = LinearInterpolation((R_a_tmp[:,s]),x_tmp[:,s], extrapolation_bc=Line()) # x(cih_2,z_2)
            v_itp = LinearInterpolation((R_a_tmp[:,s]),v_tmp[:,s], extrapolation_bc=Line()) # v(cih_2,z_2)
            u_c_itp = LinearInterpolation((R_a_tmp[:,s]),u_c_tmp[:,s], extrapolation_bc=Line()) # u_c(cih_2,z_2)
            a_itp = LinearInterpolation((R_a_tmp[:,s]),h.a_grid, extrapolation_bc=Line()) # Interpolate also a_3 for each cash-in-hand and z_2: a_3(cih_2,z_2)
            port.lst_x[s,j] = x_itp # s here indexes z_2, j indexes A_2. Each [s,j] contains the interpolation x(cih_2,z_2)
            port.lst_v[s,j] = v_itp # s here indexes z_2, j indexes A_2. Each [s,j] contains the interpolation v(cih_2,z_2)
            port.lst_u_c[s,j] = u_c_itp # Store the interpolation for s=z_2 and j={High, Low}. Here lst_u_c[s,j](cih) would give u_c(cih_2,z_2,A). 
            port.lst_a[s,j] = a_itp
        end
        port.lst_x_tmp[1,j] = x_tmp #x(z_2,a_3) for high and low productivity
        port.lst_v_tmp[1,j] = v_tmp #v(z_2,a_3) for high and low productivity 
        port.lst_R_a_tmp[1,j] = R_a_tmp #cih(z_2,a_3) for high and low productivity
    end
    # Period 1
    r = r_path_rest[1, 1]
    r = clamp(r, r_min, r_max)
    R = 1 + r # This is R_0, given by the stationary steady state
    path.r[1,1] = r
    path.r[1,2] = r_b
    theta_tmp!(h, port; R_b=R_b, R_k_h=R_k_h, R_k_l=R_k_l, probH=probH) # This updates port. In particular, port.theta_star(a_2,z_1). See household_problem_and_aggregation.jl

    for j in eachindex(k_path[1,:])
        b, k = b_path[1,j], k_path[1,j]
        bprime = b_path[2, j] 
        #r = r_path_rest[i - 1, j]
        
        r_lf = r_path_lf[1,j] # Note here that the index is i, instead of i-1, because r_path_lf has size 150x1
        w = w_path_lf[1,j] # Get w from the guess of r
        T = get_T_AR(; b, bprime, r = r, r_lf, k)
        min_T = nextfloat(minimum_feasible_transfer(h, w)) # minimum transfer
        T = max(T, min_T)

        path.T[1,j] = T
        path.w[1,j] = w

        backwards_once_t0_port!(h, path[1,j]; port, R, R_b, R_k_h, R_k_l, T, w, path[1,j].a_tmp, probH)
        interpolate_asset_policy!(path[1,j].lower_index, path[1,j].lower_weight; h.a_grid, path[1,j].a_pol) 
    end
    # These lines could go in period 2
    # Before starting the forward iteration, build the distribution of households in period 2 (the state is cash-in-hand, although we use the same a_grid for cih_grid)
    cih_H = ((1.0.-port.theta_pol).*R_b .+ port.theta_pol .* R_k_h).*path[1,1].a_pol # cih_H[z_1,a_1] in period 2
    cih_L = ((1.0.-port.theta_pol).*R_b .+ port.theta_pol .* R_k_l).*path[1,1].a_pol # cih_L[z_1,a_1] in period 2
    #cih_grid = (1+r_path_lf[1,1]).* h.a_grid # Use the same grid for cash-in-hand
    cih_grid = h.a_grid # Use the same grid for cash-in-hand
    # Interpolation: needed to build a distribution of households in the space of cash-in-hand_2 and z_2
    interpolate_asset_policy!(port.lower_index_H, port.lower_weight_H; a_grid = cih_grid, a_pol = cih_H) # Implicit here the cih_grid is the same than the a_grid
    interpolate_asset_policy!(port.lower_index_L, port.lower_weight_L; a_grid = cih_grid, a_pol = cih_L) # Implicit here the cih_grid is the same than the a_grid
    # We know from previous two lines where cih_H and cih_L live in cih_grid (where cih_H[a_1,z_1] and cih_L[a_1,z_1]).

    # Let's use the interpolations to get a_3[cih_2,z_2] for A high and low. Store them in port.a_pol_H and port.a_pol_L
    # Build a_pol_2. This is the assets chosen in period 2 to carry into period 3: a_3(cih_2,z_2) for each A 
    for i in axes(port.a_pol_H,1) # for each cash-in-hand in the grid cih_2
        for s in axes(port.a_pol_H,2) # for each z_2
            port.a_pol_H[i,s] = clamp(port.lst_a[s,1](cih_grid[i]),a_min, a_max) # lst_a[s,1]: the "1" refers to A_H. This lst_a are interpolations we build in line 145. Clamp to get a_3[cih_2,z_2] when A_2 is high
            port.a_pol_L[i,s] = clamp(port.lst_a[s,2](cih_grid[i]),a_min, a_max) # lst_a[s,2]: the "2" refers to A_L.  Clamp to get a_3[cih_2,z_2] when A_2 is low
        end
    end
    # Later, use interpolate_asset_policy with the a_grid and port.a_pol_H and port.a_pol_L to get path[,:].lower_weight and path[,:].lower_index
    # The following two lines could be used to compute asset supply in period 2
    interpolate_asset_policy!(path[2,1].lower_index, path[2,1].lower_weight; h.a_grid, a_pol = port.a_pol_H) # I know here distribution of assets in period 3 (in a_grid)
    interpolate_asset_policy!(path[2,2].lower_index, path[2,2].lower_weight; h.a_grid, a_pol = port.a_pol_L) # I know here distribution of assets in period 3 (in a_grid)

    ## Forward iteration from period 1 to period T
    k_supply = sum(port.k_pol .* path[1,1].pdf) # Supply of capital in period 1
    b_supply = sum(port.b_pol .* path[1,1].pdf) # Supply of bonds in period 1

    for j in eachindex(k_path[1,:])    
        for i in eachindex(path[:,j])        
            if i==3 && j==1 # CHECK IF a = asset_supply(h.a_grid, path[i,j].pdf) GET THE SAME a_3, AS BELOW I GET PDF IN SPACE [a_3,z_3]
                a = sum(port.a_pol_H .* port.pdf_H) #Sum of policies a_3[cih_2,z_2] across pdf in [cih_2,z_2] when A is High
            elseif i==3 && j==2
                a = sum(port.a_pol_L .* port.pdf_L) #Sum of policies a_3[cih_2,z_2] across pdf in [cih_2,z_2] when A is Low
            else
                a = aggregate_assets(h.a_grid, path[i,j].pdf) # In i=1, a is the same if j=1 or j=2. KEEP THIS LINE IF i=1 or i=2. Note it works also in i=2
            end
            path.a[i,j] = a
            path.b[i,j] = b_path[i,j]
            path.k[i,j] = k_path[i,j]

            if i==1 && j==2 # Residual for the bonds market clearing in period 2. NOTE! we are storing it in residual[1,2], even though is a. eq. condition in the second period 
                residuals[1, 2] =  b_supply - b_path[2,1] # Note b_path[2,1] must be equal to b_path[2,2], i.e, the debt policy in period 2 must be the same in High and Low 
            else
                residuals[i, j] =  a - k_path[i,j] - b_path[i,j] # Here residuals[1,1] is the resid in the stationary distribution. Should be cleared with e_lf.r
            end
            (i == lastindex(path[:,j])) && break
            if i==1
                # Update pdf in the space [a_2,z_2]. Needed to get asset supply in period 2
                forward_pdf!(path[2,j].pdf; h, path[1,j].pdf, path[1,j].lower_index, path[1,j].lower_weight) # Updates path[2,j].pdf , i.e, the distribution of households in [a_2,z_2]. Should be the same in High and Low 
                # Update pdf in the space [cih_2,z_2] for High and Low. Needed to get asset supply in period 3, in High and Low
                forward_pdf!(port.pdf_H; h, path[1,1].pdf, lower_index = port.lower_index_H, lower_weight = port.lower_weight_H) # Here we get pdf in [cih_2,z_2] when A is High starting from a pdf in the space [z_1,a_1]
                forward_pdf!(port.pdf_L; h, path[1,2].pdf, lower_index = port.lower_index_L, lower_weight = port.lower_weight_L) # Note path[1,1].pdf = path[1,2].pdf (given by the stationary distribution)
            elseif i==2
                #Update pdf in space [a_3,z_3]. Needed to get asset supply in period 4 (and could be used to get asset supply in t=3 using a = asset_supply(h.a_grid, path[i,j].pdf))
                forward_pdf!(path[3,1].pdf; h, pdf = port.pdf_H, path[2,1].lower_index, path[2,1].lower_weight) # This updates the pdf in [a_3,z_3] starting from a distribution in period 2 in [cih_2,z_2] and using the policy a_3[cih_2,z_2] to do this
                forward_pdf!(path[3,2].pdf; h, pdf = port.pdf_L, path[2,2].lower_index, path[2,2].lower_weight) # This updates the pdf in [a_3,z_3] starting from a distribution in period 2 in [cih_2,z_2] and using the policy a_3[cih_2,z_2] to do this
            else # i>=3 # Note in period 2 we do not update path[3,j].pdf, as it will not be used to compute assets in period 3
                forward_pdf!(path[i+1,j].pdf; h, path[i,j].pdf, path[i,j].lower_index, path[i,j].lower_weight)
            end
        end
    end
end