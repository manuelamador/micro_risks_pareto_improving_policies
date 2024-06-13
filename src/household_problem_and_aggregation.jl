# This file contains the functions that solve the household bellman equation, aggregate results, as well as compute the stationary solution given R, w, T.
# These functions use only information on the household problem and prices and transfers; and do not use information about technology or fiscal policy. 


######################################
# HOUSEHOLD PROBLEM
######################################


"""
    minimum_feasible_transfer(h, w)

Returns the minimum feasible transfer consistent with feasibility constraints of the household problem given a wage `w`. The struct `h` contains the household parameters.   
"""
minimum_feasible_transfer(h, w) = maximum(disutility_given_w(h.v; w = w * z) - 
    labor_income(h.v; w = w * z) for z in h.z_grid)



"""
    backwards_euler_x(u::CRRA, pp, nn)

Returns a vector of `x_t` that solves the Euler equation with equality with t+1 scaled marginal utilities given by vector `nn` with probability distribution `pp`.  

# Arguments 
- `u`: utility function.
- `pp`: probability vector over states at t+1.
- `nn`: associated vector of "scaled" marginal utility at t+1, i.e. η_t+1 = R_{t+1} x_{t+1}^(-1/ξ).
"""
function backwards_euler_x(u::CRRA, pp, nn)
    ξ = get_inverse_ies(u)
    β = get_β(u)
    Evx = zero(eltype(nn))
    for (p, η) in zip(pp, nn)
        Evx += p * (η)^(-ξ)
    end 
    x = (β * Evx)^(-1/ξ)    
    return x
end


"""
    backwards_euler_x(u::EZ, pp, vv, nn)

Returns a vector of `x_t` that solves the Euler equation with equality with t+1 scaled marginal utilities given by vector `nn` with probability distribution `pp`.

# Arguments
- `u`: utility function.
- `pp`: probability vector over states at t+1.
- `vv`: associated vector of values at t+1.
- `nn`: associated vector of "scaled" marginal utility at t+1, i.e. η_t+1 = R_{t+1} x_{t+1}^(-1/ξ).
"""    
function backwards_euler_x(u::EZ, pp, vv, nn)
    ξ = get_inverse_ies(u)
    γ = get_ra(u)
    β = get_β(u)
    Ev = Evx = zero(eltype(vv))
    for (v, η, p) in zip(vv, nn, pp)
        Evx += p * (v)^((ξ - γ)) * (η)^(-ξ)
        Ev += p * u.risk(v)
    end 
    x = (β * inverse(u.risk, Ev^((γ - ξ))) * Evx)^(-1/ξ)
    return x
end


"""
    _backwards_euler_x_helper(u, P′, i, s, next_values)

Helper function for `backwards_once!` that calls the appropriate `backwards_euler_x` function depending on the type of utility function, and the current state (i, s). 

# Arguments
- `u`: utility function.
- `P′`: transition matrix.
- `i`: index of current asset position.
- `s`: index of current employment state.
- `next_values`: struct containing next period's values, which includes the vector of `η`, and `v` (the latter for the case of Epstein-Zin). 

"""
_backwards_euler_x_helper(u::CRRA, P′, i, s, next_values) = backwards_euler_x(u, view(P′, :, s), view(next_values.η, i, :))
_backwards_euler_x_helper(u::EZ, P′, i, s, next_values) = backwards_euler_x(u, view(P′, :, s), view(next_values.v, i, :), view(next_values.η, i, :))


"""
    backwards_once!(h, current_values; next_values, R, T, w [, a_tmp])

Modifies `current_values` using the Euler equation and the continuation values in `next_values`.

# Arguments 
- `h`: household struct.
- `current_values`: struct containing the matrices `η` and `v` to be computed (the latter for the case of Epstein-Zin).
- `next_values`: struct containing the next period's values, which includes the matrices `η` and `v` (the latter for the case of Epstein-Zin).
- `R`: current gross interest rate.
- `T`: current transfer.
- `w`: current wage.
- `a_tmp`: temporary array.

The values of `current_values.η` and `current_values.v` are modified (the latter only for the case of Epstein-Zin). 
"""
function backwards_once!(h, current_values; next_values, R, T, w, a_tmp = similar(current_values.η))   
    η_now = current_values.η
    (; u) = h
    P′ = h.Pprime
    ξ = get_inverse_ies(u)
    R_ = R^(-1/ξ)

    a_min = first(h.a_grid)
    a_max = last(h.a_grid)

    Threads.@threads for s in eachindex(h.z_grid)
        z = h.z_grid[s]
        add =  labor_income(h.v; w = w * z)  + T - disutility_given_w(h.v; w = w * z)

        # Endogenous Grid Method 
        for i in eachindex(h.a_grid)  
            # find the associated x according the the Euler equation and next_values
            x = _backwards_euler_x_helper(h.u, P′, i, s, next_values)
            a_next = h.a_grid[i]  # associated next period policy
            # computes the level of assets associated with the Euler equation result
            a_tmp[i, s]  = (x + a_next - add) / R
        end

        for i in eachindex(h.a_grid)
            a = h.a_grid[i]
            a_prime = interp1D(a, view(a_tmp, :, s), h.a_grid)
            a_prime = clamp(a_prime, a_min, a_max)
            x = R * a + add - a_prime
            η_now[i, s] = R_ * x
            _backward_once_helper!(u, current_values, next_values, h, a_prime, P′, x, i, s)
        end
    end    
end


_backward_once_helper!(::CRRA, args...) = nothing

"""
    _backward_once_helper!(u::EZ, current_values, next_values, h, a_prime, P′, x, i, s)

Modifies the value function at state (i, s) contained in `current_values` using the values in `next_values` and the optimal policies `x` and `a_prime`. 
This is for the case of Epstein-Zin utility.
"""
function _backward_once_helper!(u::EZ, current_values, next_values, h, a_prime, P′, x, i, s)
    β = u.β

    Ev = zero(a_prime)
    for s2 in axes(next_values.v, 2)
        v_prime = interp1D(a_prime, h.a_grid, view(next_values.v, :, s2))
        Ev += P′[s2, s] * u.risk(v_prime)
    end
    current_values.v[i, s] =
        inverse(u.temporal, (1 - β) * u.temporal(x) + β * u.temporal(inverse(u.risk, Ev)))
end 


"""
    asset_policy_given_η!(a_pol; h, η, R, T, w)

Computes the asset policy given the scaled marginal utility of consumption `η` the current interest rate `R`, transfer `T` and wage `w`. 
Modifies the matrix `a_pol` in place.

# Arguments 
- `a_pol`: matrix to store the policy.
- `h`: household struct.
- `η`: matrix of scaled marginal utility of consumption.
- `R`: current gross interest rate.
- `T`: current transfer.
- `w`: current wage.
"""
function asset_policy_given_η!(a_pol; h, η, R, T, w)
    ξ = get_inverse_ies(h.u) 
    for s in eachindex(h.z_grid), i in eachindex(h.a_grid)
        a_pol[i, s] = R * h.a_grid[i] + T + labor_income(h.v; w = w * h.z_grid[s]) -  η[i, s] * R^(1/ξ) - disutility_given_w(h.v; w = w * h.z_grid[s])
    end 
end

"""
    asset_policy_given_η!(ws::HouseholdWorkspace; R, T, w)

Computes the asset policy over the entire state space given the scaled marginal utility of consumption `η` contained in `ws.η`, and the current interest rate `R`, transfer `T` and wage `w`.
Modifies the matrix `ws.a_pol` in place.
"""
asset_policy_given_η!(ws::HouseholdWorkspace; R, T, w) = asset_policy_given_η!(ws.a_pol; ws.h, ws.η, R, T, w) 



"""
    interpolate_asset_policy!(lower_index, lower_weight; a_grid, a_pol)

Given the policy function `a_pol` and the grid `a_grid`, computes the indices and weights for linear interpolation. 
The matrices `lower_index` and `lower_weight` are modified in place with the lower indices and their associated weights respectively.
"""
function interpolate_asset_policy!(lower_index, lower_weight; a_grid, a_pol)
    Threads.@threads for i in eachindex(a_pol, lower_index, lower_weight)
        a = a_pol[i]
        j = searchsortedfirst(a_grid, a)
        ind = clamp(j, firstindex(a_grid) + 1, lastindex(a_grid))
        l_weight =  (a_grid[ind] - a) / (a_grid[ind] - a_grid[ind - 1]) 
        lower_index[i] = ind - 1
        lower_weight[i] = clamp(l_weight, zero(l_weight), one(l_weight))
    end 
end

"""
    interpolate_asset_policy!(ws::HouseholdWorkspace)

Using the policy function in `ws.a_pol` and the grid `ws.agrid`, computes the indices and weights for linear interpolation.
The matrices `ws.lower_index` and `ws.lower_weight` are modified in place with the lower indices and their associated weights respectively.

"""
interpolate_asset_policy!(ws::HouseholdWorkspace) = interpolate_asset_policy!(ws.lower_index, ws.lower_weight; ws.h.a_grid, ws.a_pol)


"""
    forward_pdf!(pdf_next; h, pdf, lower_index, lower_weight)

Computes the next period probability distribution given the current probability distribution `pdf`, and the "interpolated" information for the asset policy function contained in `lower_index`, `lower_weight`. It updates the matrix `pdf_next` in place.
"""
function forward_pdf!(pdf_next; h, pdf, lower_index, lower_weight)
    P′ = h.Pprime
    fill!(pdf_next, zero(eltype(pdf_next)))
    for s in axes(pdf, 2)
        for i in axes(pdf, 1)
            ind = lower_index[i, s]
            weight = lower_weight[i, s]
            for s′ in axes(pdf, 2)
                mass = P′[s′, s] * pdf[i, s]
                p = mass * weight 
                pdf_next[ind, s′] += p
                pdf_next[ind + 1, s′] += mass - p
            end 
        end 
    end 
end 


######################################
#  AGGREGATE FUNCTIONS
######################################

"""
    aggregate_assets(a_grid, pdf)

Computes the aggregate asset supply.
"""
function aggregate_assets(a_grid, pdf)
    assets = zero(eltype(pdf))
    for s in axes(pdf, 2), i in axes(pdf, 1)
        assets += a_grid[i] * pdf[i, s]
    end 
    return assets
end 


""" 
    aggregate_consumption(h; R, w, T, a_pol, pdf)

Returns the aggregate consumption given the policy function `a_pol` and the probability distribution `pdf`. 

# Arguments
- `h`: household struct.
- `R`: current gross interest rate.
- `w`: current wage.
- `T`: current transfer.
- `a_pol`: policy function.
- `pdf`: probability distribution.
""" 
function aggregate_consumption(h; R, w, T, a_pol, pdf)
    tot = zero(eltype(pdf))
    for s in axes(pdf, 2), i in axes(pdf, 1)
        tot += (R * h.a_grid[i] + T + labor_income(h.v; w = w * h.z_grid[s]) - a_pol[i, s]) * pdf[i, s]
    end 
    return tot 
end 

"""
    aggregate_labor(h; w)

Returns the aggregate labor supply given a wage `w` (given the GHH utility function assumed). 
It assumed that the probability distribution for the employment state is at its stationary distribution contained in `h.Pss`.
"""
function aggregate_labor(h::Household; w)
    (; v, z_grid, Pss) = h
    n = zero(eltype(z_grid))
    for (p, z) in zip(Pss, z_grid)
        n += p * z * labor(v; w = w * z)
    end
    return n
end


"""
    consumption_alloc(h; R, w, T, a_pol)

Returns the consumption allocation given the policy function `a_pol`.

# Arguments
- `h`: household struct.
- `R`: current gross interest rate.
- `w`: current wage.
- `T`: current transfer.
- `a_pol`: policy function matrix.
"""
function consumption_alloc(h::Household; R, w, T, a_pol)
    (; z_grid,  a_grid) = h
    c = similar(a_pol)
    Threads.@threads for s in axes(c, 2)
        for i in axes(c, 1)
            c[i, s] = R * a_grid[i] + T + labor_income(h.v; w = w * z_grid[s]) - a_pol[i, s] 
        end
    end
    return c
end
consumption_alloc(e::StationaryEquilibrium) = consumption_alloc(e.h; R = 1 + e.r, e.w, e.T, e.ws.a_pol) 


is_pol_valid(a_pol, h) = all(a_pol .< last(h.a_grid))
is_pol_valid(e::StationaryEquilibrium) = is_pol_valid(e.ws.a_pol, e.h)


######################################
#  STATIONARY HOUSEHOLD PROBLEM
######################################


"""
    update_stationary_pdf!(ws[; tol, max_iters, print_every, verbose])

Computes the stationary probability distribution given the "interpolated" information for the asset policy function contained in `ws.lower_index`, `ws.lower_weight`. It updates the matrix `ws.pdf` in place with the stationary distribution. It uses the original value of `ws.pdf` as a starting guess, and modifies the matrix `ws.pdf_tmp` in temporary calculations. 
"""
function update_stationary_pdf!(ws; tol = _TOL, max_iters = 10_000, print_every = 100, verbose = false)
    (; h, pdf, pdf_tmp, lower_index, lower_weight) = ws 
    pdf_0, pdf_1 = pdf, pdf_tmp
    iter = 1 
    while true 
        forward_pdf!(pdf_1; h, pdf = pdf_0, lower_index, lower_weight)
        distance = chebyshev(pdf_0, pdf_1)
        verbose && (iter % print_every == 1) && println("iter $iter: $distance")
        if (distance < tol) 
            verbose && println("Converged. Iter $iter: $distance")
            break
        end 
        iter += 1 
        (iter > max_iters) && (@warn "PDF did not converge!"; break) 
        pdf_0, pdf_1 = pdf_1, pdf_0 
    end 
    pdf .= pdf_0
end 


# Helper functions for stationary!
_get_Euler_iterator(::EZ, ws) = ((; v = ws.v, η = ws.η), (; v = ws.v_tmp, η = ws.η_tmp)); 
_get_Euler_iterator(::CRRA, ws) = ((; η = ws.η), (; η = ws.η_tmp)); 
_get_Euler_iterator(ws) = _get_Euler_iterator(ws.h.u, ws)

_distance_Euler_iterators(::CRRA, iterator_0, iterator_1) = chebyshev(iterator_0.η, iterator_1.η)
_distance_Euler_iterators(::EZ, iterator_0, iterator_1) = (; v = chebyshev(iterator_0.v, iterator_1.v), 
          η = chebyshev(iterator_0.η, iterator_1.η)); 

_Euler_iteration_converged(::CRRA, distance, tols) = distance < tols.policy_tol  
_Euler_iteration_converged(::EZ, distance, tols) = distance.η < tols.policy_tol && distance.v < tols.value_tol  

function _stationary_update_Euler_iterator!(::EZ, ws, iterator) 
    ws.v .= iterator.v
    ws.η .= iterator.η
end 
_stationary_update_Euler_iterator!(::CRRA, ws, iterator) = (ws.η .= iterator.η)


"""
    stationary!(ws; R, T, w[, max_iters, print_every, value_tol, policy_tol, pdf_tol, verbose])

Solves the stationary household problem given the utility function `u` and the workspace `ws`.
    
# Arguments
- `ws`: workspace.
- `R`: current gross interest rate.
- `T`: current transfer.
- `w`: current wage.
"""
function stationary!(
    ws; 
    R, T, w, 
    max_iters = _MAX_ITERS, print_every = 100, 
    value_tol = _TOL,
    policy_tol = _TOL,
    pdf_tol = _TOL, 
    verbose = true
)
    (; h, a_tmp) = ws

    iterator_0, iterator_1 = _get_Euler_iterator(ws)
    i = 1
    while true
        backwards_once!(h, iterator_1; next_values = iterator_0, R, T, w, a_tmp)
        distance = _distance_Euler_iterators(h.u, iterator_0, iterator_1)
        if _Euler_iteration_converged(h.u, distance, (; policy_tol, value_tol))
            verbose && println("Converged. iter $i;  R = $R, T = $T, w = $w, distance: $distance")
            break
        end 
        if i > max_iters
            @warn("Max iters reached. DID NOT CONVERGE! R = $R, T = $T, w = $w,  distance: $distance")
            break
        end 
        verbose && (i % print_every == 1) && println("iter $i: $dis1") 
        i += 1
        iterator_0, iterator_1 = iterator_1, iterator_0
    end
    # v, η store the final answer
    _stationary_update_Euler_iterator!(h.u, ws, iterator_0)
    asset_policy_given_η!(ws; R, T, w)
    interpolate_asset_policy!(ws)
    update_stationary_pdf!(ws; tol = pdf_tol)
end 


"""
    stationary(h; R, T, w[, kwargs...])

Returns the stationary solution given the household struct `h` and the prices `R`, `T`, `w`.
"""
function stationary(h; R, T, w, kwargs...)
    ws = HouseholdWorkspace(; h, R, T, w)
    stationary!(ws; R, T, w, kwargs...)
    return ws
end 
