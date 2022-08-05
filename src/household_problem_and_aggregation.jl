
# Household problem and aggregation 

function minimum_feasible_transfer(h, w)
    (; z_grid) = h
    wi = w .* z_grid
    y = map(x -> labor_income(h.v; w = x), wi)
    vofn = map(x -> disutility_given_w(h.v; w = x), wi)
    return -minimum(y .- vofn)
end


# Household's euler equation 


function backwards_euler_x(u::CRRA, π_vec, η_vec)
    ξ = get_inverse_ies(u)
    β = get_β(u)
    Evx = zero(eltype(η_vec))
    for (π, η) in zip(π_vec, η_vec)
        Evx += π * (η)^(-ξ)
    end 
    x = (β * Evx)^(-1 / ξ)
    return x
end


function backwards_euler_x(u::EZ, π_vec, v_vec, η_vec)
    ξ = get_inverse_ies(u)
    γ = get_ra(u)
    β = get_β(u)
    Ev = Evx = zero(eltype(v_vec))
    for (v, η, π) in zip(v_vec, η_vec, π_vec)
        Evx += π * (v)^((ξ - γ)) * (η)^(-ξ)
        Ev += π * u.risk(v)
    end 
    x = (β * inverse(u.risk, Ev^((γ - ξ))) * Evx)^(-1 / ξ)
    return x
end


function backwards_once!(h, current_values; next_values, R, T, w, a_tmp = similar(current_values.η))   
    η_now = current_values.η

    (; u, P) = h
    ξ = get_inverse_ies(u)

    @batch for s in eachindex(h.z_grid)
        z = h.z_grid[s]
        dis = disutility_given_w(h.v; w = w * z)
        lab = labor_income(h.v; w = w * z)

        @inbounds for i in eachindex(h.a_grid)
            x = _backwards_euler_x_helper(h.u, P, i, s, next_values)
            c = x + dis
            a_next = h.a_grid[i]
            a_tmp[i, s]  = (c + a_next - lab - T) / R
        end

        @inbounds for i in eachindex(h.a_grid)
            a = h.a_grid[i]
            a_prime = interp1D(a, view(a_tmp, :, s), h.a_grid)
            a_prime = clamp(a_prime, first(h.a_grid), last(h.a_grid))
            max_feasible_a_prime = R * a + T + lab - dis
            x = max_feasible_a_prime - a_prime
            η_now[i, s] = R^(-1/ξ) * x
            _extra_updates!(u, current_values, next_values, h, a_prime, P, x, i, s)
        end
    end    
end

_backwards_euler_x_helper(u::CRRA, P, i, s, next_values) = backwards_euler_x(u, view(P, s, :), view(next_values.η, i, :))
_backwards_euler_x_helper(u::EZ, P, i, s, next_values) = backwards_euler_x(u, view(P, s, :), view(next_values.v, i, :), view(next_values.η, i, :))

_extra_updates!(::CRRA, args...) = nothing
function _extra_updates!(u::EZ, current_values, next_values, h, a_prime, P, x, i, s)
    β = u.β

    Ev = zero(a_prime)
    for s2 in axes(next_values.v, 2)
        v_prime = interp1D(a_prime, h.a_grid, view(next_values.v, :, s2))
        Ev += P[s, s2] * u.risk(v_prime)
    end
    current_values.v[i, s] =
        inverse(u.temporal, (1 - β) * u.temporal(x) + β * u.temporal(inverse(u.risk, Ev)))
end 


function asset_policy_given_η!(a_pol; h, η, R, T, w)
    ξ = get_inverse_ies(h.u) 
    @tullio a_pol[i, s] = R * h.a_grid[i] + T + labor_income(h.v; w = w * h.z_grid[s]) -  η[i, s] * R^(1/ξ) - disutility_given_w(h.v; w = w * h.z_grid[s])
end

asset_policy_given_η!(ws::HouseholdWorkspace; R, T, w) = asset_policy_given_η!(ws.a_pol; ws.h, ws.η, R, T, w) 


function interpolate_asset_policy!(lower_index, lower_weight; a_grid, a_pol)
    @batch for s in axes(a_pol, 2)
        @inbounds for i in axes(a_pol, 1)
            a = a_pol[i, s]
            j = searchsortedfirst(a_grid, a)
            ind = clamp(j, firstindex(a_grid) + 1, lastindex(a_grid))
            l_weight =  (a_grid[ind] - a) / (a_grid[ind] - a_grid[ind - 1]) 
            lower_index[i, s] = ind - 1
            lower_weight[i, s] = clamp(l_weight, zero(l_weight), one(l_weight))
        end
    end 
end

interpolate_asset_policy!(ws::HouseholdWorkspace) = interpolate_asset_policy!(ws.lower_index, ws.lower_weight; ws.h.a_grid, ws.a_pol)


function stationary(h; R, T, w, kwargs...)
    ws = HouseholdWorkspace(; h, R, T, w)
    stationary!(ws; R, T, w, kwargs...)
    return ws
end 


stationary!(ws; R, T, w, kwargs...) = stationary!(ws.h.u, ws; R, T, w, kwargs...)


function stationary!(
    ::EZ, ws; 
    R, T, w, 
    max_iters = _MAX_ITERS, print_every = 100, 
    value_tol = _TOL,
    policy_tol = _TOL,
    pdf_tol = _TOL, 
    verbose = true
)

    (; h, v, η, v_tmp, η_tmp, a_tmp) = ws
    v_η_0 = (; v = v, η = η)
    v_η_1 = (; v = v_tmp, η = η_tmp)
    i = 1
    while true
        backwards_once!(h, v_η_1; next_values = v_η_0, R, T, w, a_tmp = a_tmp)
        dis1 = chebyshev(v_η_0.v, v_η_1.v)
        if dis1 < value_tol 
            # checked that policy converged too 
            dis2 = chebyshev(v_η_0.η, v_η_1.η)
            if dis2 < policy_tol
                verbose && println("Converged. iter $i;  R = $R, T = $T, w = $w, |v|: $dis1 , |η| : $dis2")
                break
            end
        end 
        if i > max_iters
            dis2 = (v_η_0.η, v_η_1.η)
            @warn("Max iters reached. DID NOT CONVERGE! R = $R, T = $T, w = $w,  |v|: $dis1 , |η| : $dis2")
            break
        end 
        verbose && (i % print_every == 1) && println("iter $i: $dis1") 
        i += 1
        v_η_0, v_η_1 = v_η_1, v_η_0
    end
    # v, η store the final answer
    ws.v .= v_η_0.v
    ws.η .= v_η_0.η

    asset_policy_given_η!(ws; R, T, w)
    interpolate_asset_policy!(ws)
    update_stationary_pdf!(ws; tol = pdf_tol)

end 


function stationary!(
    ::CRRA, ws; 
    R, T, w, 
    max_iters = _MAX_ITERS, print_every = 100, 
    policy_tol = _TOL,
    pdf_tol = _TOL, 
    verbose = true
)

    (; h, η, η_tmp, a_tmp) = ws
    η_0 = η
    η_1 = η_tmp
    i = 1
    while true
        backwards_once!(h, (; η = η_1); next_values = (; η = η_0), R, T, w, a_tmp = a_tmp)
        dis1 = chebyshev(η_0, η_1)
        if dis1 < policy_tol 
            verbose && println("Converged. iter $i;  R = $R, T = $T, w = $w, |η|: $dis1")
            break
        end 
        if i > max_iters
            @warn("Max iters reached. DID NOT CONVERGE! R = $R, T = $T, w = $w,  |η|: $dis1")
            break
        end 
        verbose && (i % print_every == 1) && println("iter $i: $dis1") 
        i += 1
        η_0, η_1  = η_1, η_0
    end
    # η stores the final answer
    ws.η .= η_0

    asset_policy_given_η!(ws; R, T, w)
    interpolate_asset_policy!(ws)
    update_stationary_pdf!(ws; tol = pdf_tol)
end 


function forward_pdf!(pdf_next; h, pdf, lower_index, lower_weight)
    P = h.P
    fill!(pdf_next, zero(eltype(pdf_next)))
    for i in axes(pdf, 1)
        @inbounds for s in axes(pdf, 2)
            ind = lower_index[i, s]
            weight = lower_weight[i, s]
            for s1 in axes(pdf, 2)
                mass = P[s, s1] * pdf[i, s]
                pdf_next[ind, s1] += mass * weight 
                pdf_next[ind + 1, s1] += mass - mass * weight
            end 
        end 
    end 
end 


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


asset_supply(a_grid, pdf) = @tullio s := a_grid[i] * pdf[i, s]


aggregate_c(h; R, w, T, a_pol, pdf) = @tullio threads=true tot :=  (R * h.a_grid[i] + T + labor_income(h.v; w = w * h.z_grid[s]) - a_pol[i, s]) * pdf[i, s]


function labor_supply(h::Household; w)
    (; v, z_grid, Pss) = h
    n = zero(eltype(z_grid))
    for (p, z) in zip(Pss, z_grid)
        n += p * z * labor(v; w = w * z)
    end
    return n
end


# Given a policy function, compute the individual consumption matrix
function consumption_alloc(h; R, w, T, a_pol)
    (; z_grid,  a_grid) = h
    c = similar(a_pol)
    @batch for s in axes(c, 2)
        for i in axes(c, 1)
            @inbounds c[i, s] = R * a_grid[i] + T + labor_income(h.v; w = w * z_grid[s]) - a_pol[i, s] 
        end
    end
    return c
end


is_pol_valid(a_pol, h) = all(a_pol .< last(h.a_grid))
is_pol_valid(e::StationaryEquilibrium) = is_pol_valid(e.ws.a_pol, e.h)

