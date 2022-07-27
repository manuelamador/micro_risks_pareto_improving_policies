
# Household problem and aggregation 

function minimum_feasible_transfer(h, w)
    (; z_grid) = h
    wi = w .* z_grid
    y = map(x -> labor_income(h.v; w = x), wi)
    vofn = map(x -> disutility_given_w(h.v; w = x), wi)
    return -minimum(y .- vofn)
end


# Household's euler equation 

function backwards_euler_x(u::EZ, π_vec, v_vec, η_vec)
    ξ = get_inverse_ies(u)
    γ = get_ra(u)
    β = get_β(u)
    Ev = Evx = zero(eltype(v_vec))
    @inbounds for i in eachindex(π_vec)
        v = v_vec[i]
        η = η_vec[i]
        Evx += π_vec[i] * (v)^((ξ - γ)) * (η)^(-ξ)
        Ev += π_vec[i] * u.risk(v)
    end 
    x = (β * inverse(u.risk, Ev^((γ - ξ))) * Evx)^(-1 / ξ)
    return x
end


function asset_policy_given_η!(a_pol; h, η, R, T, w)
    ξ = get_inverse_ies(h.u) 
    @tullio a_pol[i, s] = R * h.a_grid[i] + T + labor_income(h.v; w = w * h.z_grid[s]) -  η[i, s] * R^(1/ξ) - disutility_given_w(h.v; w = w * h.z_grid[s])
end


# function asset_policy_given_η!(a_pol; h, η, R, T, w)
#     ξ = get_inverse_ies(h.u) 
#     @batch for s in eachindex(h.z_grid)
#         z = h.z_grid[s]
#         dis = disutility_given_w(h.v; w = w * z)
#         lab = labor_income(h.v; w = w * z)
#         @inbounds for i in eachindex(h.a_grid)
#             x = η[i, s] * R^(1/ξ)
#             c = x + dis
#             a = h.a_grid[i]
#             a_prime = R * a + T + lab - c 
#             a_pol[i, s] = a_prime 
#         end 
#     end 
# end


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


function backwards_once!(v_now, η_now; v, η, h, R, T, w, a_tmp = similar(v_now))   
    v_next, η_next = v, η
    (; u, P) = h
    ξ = get_inverse_ies(u)
    β = get_β(u)

    @batch for s in eachindex(h.z_grid)
        z = h.z_grid[s]
        dis = disutility_given_w(h.v; w = w * z)
        lab = labor_income(h.v; w = w * z)

        @inbounds for i in eachindex(h.a_grid)
            x = backwards_euler_x(h.u, view(P, s, :), view(v_next, i, :), view(η_next, i, :))
            c = x + dis
            a_next = h.a_grid[i]
            a_tmp[i, s]  = (c + a_next - lab - T) / R
        end

        @inbounds for i in eachindex(h.a_grid)
            a = h.a_grid[i]
            a_prime = interp1D(a, view(a_tmp, :, s), h.a_grid)
            a_prime = clamp(a_prime, first(h.a_grid), last(h.a_grid))
            max_feasible_a_prime = R * a + T + lab - dis
            
            # a_prime = min(max_feasible_a_prime, a_prime)
            x = max_feasible_a_prime - a_prime
            Ev = zero(a)
            for s2 in axes(v_next, 2)
                v_prime = interp1D(a_prime, h.a_grid, view(v_next, :, s2))
                Ev += P[s, s2] * u.risk(v_prime)
            end
            v_now[i, s] =
                inverse(u.temporal, (1 - β) * u.temporal(x) + β * u.temporal(inverse(u.risk, Ev)))

            η_now[i, s] = R^(-1/ξ) * x
        end
    end    
end


function stationary( 
    h; R, T, w, 
    max_iters = _MAX_ITERS, 
    print_every = 100, 
    value_tol = _TOL,
    policy_tol = _TOL,
    pdf_tol = _TOL, 
    verbose = true
)
    ws = HouseholdWorkspace(; h, R, T, w)
    stationary!(ws; R, T, w, max_iters, print_every, value_tol, policy_tol, pdf_tol, verbose)
    return ws
end 


function stationary!(
    ws; 
    R, T, w, 
    max_iters = _MAX_ITERS, print_every = 100, 
    value_tol = _TOL,
    policy_tol = _TOL,
    pdf_tol = _TOL, 
    verbose = true
)

    (; h, v, η, v_tmp, η_tmp, a_tmp) = ws
    v_0, η_0 = v, η
    v_1, η_1 = v_tmp, η_tmp
    i = 1
    while true
        backwards_once!(v_1, η_1; h, v = v_0, η = η_0, R, T, w, a_tmp = a_tmp)
        dis1 = chebyshev(v_0, v_1)
        if dis1 < value_tol 
            # checked that policy converged too 
            dis2 = chebyshev(η_0, η_1)
            if dis2 < policy_tol
                verbose && println("Converged. iter $i;  R = $R, T = $T, w = $w, |v|: $dis1 , |η| : $dis2")
                break
            end
        end 
        if i > max_iters
            dis2 = chebyshev(η_0, η_1) 
            @warn("Max iters reached. DID NOT CONVERGE! R = $R, T = $T, w = $w,  |v|: $dis1 , |η| : $dis2")
            break
        end 
        verbose && (i % print_every == 1) && println("iter $i: $dis1") 
        i += 1

        v_0, v_1  = v_1, v_0
        η_0, η_1  = η_1, η_0
    end
    # v, η store the final answer
    ws.v .= v_0
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


# function asset_supply(a_grid, pdf)
#     tot = zero(eltype(pdf))
#     @inbounds for s in axes(pdf, 2)
#         for i in axes(pdf, 1)
#             tot += a_grid[i] * pdf[i, s]
#         end
#     end 
#     return tot
# end

asset_supply(a_grid, pdf) = @tullio s := a_grid[i] * pdf[i, s]


aggregate_c(h; R, w, T, a_pol, pdf) = @tullio threads=true tot :=  (R * h.a_grid[i] + T + labor_income(h.v; w = w * h.z_grid[s]) - a_pol[i, s]) * pdf[i, s]


# function aggregate_c(h; R, w, T, a_pol, pdf)
#     (; z_grid,  a_grid) = h
#     tot = zero(eltype(pdf))
#     @inbounds for s in axes(pdf, 2)
#         for i in axes(pdf, 1)
#             tot += (R * a_grid[i] + T + labor_income(h.v; w = w * z_grid[s]) - a_pol[i, s]) * pdf[i, s]
#         end
#     end
#     return tot
# end


# labor_supply(h::Household; w) = @tullio s := h.Pss[i] * h.z_grid[i] * labor(h.v, w = w * h.z_grid[i])

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

