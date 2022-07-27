

function _jacobian_item(h::Household)  
    v = Array{eltype(h.a_grid)}(undef, length(h.a_grid), length(h.z_grid))
    η = similar(v) 
    a_pol = similar(v)
    pdf = similar(v)
    lower_index = similar(v, Int)
    lower_weight = similar(v)
    a_tmp = similar(v)
    return (; v, η, a_pol, pdf, lower_index, lower_weight, a_tmp)
end 


function JacobianCache(ws::HouseholdWorkspace; cap_s, cap_t, R, T, w, ΔR, ΔT)
    path_up = StructArray(_jacobian_item(ws.h) for _ in 1:cap_t)
    path_down = StructArray(_jacobian_item(ws.h) for _ in 1:cap_t)
    cache = JacobianCache(path_up, path_down, cap_s, cap_t, ws, R, T, ΔR, ΔT)
    reset!(cache, ws; R, T, w, ΔR, ΔT)
    return cache
end
JacobianCache(e::StationaryEquilibrium; cap_s, cap_t, ΔR, ΔT) = JacobianCache(e.ws; R = 1 + e.r, e.T, e.w, cap_s, cap_t, ΔR, ΔT)


function reset!(cache::JacobianCache, ws::HouseholdWorkspace; R, T, w, ΔR, ΔT)
    (; cap_s, cap_t) = cache
    h = ws.h

    cache.ws = ws
    cache.R, cache.ΔR = R, ΔR
    cache.T, cache.ΔT = T, ΔT
    cache.up[1].pdf .= ws.pdf
    cache.down[1].pdf .= ws.pdf
    for (path, R1, T1) in [ (cache.up, R + ΔR, T + ΔT), (cache.down, R - ΔR, T - ΔT)]
        for i in cap_t:-1:1 
            v = (i >= cap_t) ? ws.v : path[i + 1].v 
            η = (i >= cap_t) ? ws.η : path[i + 1].η
            R0 = (i == cap_s) ? R1 : R
            T0 = (i == cap_s) ? T1 : T  
            backwards_once!(path[i].v, path[i].η; v, η, h, R = R0, T  = T0, w, path[i].a_tmp)
            asset_policy_given_η!(path[i].a_pol; h, path[i].η, R = R0, T = T0, w)
            interpolate_asset_policy!(path[i].lower_index, path[i].lower_weight; h.a_grid, path[i].a_pol)
        end
    end
end 
reset!(cache::JacobianCache, e::StationaryEquilibrium; ΔR, ΔT) = reset!(cache, e.ws; R = 1 + e.r, e.T, e.w, ΔR, ΔT)


function jacobian_column(s, ws::HouseholdWorkspace; R, T, w, cap_s, cap_t, ΔR, ΔT)
    cache = JacobianCache(ws; cap_s, cap_t, R, T, w, ΔR, ΔT)
    f = x -> jacobian_column!(cache, x)
    return f.(s)
end 
jacobian_column(s, e::StationaryEquilibrium; cap_s, cap_t, ΔR, ΔT) = jacobian_column(s, e.ws; R = 1 + e.r, e.T, e.w, cap_s, cap_t, ΔR, ΔT)
 

function jacobian_column!(cache::JacobianCache, s)
    (; ws, cap_s, cap_t, ΔR, ΔT) = cache 
    h = ws.h
    (s > cap_s || s < 1) && @error " s out of valid range .. exiting."

    jac = similar(ws.a_pol, length(cache.up))

    for t in axes(cache.up, 1)
        adiff = zero(eltype(h.a_grid)) 
        @inbounds for (path, signΔ) in [ (cache.up, 1), (cache.down, -1)]
            lower_weight = t > s ? ws.lower_weight : path[cap_s - (s - t)].lower_weight
            lower_index = t > s ? ws.lower_index : path[cap_s - (s - t)].lower_index
            a = asset_supply(h.a_grid, path[t].pdf)
            adiff += a * signΔ 
            (t < cap_t) && forward_pdf!(path[t+1].pdf; h, path[t].pdf, lower_index, lower_weight)
        end
        jac[t] = adiff / (2 * (ΔR + ΔT))
    end 
    return jac 
end 


function jacobian(ws::HouseholdWorkspace; R, T, w, ΔR, ΔT, cap_t)
    cache = JacobianCache(ws; cap_t, cap_s = cap_t, R, T, w, ΔR, ΔT)
    h = ws.h

    jac = similar(ws.a_pol, length(cache.up), length(cache.up))

    for s in axes(cache.up, 1)
        for t in axes(cache.up, 1)
            adiff = zero(eltype(h.a_grid)) 
            @inbounds for (path, signΔ) in ((cache.up, 1), (cache.down, -1))
                lower_weight = t > s ? ws.lower_weight : path[cap_t - (s - t)].lower_weight
                lower_index = t > s ? ws.lower_index : path[cap_t - (s - t)].lower_index
                a = asset_supply(h.a_grid, path[t].pdf)
                adiff += a * signΔ 
                (t < cap_t) && forward_pdf!(path[t+1].pdf; h, path[t].pdf, lower_index, lower_weight)
            end
            jac[s, t] = adiff / (2 * (ΔR + ΔT))
        end 
    end 
    return jac 
end 
jacobian(e::StationaryEquilibrium; cap_t, ΔR, ΔT) = jacobian(e.ws; R = 1 + e.r, e.T, e.w, cap_t, ΔR, ΔT)


function pv_elasticities!(cache::JacobianCache, s; Rk)
    (; cap_s, cap_t, ΔR, ws, R) = cache 
    h = ws.h
    (s > cap_s || s < 1) && @error " s out of valid range .. exiting."

    pv = zero(R)
    for t in axes(cache.up, 1)
        adiff = zero(R) 
        @inbounds for (path, signΔ) in [ (cache.up, 1), (cache.down, -1)]
            lower_weight = t > s ? ws.lower_weight : path[cap_s - (s - t)].lower_weight
            lower_index = t > s ? ws.lower_index : path[cap_s - (s - t)].lower_index
            a = asset_supply(h.a_grid, path[t].pdf)
            adiff += a * signΔ
            (t < cap_t) && forward_pdf!(path[t+1].pdf; h, path[t].pdf, lower_index, lower_weight)
        end
        pv += adiff * Rk^(- t)
    end 
    aa = asset_supply(h.a_grid, ws.pdf)
    return pv * (Rk - R) / (2 * ΔR) / aa * Rk^s
end 


function pv_elasticities(s, e::StationaryEquilibrium; cap_s, cap_t, ΔR = 1e-4)
    cache = JacobianCache(e; cap_s, cap_t, ΔR, ΔT = 0.0)
    Rk = 1 + mpk(e.t; k = e.k, n = e.n) - e.t.δ
    f = x -> pv_elasticities!(cache, x; Rk)
    return f.(s)
end 

