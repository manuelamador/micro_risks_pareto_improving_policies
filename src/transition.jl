
function _transition_item(h::Household)  
    v = Array{eltype(h.a_grid)}(undef, length(h.a_grid), length(h.z_grid))
    η = similar(v) 
    a_pol = similar(v)
    pdf = similar(v)
    lower_index = similar(v, Int)
    lower_weight = similar(v)
    a_tmp = similar(v)
    r = T = a = k = b = zero(eltype(h.a_grid)) 
    return (; h, v, η, a_pol, pdf, lower_index, lower_weight, a_tmp, r, T, a, k, b)
end 


function _initialize_path(e_init, e_final; k_path)
    h = e_init.h

    lst = StructArray(_transition_item(h) for _ in k_path)
    lst[1].pdf .= e_init.ws.pdf

    lst[end].v .= e_final.ws.v 
    lst[end].η .= e_final.ws.η
    lst[end].a_pol .= e_final.ws.a_pol
    lst[end].lower_index .= e_final.ws.lower_index
    lst[end].lower_weight .= e_final.ws.lower_weight

    lst.r[end] = e_final.r
    lst.b[end] = e_final.b
    lst.k[end] = e_final.k
    lst.a[end] = e_final.a

    return lst 
end 

# OLD SOLVER -- COULD SERVE AS BACKUP 
# 
# function solve_transition(
#     e_init, e_final; k_path, b_path, r_path = nothing, verbose = true,
#     nlsolve_kwargs = nothing
# )
#     nlsolve_baseline_kwargs = (
#         ftol = _ZERO_FTOL,
#         show_trace = verbose,
#         method = :anderson,
#         m = 10,
#         iterations = 300,
#         beta = -0.5  # adjustment parameter for the fixed point algorithm
#     )
    
#     nlsolve_kwargs_merged = isnothing(nlsolve_kwargs) ?             
#         nlsolve_baseline_kwargs :
#         merge(nlsolve_baseline_kwargs, nlsolve_kwargs)

#     path = _initialize_path(e_init, e_final; k_path)

#     if isnothing(r_path) 
#         r_path_rest = [(e_final.r) for _ in 1:length(k_path) - 1]
#     else 
#         r_path_rest = r_path[2:end]
#     end 
#     r_path_rest[end] = e_final.r

#     f! = (F, x) -> _residuals!(F, x; path, e_init, k_path, b_path)
#     r_sol = nlsolve(f!, r_path_rest; nlsolve_kwargs_merged...)
#     _residuals!(similar(r_sol.zero), r_sol.zero; path, e_init, k_path, b_path)
    
#     return path 
# end



function solve_transition(
    e_init, e_final; k_path, b_path, r_path = nothing, verbose = true,
    print_every = 10, max_iters = _ZERO_MAX_ITERS, tol = _ZERO_FTOL
)    
    path = _initialize_path(e_init, e_final; k_path)

    if isnothing(r_path) 
        r_sol = [(e_final.r) for _ in 1:length(k_path) - 1]
    else 
        r_sol = r_path[1:length(k_path) - 1]
    end 

    r_sol[1] = e_init.r

    jac_R = jacobian(e_init; cap_t = length(r_sol) + 1, ΔR = 1e-4, ΔT = 0.0)[2:end, 2:end]
    jac_T = jacobian(e_init; cap_t = length(r_sol) + 1, ΔT = 1e-4, ΔR = 0.0)[2:end, 2:end]
    jac = similar(jac_T)
    for s in axes(jac, 1)
        for t in axes(jac, 2)
            jac[s, t] = jac_R[s, t] - jac_T[s, t] * (b_path[t] + k_path[t]) 
        end
    end 

    residuals = fill(zero(e_final.a), length(r_sol))
    iter = 1

    f! = (F, x) -> _residuals!(F, x; path, e_init, e_final, k_path, b_path)
    while true 
        f!(residuals, r_sol)
        r_sol .= r_sol .- inv(jac) * residuals
        dis = maximum((abs(x) for x in residuals))
        verbose && (iter % print_every == 1) && println("Iter: $iter, error: $dis")
        dis < tol && break
        iter > max_iters && (warn("Did not converge"); break)
        iter += 1 
    end  

    return path 
end


function _residuals!(residuals, r_path_rest; path, e_init, e_final, k_path, b_path)
    h = e_init.h
    t = e_init.t
    w = e_init.w
    n0 = e_init.n
    r0 = e_init.r 
    k0 = e_init.k

    r_min = - t.δ + eps()  # minimum interest rate
    r_max = 1 / h.u.β  # a maximum -- could be relaxed
    min_T = minimum_feasible_transfer(h, w) + eps() # minimum transfer
    
    @inbounds for i in reverse(eachindex(path))
        r = (i == firstindex(path)) ? r0 : r_path_rest[i - 1]
        r = clamp(r, r_min, r_max)
        T_ = get_T(t; 
            b = b_path[i], bprime = b_path[i+1],
            r = r, k = k_path[i], 
            k0, r0, n0)
        T = max(T_, min_T)   
        path.T[i] = T
        path.r[i] = r
        R = 1 + r
        next_values = i == lastindex(path) ? e_final.ws : path[i+1]
        backwards_once!(h, path[i]; next_values = next_values, R, T, w, path[i].a_tmp)
        asset_policy_given_η!(path[i].a_pol; h, path[i].η, R, T, w)
        interpolate_asset_policy!(path[i].lower_index, path[i].lower_weight; h.a_grid, path[i].a_pol)
    end

    @inbounds for i in eachindex(path) 
        a = asset_supply(h.a_grid, path[i].pdf)
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