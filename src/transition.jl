
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


function solve_transition(
    e_init, e_final; k_path, b_path, r_path = nothing, verbose = true,
    nlsolve_kwargs = nothing
)
    nlsolve_baseline_kwargs = (
        ftol = _ZERO_FTOL,
        show_trace = verbose,
        method = :anderson,
        m = 10,
        iterations = 300,
        beta = -0.5  # adjustment parameter for the fixed point algorithm
    )
    
    nlsolve_kwargs_merged = isnothing(nlsolve_kwargs) ?             
        nlsolve_baseline_kwargs :
        merge(nlsolve_baseline_kwargs, nlsolve_kwargs)

    path = _initialize_path(e_init, e_final; k_path)

    if isnothing(r_path) 
        r_path_rest = [(e_final.r) for _ in 1:length(k_path) - 1]
    else 
        r_path_rest = r_path[2:end]
    end 
    r_path_rest[end] = e_final.r

    f! = (F, x) -> _residuals!(F, x; path, e_init, k_path, b_path)
    r_sol = nlsolve(f!, r_path_rest; nlsolve_kwargs_merged...)
    _residuals!(similar(r_sol.zero), r_sol.zero; path, e_init, k_path, b_path)
    
    return path 
end


function _residuals!(residuals, r_path_rest; path, e_init, k_path, b_path)
    h = e_init.h
    t = e_init.t
    w = e_init.w
    n0 = e_init.n
    r0 = e_init.r 
    k0 = e_init.k

    r_min = - t.δ + eps()  # minimum interest rate
    r_max = 1 / h.u.β  # a maximum -- could be relaxed
    min_T = minimum_feasible_transfer(h, w) + eps() # minimum transfer
    
    @inbounds for i in reverse(1:length(r_path_rest))
        r = (i == 1) ? r0 : r_path_rest[i-1]
        r = clamp(r, r_min, r_max)
        T_ = get_T(t; 
            b = b_path[i], bprime = b_path[i+1],
            r = r, k = k_path[i], 
            k0, r0, n0)
        T = max(T_, min_T)   
        path.T[i] = T
        path.r[i] = r
        R = 1 + r
        backwards_once!(path[i].v, path[i].η; path[i+1].v, path[i+1].η, h, R, T, w)
        asset_policy_given_η!(path[i].a_pol; h, path[i].η, R, T, w)
        interpolate_asset_policy!(path[i].lower_index, path[i].lower_weight; h.a_grid, path[i].a_pol)
    end

    @inbounds for i in 1:length(r_path_rest) + 1 
        a = asset_supply(h.a_grid, path[i].pdf)
        path.a[i] = a
        path.b[i] = b_path[i]
        path.k[i] = k_path[i]
        if i > 1
            residuals[i-1] =  a - k_path[i] - b_path[i]
        end
        (i == length(r_path_rest) + 1) && break
        forward_pdf!(path[i+1].pdf; h, path[i].pdf, path[i].lower_index, path[i].lower_weight)
    end
end 