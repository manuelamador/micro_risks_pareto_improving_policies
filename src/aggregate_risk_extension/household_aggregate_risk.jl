
#######################################################
#            AGGREGATE RISK FUNCTIONS                 #
#######################################################

# This file contains the functions that solve the household bellman equation, aggregate results, as well as compute the stationary solution given R, w, T.
# These functions use only information on the household problem and prices and transfers; and do not use information about technology or fiscal policy. 


"""
    backwards_once_t0!(h, current_values; next_values_H, next_values_L, R, T, w, a_tmp, probH)

Modifies `current_values` (i.e period-0 values) using the Euler equation and the continuation values in `next_values_H` and `next_values_L`.

# Arguments 
- `h`: household struct.
- `current_values`: struct containing the matrices `η` and `v` to be computed (the latter for the case of Epstein-Zin).
- `next_values_H`: struct containing the next period's values when the realization of the aggregate prod. shock is high. It includes the matrices `η` and `v` (the latter for the case of Epstein-Zin).
- `next_values_L`: struct containing the next period's values when the realization of the aggregate prod. shock is low. It includes the matrices `η` and `v` (the latter for the case of Epstein-Zin).
- `R`: current gross interest rate.
- `T`: current transfer.
- `w`: current wage.
- `a_tmp`: temporary array.
- `probH`: Probability of the high aggregate productivity shock

The values of `current_values.η` and `current_values.v` are modified.
"""

function backwards_once_t0!(h, current_values; next_values_H, next_values_L, R, T, w, a_tmp = similar(current_values.η), probH)   
    η_now = current_values.η
    (; u) = h
    P′ = h.Pprime
    ξ = get_inverse_ies(u)
    γ = get_ra(u)
    β = get_β(u)
    R_ = R^(-1/ξ)
    a_min = first(h.a_grid)
    a_max = last(h.a_grid)

    Threads.@threads for s in eachindex(h.z_grid)
        z = h.z_grid[s]
        add =  labor_income(h.v; w = w * z)  + T - disutility_given_w(h.v; w = w * z)
        for i in eachindex(h.a_grid)
            p_vec = view(P′, :, s)
            v_vec_H = view(next_values_H.v, i, :)
            η_vec_H = view(next_values_H.η, i, :)
            v_vec_L = view(next_values_L.v, i, :)
            η_vec_L = view(next_values_L.η, i, :)
            Ev = Evx = zero(eltype(v_vec_H))
            for (v_H, η_H, v_L, η_L, p) in zip(v_vec_H, η_vec_H, v_vec_L, η_vec_L, p_vec)
                Evx += p * (probH * (v_H)^((ξ - γ)) * (η_H)^(-ξ) + (1-probH) * (v_L)^((ξ - γ)) * (η_L)^(-ξ))
                Ev += p * (probH * u.risk(v_H) + (1-probH) * u.risk(v_L))
            end 
            x = (β * inverse(u.risk, Ev^((γ - ξ))) * Evx)^(-1/ξ)
            a_next = h.a_grid[i]
            a_tmp[i, s]  = (x + a_next - add) / R
        end

        for i in eachindex(h.a_grid)
            a = h.a_grid[i]
            a_prime = interp1D(a, view(a_tmp, :, s), h.a_grid)
            a_prime = clamp(a_prime, a_min, a_max)
            x = R * a + add - a_prime
            η_now[i, s] = R_ * x
            Ev = zero(a_prime)
            for s2 in axes(next_values_H.v, 2) # from 1 to 10
                v_prime_H = interp1D(a_prime, h.a_grid, view(next_values_H.v, :, s2))
                v_prime_L = interp1D(a_prime, h.a_grid, view(next_values_L.v, :, s2))
                Ev += P′[s2, s] * (probH * u.risk(v_prime_H) + (1-probH) * u.risk(v_prime_L))
            end
            current_values.v[i, s] = inverse(u.temporal, (1 - β) * u.temporal(x) + β * u.temporal(inverse(u.risk, Ev)))
        end
    end    
end

#######################################################
# Portfolio problem
#######################################################

"""
    shadow_rates!(h,t; path)

Calculates the Euler equation residuals and the "shadow rate" (i.e, the ratio between the period-0 and expected discounted period-1 marginal utilities of consumption) for each household. 
"""

function shadow_rates!(h,t; path)    
    (; u) = h
    P′ = h.Pprime
    ξ = get_inverse_ies(u)
    γ = get_ra(u)
    β = get_β(u)

    r = copy(path[1,1].r)
    R = 1 + r 
    T = copy(path.T[1,1])
    w = mpl_from_mpk(t; mpk = r + t.δ)
    R_ = R^(-1/ξ)
    a_min = first(h.a_grid)
    a_max = last(h.a_grid)
    probH = 0.5
    current_values_p = deepcopy(path[1,1])
    next_values_H = deepcopy(path[2,1]) 
    next_values_L = deepcopy(path[2,2])

    a_tmp = zeros(length(h.a_grid),length(h.z_grid))
    η_now = similar(a_tmp)
    x_now = similar(a_tmp)
    rhs_euler = similar(a_tmp)
    rhs_euler_s = similar(a_tmp)
    a_prime_pol = similar(a_tmp)

    for s in eachindex(h.z_grid)
        z = h.z_grid[s]
        add =  labor_income(h.v; w = w * z)  + T - disutility_given_w(h.v; w = w * z)

        for i in eachindex(h.a_grid)
            p_vec = view(P′, :, s)
            v_vec_H = view(next_values_H.v, i, :)
            η_vec_H = view(next_values_H.η, i, :)
            v_vec_L = view(next_values_L.v, i, :)
            η_vec_L = view(next_values_L.η, i, :)
            Ev = Evx = zero(eltype(v_vec_H))
            for (v_H, η_H, v_L, η_L, p) in zip(v_vec_H, η_vec_H, v_vec_L, η_vec_L, p_vec)
                Evx += p * (probH * (v_H)^((ξ - γ)) * (η_H)^(-ξ) + (1-probH) * (v_L)^((ξ - γ)) * (η_L)^(-ξ))
                Ev += p * (probH * u.risk(v_H) + (1-probH) * u.risk(v_L))
            end 
            x = (β * inverse(u.risk, Ev^((γ - ξ))) * Evx)^(-1/ξ)
            a_next = h.a_grid[i]
            a_tmp[i, s]  = (x + a_next - add) / R
        end

        for i in eachindex(h.a_grid)
            a = h.a_grid[i]
            a_prime = interp1D(a, view(a_tmp, :, s), h.a_grid)
            a_prime = clamp(a_prime, a_min, a_max)
            a_prime_pol[i,s] = a_prime
            x = R * a + add - a_prime
            x_now[i, s] = x
            η_now[i, s] = R_ * x

            Ev = zero(a_prime)
            Evx = zero(a_prime)
            Evx_s = zero(a_prime)
            for s2 in axes(next_values_H.v, 2)
                v_prime_H = interp1D(a_prime, h.a_grid, view(next_values_H.v, :, s2))
                v_prime_L = interp1D(a_prime, h.a_grid, view(next_values_L.v, :, s2))
                η_prime_H = interp1D(a_prime, h.a_grid, view(next_values_H.η, :, s2))
                η_prime_L = interp1D(a_prime, h.a_grid, view(next_values_L.η, :, s2))
                x_prime_H = η_prime_H ./ (1+path[2,1].r)^(-1/ξ)
                x_prime_L = η_prime_L ./ (1+path[2,2].r)^(-1/ξ)
                Ev += P′[s2, s] * (probH * u.risk(v_prime_H) + (1-probH) * u.risk(v_prime_L))
                Evx += P′[s2, s] * (probH * (v_prime_H)^((ξ - γ)) * (η_prime_H)^(-ξ) + (1-probH)* (v_prime_L)^((ξ - γ)) * (η_prime_L)^(-ξ))
                Evx_s += P′[s2, s] * (probH * (v_prime_H)^((ξ - γ)) * (x_prime_H)^(-ξ) + (1-probH)* (v_prime_L)^((ξ - γ)) * (x_prime_L)^(-ξ))
            end
            rhs_euler[i,s] = (β * inverse(u.risk, Ev^((γ - ξ))) * Evx)^(-1/ξ)
            rhs_euler_s[i,s] = β * inverse(u.risk, Ev^((γ - ξ))) * Evx_s
            current_values_p.v[i, s] = inverse(u.temporal, (1 - β) * u.temporal(x) + β * u.temporal(inverse(u.risk, Ev)))
        end
    end
    euler_resids = x_now .- rhs_euler
    shadow_rates = (x_now.^(- ξ))./rhs_euler_s

    return euler_resids, shadow_rates
end

"""
    resid_portfolio!(theta,a_2; h, port, R_b, R_k_h, R_k_l, s1, probH)

Calculates the residual of the portfolio equation (see Section D.7 in the Appendix)

# Arguments 
- `theta`: Share of total savings allocated to risky capital
- `a_2`: Total savings in the second period
- `h`: Household struct.
- `port`: Portfolio problem struct.
- `R_b`: Bond "risk-free" gross interest rate.
- `R_k_h`: Return on capital when aggregate productivity is high.
- `R_k_l`: Return on capital when aggregate productivity is low.
- `s1`: Index of the idiosyncratic productivity in the first period
- `probH`: Probability of the high aggregate productivity shock
"""

function resid_portfolio!(theta,a_2; h, port, R_b, R_k_h, R_k_l, s1, probH)
    P′ = h.Pprime
    #s1 is the position of z in period 1
    cih_H = ((1-theta)*R_b + theta*R_k_h)*a_2
    cih_L = ((1-theta)*R_b + theta*R_k_l)*a_2

    Ev_k = 0.0
    Ev_b = 0.0
    for s2 in axes(P′, 1) # for each z tomorrow
        u_prime_H = port.lst_u_c[s2,1](cih_H) #u'(c_2) where c_2(z_2,cih_2). Here s2 refers to z_2, 1 refers to high, and cih_H is the cash-in-hand in period 2
        u_prime_L = port.lst_u_c[s2,2](cih_L) #u'(c_2) where c_2(z_2,cih_2). Here s2 refers to z_2, 2 refers to low, and cih_L is the cash-in-hand in period 2
        Ev_k += P′[s2, s1] * (probH * u_prime_H * R_k_h + (1-probH) * u_prime_L *R_k_l) #Ev_k(z_2,cih_2) where the expectation is taken when z_1 in period 1
        Ev_b += P′[s2, s1] * (probH * u_prime_H + (1-probH) * u_prime_L) #Ev_b(z_2,cih_2)
    end
    return R_b*Ev_b-Ev_k 
end

"""
    get_theta(a_2; h, port, R_b, R_k_h, R_k_l, s1, probH)

Finds the share of total savings allocated to risky capital (i.e. "theta") that makes the residual of the portfolio equation equal to zero.

# Arguments 
- `a_2`: Total savings in the second period
- `h`: Household struct.
- `port`: Portfolio problem struct.
- `R_b`: Bond "risk-free" gross interest rate.
- `R_k_h`: Return on capital when aggregate productivity is high.
- `R_k_l`: Return on capital when aggregate productivity is low.
- `s1`: Index of the idiosyncratic productivity in the first period
- `probH`: Probability of the high aggregate productivity shock
"""

function get_theta(a_2; h, port, R_b, R_k_h, R_k_l, s1, probH)
    f = (theta) -> resid_portfolio!(theta,a_2; h, port, R_b, R_k_h, R_k_l, s1, probH)
    theta = find_zero(f, 0.5)
    return theta
end

"""
    theta_tmp!(h, port; R_b, R_k_h, R_k_l, probH)

Solves for the share of total savings allocated to risky capital (i.e. "theta") for each point in the state-space (z,a,Z).
It takes into account that the share must be between 0 and 1.    

# Arguments 
- `h`: Household struct.
- `port`: Portfolio problem struct.
- `R_b`: Bond "risk-free" gross interest rate.
- `R_k_h`: Return on capital when aggregate productivity is high.
- `R_k_l`: Return on capital when aggregate productivity is low.
- `probH`: Probability of the high aggregate productivity shock

The value of `port.theta_star` is modified.
"""
function theta_tmp!(h, port; R_b, R_k_h, R_k_l, probH)    
    for i in eachindex(h.a_grid)
        for m in eachindex(h.z_grid)
            # Here resids_mat includes R_b_2*Ev_b(z_2,cih_2)-Ev_k(z_2,cih_2). Note that resids_mat[i,m,1] is resid(i=a_2,m=z_1,high). From the perspective of period 1, with state z_1, and choosing a_2.
            port.resids_mat[i,m,1] = resid_portfolio!(0.0, h.a_grid[i] ; h, port, R_b=R_b , R_k_h=R_k_h , R_k_l=R_k_l , s1 = m, probH) #Evaluate theta=0
            port.resids_mat[i,m,2] = resid_portfolio!(1.0, h.a_grid[i] ; h, port, R_b=R_b , R_k_h=R_k_h , R_k_l=R_k_l , s1 = m, probH) #Evaluate theta=1
        end
    end
    for i in eachindex(h.a_grid)
        for m in eachindex(h.z_grid)
            if minimum(port.resids_mat[i,m,:]) > 0.0
                port.theta_star[i,m] = 0.0
            elseif maximum(port.resids_mat[i,m,:]) < 0.0
                port.theta_star[i,m] = 1.0
            else
                port.theta_star[i,m]= get_theta(h.a_grid[i] ; h, port, R_b=R_b , R_k_h=R_k_h , R_k_l=R_k_l , s1 = m, probH) # Here, theta_star(a_2,z_1). 
            end
        end
    end
end

"""
    backwards_once_t0_port!(h, current_values; port, R, R_b, R_k_h, R_k_l, T, w, a_tmp, probH)

Modifies `current_values` (i.e period-0 values) using the Euler equation and the continuation values that results from solving the portfolio problem.

# Arguments 
- `h`: household struct.
- `current_values`: struct containing the matrices `η` and `v` to be computed (the latter for the case of Epstein-Zin).
- `port`: Portfolio problem struct.
- `R`: Period-0 (initial steady state) gross interest rate.
- `R_b`: Bond "risk-free" gross interest rate.
- `R_k_h`: Return on capital when aggregate productivity is high.
- `R_k_l`: Return on capital when aggregate productivity is low.
- `T`: current transfer.
- `w`: current wage.
- `a_tmp`: temporary array.
- `probH`: Probability of the high aggregate productivity shock

The value of `current_values.v` is modified.
"""

function backwards_once_t0_port!(h, current_values; port, R, R_b, R_k_h, R_k_l, T, w, a_tmp = similar(current_values.η), probH)   
    (; u) = h
    P′ = h.Pprime
    ξ = get_inverse_ies(u)
    γ = get_ra(u)
    β = get_β(u)

    a_min = first(h.a_grid)
    a_max = last(h.a_grid)

    theta_star = port.theta_star # theta_star(a_2,z_1)
    lst_v = port.lst_v # Interpolation v(cih_2,z_2). Sintaxis: port.lst_v[l,j](cih_2) gives the interpolated value v_2(cih_2,z_2,A_2) (at the particular point cih_2, where l=z_2 and j=A_2) 
    lst_x = port.lst_x # Interpolation x(cih_2,z_2). Sintaxis: port.lst_x[l,j](cih_2) gives the interpolated value x_2(cih_2,z_2,A_2) (at the particular point cih_2, where l=z_2 and j=A_2)

    for s in eachindex(h.z_grid) # For each z_1. Here @batch must not be used
        z = h.z_grid[s]
        add =  labor_income(h.v; w = w * z)  + T - disutility_given_w(h.v; w = w * z)

        for i in eachindex(h.a_grid) # For each a_2
            #x = _backwards_euler_x_helper(h.u, P′, i, s, next_values)
            R_H = (1-theta_star[i,s])*R_b + theta_star[i,s]*R_k_h # Note theta_star[i,s] here is theta_star(a_2,z_1)
            R_L = (1-theta_star[i,s])*R_b + theta_star[i,s]*R_k_l # Note theta_star[i,s] here is theta_star(a_2,z_1)
            cih_H = R_H * h.a_grid[i] #cash-in-hand in t=2 (cih_2) in z_1 and a_2, if high
            cih_L = R_L * h.a_grid[i] #cash-in-hand in t=2 (cih_2) in z_1 and a_2, if low
            Ev = Evx = zero(eltype(R_H))
            for s2 in eachindex(h.z_grid) # For each z_2
                v_H = lst_v[s2,1](cih_H) # This is value v_2(cih_2,z_2,A_2) when A_2 is high, z_1 = s, a_2 = i (and thus, cih_2=cih_H)
                v_L = lst_v[s2,2](cih_L) # This is value v_2(cih_2,z_2,A_2) when A_2 is low, z_1 = s, a_2 = i (and thus, cih_2=cih_L)
                x_H = lst_x[s2,1](cih_H) # This is x_2(cih_2,z_2,A_2) when A_2 is high, z_1 = s, a_2 = i (and thus, cih_2=cih_H)
                x_L = lst_x[s2,2](cih_L) # This is x_2(cih_2,z_2,A_2) when A_2 is low, z_1 = s, a_2 = i (and thus, cih_2=cih_L)
                Evx += P′[s2, s] * (probH * (v_H)^((ξ - γ)) * R_H * (x_H)^(-ξ) + (1-probH) * (v_L)^((ξ - γ)) * R_L * (x_L)^(-ξ)) 
                Ev += P′[s2, s] * (probH * u.risk(v_H) + (1-probH) * u.risk(v_L))
            end 
            x = (β * inverse(u.risk, Ev^((γ - ξ))) * Evx)^(-1/ξ) # This is x_1 = c_1-v(n_1) given a_2 and z_1. Note (a_2,z_1) implies theta_star(a_2,z_1), which implies a return (r_H or r_L) and a cash-in-hand (cih_H or cih_L)
            a_next = h.a_grid[i] # Level of assets a_2
            a_tmp[i, s]  = (x + a_next - add) / R # c_1 - v(n_1) + a_2 - w_1*z_1*n_1 - T_1 + v(n_1) = R_1 a_1. I recover a_1(a_2,z_1) where a_2 is in the grid (and a_1 not necessarily)
        end

        for i in eachindex(h.a_grid)
            a = h.a_grid[i] # a here refers to a_1 (now in the grid) 
            a_prime = interp1D(a, view(a_tmp, :, s), h.a_grid) # Given s=z_1, take the vector a_1(a_2,z_1) and interpolate a_2. 
            a_prime = clamp(a_prime, a_min, a_max) # Given a_1, I know policy a_2 (standard policy function)
            theta_itp = interp1D(a_prime, h.a_grid, view(theta_star, :, s)) # From theta_star(a_2,z_1) and a_2(a_1,z_1) from previous line, get theta(a_1,z_1) through interpolation  
            x = R * a + add - a_prime # Recover x_1(a_1,z_1) from the budget constraint
            current_values.a_pol[i,s] = a_prime # Substitute for "asset_policy_given_η!".
            port.theta_pol[i,s] = theta_itp 
            port.k_pol[i,s] = theta_itp*a_prime
            port.b_pol[i,s] = (1.0 -theta_itp)*a_prime
            
            #_extra_updates!(u, current_values, next_values, h, a_prime, P′, x, i, s)
            R_H = (1-theta_itp)*R_b + theta_itp*R_k_h # Note this is R_H when the current state is (a_1,z_1)
            R_L = (1-theta_itp)*R_b + theta_itp*R_k_l # Note this is R_L when the current state is (a_1,z_1)
            cih_H = R_H * a_prime #cash-in-hand in t=2 (cih_2) in z_1 and a_1, if high
            cih_L = R_L * a_prime #cash-in-hand in t=2 (cih_2) in z_1 and a_1, if low
            Ev = Evx = zero(a_prime)
            for s2 in eachindex(h.z_grid) # For each z_2
                v_H = lst_v[s2,1](cih_H) # This is value v_2(cih_2,z_2,A_2) when A_2 is high, z_1 = s, a_1 = i (and thus, cih_2=cih_H)
                v_L = lst_v[s2,2](cih_L) # This is value v_2(cih_2,z_2,A_2) when A_2 is low, z_1 = s, a_1 = i (and thus, cih_2=cih_L)
                x_H = lst_x[s2,1](cih_H) # This is x_2(cih_2,z_2,A_2) when A_2 is high, z_1 = s, a_1 = i (and thus, cih_2=cih_H)
                x_L = lst_x[s2,2](cih_L) # This is x_2(cih_2,z_2,A_2) when A_2 is low, z_1 = s, a_1 = i (and thus, cih_2=cih_L)
                Evx += P′[s2, s] * (probH * (v_H)^((ξ - γ)) * R_H * (x_H)^(-ξ) + (1-probH) * (v_L)^((ξ - γ)) * R_L * (x_L)^(-ξ)) 
                Ev += P′[s2, s] * (probH * u.risk(v_H) + (1-probH) * u.risk(v_L))
            end
            current_values.v[i, s] = inverse(u.temporal, (1 - β) * u.temporal(x) + β * u.temporal(inverse(u.risk, Ev)))
        end
    end
end