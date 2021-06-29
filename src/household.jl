
# Helper function, iterates the value function using a policy function.
function _iterate_policy!(arrays, h, r, w, transfer)
    @unpack β, P, z_grid, u, a_grid = h
    R = 1 + r
    v0, v1, pol = arrays 
    distance = 0.0
    Threads.@threads for s in axes(v0, 2)
        @inbounds for i in indices((v0, pol, v1), 1)
            j = pol[i, s]
            x = get_x(h, w * z_grid[s], R * a_grid[i] - a_grid[j] + transfer)
            x = max(x, 0.0)
            Ev = 0.0
            @turbo for sprime in axes(v0, 2)
                Ev += P[s, sprime] * inside(u, v0[j, sprime])
            end
            v = outside(u, x, inverse_inside(u, Ev), β)
            distance = max(distance, abs(v - v0[i, s]))
            v1[i, s] = v
        end
    end
    return (v0, v1, distance)
end 


# Inner loop of the value function iteration for the solution of the household's
# problem. 
# Modifies in place arrays = (v0, v1, pol). v0 is the initial value. v1 is the 
# new one. pol stores the optimal policy.
function optimize_one_period!(arrays, h, r, w, transfer)
    @unpack β, P, z_grid, grid_points, u, a_grid = h
    R = 1 + r
    v0, v1, pol1 = arrays 
    Threads.@threads for s in axes(v1, 2)
        # Doing the optimization.
        pol = 1  # starting policy
        vmax = -Inf
        @inbounds for i in axes(v1, 1)
            improvement = false
            j = pol
            cond = true
            while cond 
                x = get_x(h, w * z_grid[s], R * a_grid[i] - a_grid[j] + transfer)
                x = max(x, 0.0)
                Ev = 0.0
                @turbo for sprime in axes(v1, 2)
                    Ev += P[s, sprime] * inside(u, v0[j, sprime])
                end
                temp = outside(u, x, inverse_inside(u, Ev), β)
                if temp >= vmax
                    improvement = true
                    vmax = temp 
                    pol = j
                elseif improvement
                    # value is not increasing anymore -- by concavity, stop 
                    cond = false
                end
                if j >= grid_points 
                    cond = false
                end 
                j += 1 
            end
            v1[i, s] = vmax
            pol1[i, s] = pol
        end
    end
end 


# Iterate the value function once.  
# This is the same as the inner loop of solve_stationary_household (optimize_one_period!),
# except that it does not compute distances and allocates new arrays.  
function optimize_one_period(v0, h, r, w, transfer)
    v1 = similar(v0)
    pol1 = similar(v0, Int)
    optimize_one_period!((v0, v1, pol1), h, r, w, transfer)
    return v1, pol1
end 


# Solves the household's problem given constant r,w, and transfer by value
# function iteration. Uses preallocated arrays to store the computation.
#
# pol_iter_trigger is the trigger policy iteration cutoff. If th sup distance of
# the policy functions is less than or equal to this level (in grid_points), then
# a policy iteration cycle is triggered. Set to -1 to ignore.  
function solve_stationary_household!(
    arrays,  # pre-allocated arrays 
    h::Household, r, w; 
    transfer = 0.0,
    tol = _TOL_VALUE, 
    max_iter = 10_000,
    transfer_check = true, 
    pol_iter_trigger::Int = 2  
)
    if transfer_check
        is_transfer_feasible(h, w, transfer) ||  
            @warn "Transfer $transfer too negative. Feasibility will break. Proceeding anyway."
    end

    v0, v1, pol0, pol1 = arrays
    iter = 0
    while true
        iter += 1 
        optimize_one_period!((v0, v1, pol1), h, r, w, transfer)
        # geting the distances
        distance = zero(eltype(v0))
        distance_pol = zero(eltype(pol1))
        @turbo for i in indices((v0, v1, pol0, pol1))
            distance = max(distance, abs(v0[i] - v1[i]))
            distance_pol = max(distance_pol, abs(pol0[i] - pol1[i]))
        end
        (distance < tol) && break
        (iter > max_iter) && error("Maximum number of iterations reached")
        if distance_pol <= pol_iter_trigger # policy almost converged
            # iterate the value function without optimizing
            while true
                v0, v1, dis = _iterate_policy!((v1, v0, pol1), h, r, w, transfer)
                (dis < tol) && break
            end
        end 
        v1, v0 = v0, v1
        pol1, pol0 = pol0, pol1
    end
    return (v = v1, pol = pol1)
end 


# Solves the household's problem given constant r,w, and transfer by 
# value function iteration. 
# It allocates new arrays for the computation.
function solve_stationary_household(
    h::Household, r, w; 
    transfer = 0.0,
    transfer_check = true, 
    pol_iter_trigger::Int = 2, 
    tol = _TOL_VALUE, 
    max_iter = 10_000
)
    arrays = prealloc_and_initialize_ss_PE(h; r, w, transfer)
    return solve_stationary_household!(arrays, h, r, w; transfer, transfer_check, tol, max_iter, pol_iter_trigger)
end 


# Preallocates the arrays for the ss_PE computation
function prealloc_ss_PE(h)
    v_0 = Array{Float64, 2}(undef, (h.grid_points, h.n))
    v_1 = similar(v_0)
    pol_0 = similar(v_0, Int64) 
    pol_1 = similar(pol_0)
    return (v_0=v_0, v_1=v_1, pol_0 = pol_0, pol_1=pol_1)
end


# preallocate and initializes the array for the ss_PE computation
function prealloc_and_initialize_ss_PE(h; r, w, transfer)
    @unpack z_grid, a_grid, u, β = h

    arrays = prealloc_ss_PE(h)
    @unpack v_0 = arrays

    # Doing an educated guess for the initial value function
    @inbounds for s in axes(v_0, 2), i in axes(v_0, 1) 
        x = get_x(h, w * z_grid[s], r * a_grid[i]  + transfer)
        v_0[i, s] = ss_value(u, x, β)
    end
    return arrays
end 


# Iterates the houseohld policy function and modifies/returns pdf_1 with the 
# new distribution. 
function iterate_pdf!(
    pdf_1, # pre-allocated array
    pdf_0,
    h, 
    pol
)
    @unpack P = h 

    fill!(pdf_1, 0.0)
    @tullio  pdf_1[pol[i, s], sprime] += P[s, sprime] * pdf_0[i, s]
    # @inbounds for s in axes(pol, 2), i in axes(pol, 1), sprime in axes(pol, 2)
    #     pdf_1[pol[i, s], sprime] += P[s, sprime] * pdf_0[i, s]
    # end
    return pdf_1
end 


# Wrapper for iterate_pdf! function.
iterate_pdf(pdf_0, h, pol) = iterate_pdf!(similar(pdf_0), pdf_0, h, pol)


# Solve for the ergodic distribution given a stationary solution to the
# household problem.
function stationary_pdf!(
    arrays, # pre-allocated arrays
    h, sol;
    tol = _TOL_VALUE 
)
    @unpack pol = sol

    pdf_0, pdf_1 = arrays

    # Initializing the pdf
    mass = 0.0 
    @turbo for i in eachindex(pdf_0, sol.v)
        # putting mass only on states that are feasible
        pdf_0[i] = isfinite(sol.v[i]) ? one(eltype(pdf_0)) : zero(eltype(pdf_0))
        mass += pdf_0[i]
    end
    pdf_0 .= pdf_0 ./ mass

    while true
        iterate_pdf!(pdf_1, pdf_0, h, pol)
        distance = my_sup_norm(pdf_0, pdf_1)
        (distance < tol) && break
        pdf_0, pdf_1 = pdf_1, pdf_0
    end
    return pdf_1
end 


# Wrapper for the stationary_pdf function.
function stationary_pdf(h, sol; tol = _TOL_PDF)
    stationary_pdf!(prealloc_ss_pdf(h), h, sol, tol = tol)
end 


# Preallocation for the ergodic distribution calculations.
function prealloc_ss_pdf(h)
    pdf_0 = Array{Float64, 2}(undef, (h.grid_points, h.n))
    pdf_1 = similar(pdf_0)
    return (pdf_0, pdf_1)
end 


# Checks whether a given level of transfer is consistent with the feasbility 
# for the households. 
function is_transfer_feasible(h, w, transfer)
    return minimum_feasible_transfer(h, w) <= transfer
end 


function minimum_feasible_transfer(h, w)
    @unpack v, z_grid = h
    wi = w .* z_grid
    y = map(x -> labor_income_given_w(v, x), wi)
    vofn = map(x -> disutility_given_w(v, x), wi) 
    return -minimum(y .- vofn)
end 


# checks that a_max is not binding 
amax_binding(sol) = any(last(col) == length(col) for col in eachcol(sol.pol))


# Checks that a household solution is valid. The values are finite and the amax is not binding. 
function is_valid(sol; verbose=true)
    are_values_valid = all(isfinite, sol.v)
    amax_OK = !amax_binding(sol)
    verbose && !are_values_valid && @warn "Some values are NaN, Inf or -Inf"
    verbose && !amax_OK && @warn "amax is binding"
    return amax_OK && are_values_valid
end


# Given a policy function, compute the individual consumption matrix
function consumption(h, r, w, transfer, pol)
    @unpack z_grid,  a_grid = h
    R = 1 + r
    c = similar(pol, Float64)
    for s in axes(c, 2), i in axes(c, 1)
        j = pol[i, s]
        c[i, s] = get_c(h, w * z_grid[s], R * a_grid[i] - a_grid[j] + transfer)
    end
    return c
end

consumption(e, alloc) = consumption(get_h(e), alloc.r, alloc.w, alloc.transfer, alloc.pol)


# Given a distribution compute the aggregate allocation. 
function aggregate(alloc, pdf) 
    aggregate = 0.0
    @turbo for i in eachindex(pdf)
        aggregate += pdf[i] * alloc[i]
    end 
    return aggregate
end


# Given a distribution compute the aggregate asset supply. 
function asset_supply(h, pdf)
    @unpack a_grid = h 
    assets = 0.0
    @turbo for s in axes(pdf, 2), i in axes(pdf, 1) 
        assets += pdf[i, s] * a_grid[i]
    end
    return assets
end


# Given a wage for the hh, compute the total labor supply.
function labor_supply(h::Household, w)
    @unpack v, z_grid, Pss = h
    n = 0.0
    for (p, z) in zip(Pss, z_grid)
        n += p * z * n_of_w(v, w * z)
    end
    return n 
end

