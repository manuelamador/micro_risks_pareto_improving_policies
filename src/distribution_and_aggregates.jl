
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
    mass = zero(eltype(pdf_0))
    for i in eachindex(pdf_0, sol.v)
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
stationary_pdf(h, sol; tol = _TOL_PDF) = stationary_pdf!(prealloc_pdfs(h), h, sol, tol = tol)


# Preallocation for the ergodic distribution calculations.
function prealloc_pdfs(h)
    pdf_0 = Array{Float64, 2}(undef, (h.grid_points, h.n))
    pdf_1 = similar(pdf_0)
    return (pdf_0, pdf_1)
end


# Given a distribution and allocation compute the aggregate allocation.
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
        n += p * z * labor(v, w * z)
    end
    return n
end
