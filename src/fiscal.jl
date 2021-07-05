
### Get tax rates given allocation and laissez_faire
function taxes(alloc, laissez_faire)
    t = get_t(laissez_faire)
    μ = get_μ(t)

    k1, n1 = alloc.k, alloc.n
    k0, n0 = laissez_faire.k, laissez_faire.n
    y0 = output(t, k = k0, n = n0)
    y1 = output(t, k = k1, n = n1)

    τn = mpn_from_factors(t, k = k1, n = n1) / μ / laissez_faire.w - 1
    τk = mpk_from_factors(t, k = k1, n = n1) / μ / (alloc.r + t.δ) - 1
    τπ = 1 - y0 / y1

    return LinearIncomeTaxes(τn, τk, τπ)
end
