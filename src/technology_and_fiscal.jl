################################################
# Technology Methods
###############################################

output(t; k, n) = t.A * k^t.α * n^(1- t.α)

mpk(t; k, n) = t.α * t.A * (n / k)^(1 - t.α)

mpl_from_mpk(t; mpk) = (1 - t.α) * (t.A / t.μ)^(1/(1-t.α)) * (t.α / mpk)^(t.α / (1 - t.α))

k_from_mpk_n(t; mpk, n) = (t.α * t.A / (t.μ * mpk))^(1/(1 - t.α)) * n 

golden_rule_k(t; n) = (t.α * t.A / t.δ) ^ (1 / (1 - t.α)) * n

y(e) = output(e.t; e.k, e.n)

################################################
# Fiscal policy methods
###############################################

"""
    get_T(t; b, bprime, r, k, r0, k0, n0)

Use Lemma 1 to obtain transfers `T` assuming goverment BC holds with equality. `(r0, k0, n0)` are the initial values of `(r, k)`. 
The value of `b` and `bprime` are the values of the government debt in the current and next period, respectively.
The wage `w` is assumed to be fixed at the initial equilibrium value, so is the labor supply, given no wealth effects.
"""
get_T(t; b, bprime, r, k, r0, k0, n0) = output(t; k, n = n0) - output(t; k = k0, n = n0) + (r0 + t.δ) * k0 - (r + t.δ) * k - (1 + r) * b + bprime 

