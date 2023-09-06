


################################################
# Technology Methods
###############################################


output(t::CobbDouglasTechnology; k, n) = t.A * k^t.α * n^(1- t.α)

mpk(t::CobbDouglasTechnology; k, n) = t.α * t.A * (n / k)^(1 - t.α)

mpl_from_mpk(t::CobbDouglasTechnology; mpk) = (1 - t.α) * (t.A / t.μ)^(1/(1-t.α)) * (t.α / mpk)^(t.α / (1 - t.α))

k_from_mpk_n(t::CobbDouglasTechnology; mpk, n) = (t.α * t.A / (t.μ * mpk))^(1/(1 - t.α)) * n 

golden_rule_k(t::CobbDouglasTechnology; n) = (t.α * t.A / t.δ) ^ (1 / (1 - t.α)) * n

y(e) = output(e.t; e.k, e.n)


################################################
# Fiscal policy methods
###############################################

function get_T(t; b, bprime, r, k, k0, r0, n0) 
    # Use Lemma 1 to obtain T assuming goverment BC holds with equality
    return output(t; k, n = n0) - output(t; k = k0, n = n0) + (r0 + t.δ) * k0 - (r + t.δ) * k - (1 + r) * b + bprime 
end 
