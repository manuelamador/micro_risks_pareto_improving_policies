###
# # AGGREGATE RISK TRANSFER FUNCTIONS


"""
    get_T_AR(; b, bprime, r, r_lf, k)

Transfers `T` assuming goverment BC holds with equality. The value of `b` and `bprime` are the values of the government debt in the current and next period, respectively.
`r_lf` is the return on capital in the laissez-faire transition path.
"""
get_T_AR(; b, bprime, r, r_lf, k) = bprime - (1 + r) * b - (r - r_lf) * k 


"""
    get_T_portfolio_AR(; b, bprime, r_b, r_k, r_lf, k)

Transfers `T` assuming goverment BC holds with equality. The value of `b` and `bprime` are the values of the government debt in the current and next period, respectively.
`r_b` is the (non-state contingent) interest rate on bonds, `r_k` is the (state-contingent) return on capital, and `r_lf` is the return on capital in the laissez-faire transition path. 
"""
get_T_portfolio_AR(; b, bprime, r_b, r_k, r_lf, k) = bprime - (1 + r_b) * b - (r_k - r_lf) * k



