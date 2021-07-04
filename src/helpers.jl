
function my_sup_norm(a, b)
    s = zero(eltype(a))
    @turbo for i in indices((a, b))
        s = max(s, abs(a[i] - b[i]))
    end
    return s
end


# helper function to compute the ergodic distribution
function ergodic(P; tol = 10^(-12))
    all(sum(P, dims=2) .≈ 1.0) || error("Matrix is not valid. Rows don't sum to one")
    aa = ones(size(P, 1))
    while true
        bb = P' * aa
        norm(bb - aa) < tol && break
        aa = bb
    end
    return aa ./ length(aa)
end
