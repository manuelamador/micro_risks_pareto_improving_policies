

# helper function to compute the ergodic distribution from transition matrix P
function ergodic(P)
    for i in axes(P, 1)
        !(sum(@view P[i, :]) ≈ 1.0) && error("Matrix is not valid. Rows don't sum to one")
    end 
    # from the calibration function
    invariant_dist = real(inv(eigvecs(P))[size(P, 1), :])
    invariant_dist .= invariant_dist ./ sum(invariant_dist)
    return invariant_dist
end


# Linear interpolation function
function interp1D(x::Real, xvec, yvec)
    i = searchsortedfirst(xvec, x)
    ind = clamp(i, firstindex(xvec) + 1, lastindex(xvec))
    @inbounds val = yvec[ind - 1] + (yvec[ind] - yvec[ind - 1]) / (xvec[ind] - xvec[ind - 1]) * (x - xvec[ind - 1])
    return val
end 


chebyshev(a, b) = @tullio (max) x = abs(a[i] - b[i])


logrange(start, stop, shifter, length) = (10^y - shifter for y in range(log10(start + shifter), log10(stop + shifter); length))


function grid(; start = 0.0, stop = 15.0, length = 500, scale = :log, shifter = 1 - start)  
    if scale == :linear 
        grid = collect(range(start, stop, length = length))
    elseif scale == :log 
        grid = collect(logrange(start, stop, shifter, length))
    else 
        @error "Scale neither linear nor log."
    end 
    return grid
end 
