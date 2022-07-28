
# Displays a statistical summary of the economy.

function summary_statics(
    e;
    iobuffer = stdout,
    laissez_faire = nothing, path = nothing
)
    h = e.h
    a_grid = h.a_grid

    pdf = sum(e.ws.pdf, dims=2)[:, 1]

    println(iobuffer, "\nFISCAL")
    println(iobuffer, "debt: ", e.b)
    println(iobuffer, "debt over y (%): ", e.b / y(e) * 100)
    println(iobuffer, "transfer over y (%): ", 100 * e.T / y(e))
    println(iobuffer, "interest rate: ", e.r)

    println(iobuffer, "\nAGGREGATES")
    println(iobuffer, "capital over y: ", e.k / y(e))

    println(iobuffer, "\nHOUSEHOLDS")

    con = consumption_alloc(h; R = 1 + e.r, e.w, e.T, e.ws.a_pol)
    con = reshape(con, 1, :)[1, :]
    allpdf = reshape(e.ws.pdf, 1, :)[1, :]
    std_con = 100 * StatsBase.std(log.(con), StatsBase.ProbabilityWeights(allpdf))
    println(iobuffer, "standard deviation of log c (*100): ", std_con, "\n")
    println(iobuffer, "mass of constrained households (%): ", 100 * sum(e.ws.pdf[1, :]))

    s = asset_supply(a_grid, e.ws.pdf)
    println(iobuffer, "mean wealth (over y): " , s / y(e))
    println(iobuffer, "median wealth (over y): ",
        quantile(a_grid , StatsBase.ProbabilityWeights(pdf), 0.5) / y(e))

    qu = [quantile(a_grid, StatsBase.ProbabilityWeights(pdf), x) for x in 0.0:0.20:1.0]

    total = sum(a_grid .* pdf)
    dist = [sum((a_grid .<= qi) .* (a_grid .* pdf)) * 100 / total for qi in qu] |> diff
    println(iobuffer, "share of wealth per asset quintiles:\n   $dist")

    if path !== nothing && laissez_faire !== nothing
        println(iobuffer, "\nTRANSITION")

        w_gains = 100 .* (path.v[1] ./ laissez_faire.ws.v .- 1)

        m1 = minimum(w_gains[laissez_faire.ws.pdf .> 0])
        m2 = maximum(w_gains[laissez_faire.ws.pdf .> 0])
        m3 = sum(w_gains .* laissez_faire.ws.pdf)
        println(iobuffer, "welfare gain in transition: min = $m1,  max = $m2,  mean = $m3")

        wgains_summary = Float64[]

        w_gains_rs = reshape(w_gains, 1, :)[1, :]

        weight = similar(laissez_faire.ws.pdf)

        qu_10 = [quantile(a_grid, StatsBase.ProbabilityWeights(pdf), x) for x in 0.0:0.1:1.0]

        for (qi, qii) in zip(qu_10, qu_10[2:end])
            if qii == qu_10[end]
                qii_ = qii + 10 * eps()
            else
                qii_ = qii
            end
            for i in axes(weight, 1), j in axes(weight, 2)
                weight[i, j] = (qi <= a_grid[i] < qii_) ? laissez_faire.ws.pdf[i, j] : 0.0
            end
            weight_rs = reshape(weight, 1, :)[1, :]
            push!(wgains_summary, mean(w_gains_rs, ProbabilityWeights(weight_rs)))
        end

        println(iobuffer, "mean welfare gain per asset decile:\n $wgains_summary")
    end
end

