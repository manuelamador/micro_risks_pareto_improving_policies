
# Displays a statistical summary of the economy.

function generate_tables(laissez_faire, constant_k_final, constant_k_path, golden_k_final, golden_k_path; io = stdout)

    table1 = zeros(5, 2)
    
    pdf_marginal = sum(laissez_faire.ws.pdf, dims=2)[:, 1]
    a_grid = laissez_faire.h.a_grid
    qu_10 = [quantile(a_grid, StatsBase.ProbabilityWeights(pdf_marginal), x) for x in 0.0:0.1:1.0]

    weight = similar(laissez_faire.ws.pdf)

    for (i, path) in enumerate([constant_k_path, golden_k_path]) 
        w_gains = 100 .* (path.v[1] ./ laissez_faire.ws.v .- 1)

        _min = minimum(w_gains[laissez_faire.ws.pdf .> 0])
        _mean = sum(w_gains .* laissez_faire.ws.pdf)
     
        w_gains_rs = reshape(w_gains, 1, :)[1, :]
        wgains_summary = Float64[]

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
    
        _poor = first(wgains_summary)
        _rich = last(wgains_summary)
        _middle = sum(wgains_summary[5:6])/2 

       table1[:, i] .= [_mean, _min, _poor, _middle, _rich]
    end 

    table2 = zeros(10, 3)
    for (i, e) in enumerate([laissez_faire, constant_k_final, golden_k_final])    
        final_pdf_marginal = sum(e.ws.pdf, dims=2)[:, 1]
        qu = [quantile(a_grid, StatsBase.ProbabilityWeights(final_pdf_marginal), x) for x in 0.0:0.20:1.0]
    
        total = sum(a_grid .* final_pdf_marginal)
        dist = [sum((a_grid .<= qi) .* (a_grid .* final_pdf_marginal)) * 100 / total for qi in qu] |> diff
        
        table2[:, i] .= [
            [e.b / y(e) * 100,
            e.r * 100,
            e.k / y(e),
            e.T / y(e) * 100, 
            100 * sum(e.ws.pdf[1, :])];   # constrained households
            dist ]            
        
    end 

    out = [ 
        pretty_table(io, table2, header = ["Laissez Faire", "Constant K", "Golden K"], 
            formatters = (v, i, j) -> i ∉ [2, 3, 4] ? round(Int, v) : round(v, digits = 1), 
            row_names = ["b/y (%)", "r (%)", "k/y", "T/y (%)", "Shared constrained hh", "Q1 Wealth Share", "Q2 Wealth Share", "Q3 Wealth Share" , "Q4 Wealth Share", "Q5 Wealth Share"], title = "Steady State Statistics"
        ),
    
        pretty_table(io, table1, header = ["Constant K", "Capital Expansion"], formatters = ft_round(1), row_names = ["mean", "min", "poor (≤10)", "middle (40-60)", "rich (>90)"], title = "Welfare Gains at announcement")
    ]

    (io == String) && return out 
    return nothing
end 