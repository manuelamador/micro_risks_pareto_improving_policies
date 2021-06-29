
# Displays a statistical summary of the economy. 

function summary_statics(
    equilibrium; 
    laissez_faire = nothing, path = nothing
)
    e = get_e(equilibrium)
    h = get_h(equilibrium)
    a_grid = h.a_grid

    pdf = sum(equilibrium.pdf, dims=2)[:, 1]

    println("\nFISCAL")
    println("debt: ", equilibrium.b)
    println("debt over y (%): ", equilibrium.b / equilibrium.y * 100)
    println("transfer over y (%): ", 100 * equilibrium.transfer / equilibrium.y)
    println("interest rate: ", equilibrium.r)

    println("\nAGGREGATES")
    println("capital over y: ", equilibrium.k / equilibrium.y)
    
    println("\nHOUSEHOLDS")
        
    con = consumption(h, equilibrium.r, equilibrium.w, equilibrium.transfer, equilibrium.pol) 
    con = reshape(con, 1, :)[1, :]
    allpdf = reshape(equilibrium.pdf, 1, :)[1, :]
    std_con = 100 * StatsBase.std(log.(con), StatsBase.ProbabilityWeights(allpdf))
    println("standard deviation of log c (*100): ", std_con, "\n")
    println("mass of constrained households (%): ", 100 * sum(equilibrium.pdf[1, :]))
    
    println("mean wealth (over y): " , equilibrium.s / equilibrium.y)
    println("median wealth (over y): ", 
        quantile(a_grid , StatsBase.ProbabilityWeights(pdf), 0.5) / equilibrium.y) 
    
    qu = [
        quantile(a_grid, StatsBase.ProbabilityWeights(pdf), x) for 
            x in 0.0:0.20:1.0
    ]
    
    total = sum(a_grid .* pdf)
    dist = [sum((a_grid .<= qi) .* (a_grid .* pdf)) * 100 / total for qi in qu] |> diff
    println("share of wealth per asset quintiles:\n   $dist")

    if path !== nothing && laissez_faire !== nothing
        println("\nTRANSITION")

        w_gains = 100 .* (path.v[1] ./ laissez_faire.v .- 1)
        
        m1 = minimum(w_gains[laissez_faire.pdf .> 0])
        m2 = maximum(w_gains[laissez_faire.pdf .> 0])
        m3 = sum(w_gains .* laissez_faire.pdf)
        println("welfare gain in transition: min = $m1,  max = $m2,  mean = $m3")
        
        wgains_summary = Float64[]

        w_gains_rs = reshape(w_gains, 1, :)[1, :]

        weight = similar(laissez_faire.pdf)

        qu_10 = [
            quantile(a_grid, StatsBase.ProbabilityWeights(pdf), x) for 
            x in 0.0:0.1:1.0
        ]

        for (qi, qii) in zip(qu_10, qu_10[2:end])
            if qii == qu_10[end]
                qii_ = qii + 10 * eps()
            else 
                qii_ = qii
            end
            for i in axes(weight, 1), j in axes(weight, 2)
                weight[i, j] = (qi <= a_grid[i] < qii_) ? laissez_faire.pdf[i, j] : 0.0
            end
            weight_rs = reshape(weight, 1, :)[1, :]
            push!(wgains_summary, mean(w_gains_rs, ProbabilityWeights(weight_rs)))
        end

        println("mean welfare gain per asset decile:\n $wgains_summary")
    end
end 


function do_plots(tran, laissez_faire)
    t = get_t(laissez_faire)
    h = get_h(laissez_faire)
    n = laissez_faire.n
        
    ylist  = map(k -> get_y(t; k , n), tran.path.k)
    size = (400, 300)
    y0 = laissez_faire.y
        
    flist = []
    
    #################
    #  B 
    
    push!(flist, plot(100 .* tran.path.b ./ y0; size, 
            title = L"$b/y_{0}$", 
            ylabel = L"\%", ylims = (0, 100)))

    ##################
    # K 
    
    push!(flist, begin
            plot(tran.path.k ./ y0; size, 
                title = L"$k/y_{0}$" , 
                ylims = (0, 1.1 * maximum(tran.path.k ./ y0))
            )
        end
    )

    ##################
    # TRANSFERS
 
    push!(flist, 
        begin 
            plot(100 .* tran.path.transfer ./ y0; size, 
                title = "Transfers and bond revenue", 
                ylabel = L"\%", label = L"$T/y_0$",
                fill = 0, alpha = 0.2, 
                legend = (.5, .5))
            plot!(100 .* tran.path.transfer ./ y0; color = 1)
            plot!(
                [100 * (bprime - (1 + r) * b) / y0 for (b, bprime, r) in 
                        zip(tran.path.b, tran.path.b[2:end], tran.path.r, ylist)]; 
                size, 
                label = L"$\frac{b' - (1+r)b}{y_0}$")
        end
    )
    
    ######################
    #  R
    
    push!(flist, plot(tran.path.r .* 100; size, title = L"$r$", ylabel = L"\%"))

    ######################
    # CONSUMPTION 
    
    agg_c = [y + (1 - t.δ) * k - kprime for (y, k, kprime) in zip(ylist, tran.path.k, tran.path.k[2:end])]
    c0 = laissez_faire.y - t.δ * laissez_faire.k
    push!(flist, begin 
            plot( 100 .* (agg_c ./ c0  .- 1), fill = 0, alpha = 0.2, title = L"$c/c_0$", ylabel = L"\%")
            plot!(100 .* (agg_c ./ c0  .- 1), color = 1)
            ylims!(-5, 5)
        end
    )

    #########################
    # STD DEV 
    
    std_con = map(tran.path) do alloc  
        con = reshape(
                consumption(h, alloc.r, alloc.w, alloc.transfer, alloc.pol),
                1, :)[1, :]
        allpdf = reshape(alloc.pdf, 1, :)[1, :]
        100 * StatsBase.std(log.(con), StatsBase.ProbabilityWeights(allpdf))
    end 
    
    con_init = consumption(h, laissez_faire.r, laissez_faire.w, laissez_faire.transfer, laissez_faire.pol)
    con_init_ = reshape(con_init,1, :)[1, :]
    allpdf_init = reshape(laissez_faire.pdf, 1, :)[1, :]
    std_con_init = 100 * StatsBase.std(log.(con_init_), StatsBase.ProbabilityWeights(allpdf_init))

    push!(flist, 
        begin 
            plot(std_con .- std_con_init; size, 
                fill = 0, alpha = 0.2, 
                title = L"$\sigma(\log c) - \sigma(\log c_0)$", ylabel = L"\%")
            plot!(std_con .- std_con_init; color = 1)
        end
    )
    
    for f in flist
        xlims!(f, 0, 100)  # only the 100 first periods
    end
    
    return flist
end
