using CairoMakie#, LaTeXStrings
using LsqFit
using CSV, DataFrames


function power_law_fit(physical_noise, logical_noise)


    m(t,p) = p[1] .* t .^p[2]# exp.(p[2] * t)
    p0 = [1.0, 2.0]
    fit = curve_fit(m, physical_noise, logical_noise, p0)
    a, b = fit.param
    σ_a, σ_b = stderror(fit)

    println("Assuming power law scaling p_log = a p_phys^b")
    println("a = $(round(a, digits=4)) ± $(round(σ_a, digits=4))")
    println("b = $(round(b, digits=3)) ± $(round(σ_b, digits=3))")
    return a, b
end

function log_log_plot(data_path)#telegate_error_rates, logical_error_rate, acceptance_ratio)#, power_exp, power_interc)

    # read from file
    
    df = CSV.read(data_path, DataFrame)
    local_error_rates = Float64.(df.p)
    error_rates_bell = Float64.(df.p_bell)#[3:60]#end]#50]
    logical_error_rates = Float64.(df.logical_error_rate)#[3:60]#end]#50]
    phys_v_log = local_error_rates - logical_error_rates
    println(phys_v_log)
    pseudo_thresh_ind = findfirst(i -> phys_v_log[i] * phys_v_log[i+1] < 0, eachindex(phys_v_log)[1:end-1])
    pseudo_thresh_local_error = local_error_rates[pseudo_thresh_ind]

    acceptance_ratio = Float64.(df.acceptance_ratios)[3:60]#end]#50]
    println(error_rates_bell)


    # ── tick formatter ────────────────────────────────────────────────────────────
    function power_of_10_label(val)
        n = round(Int, log10(val))
        return L"10^{%$n}"
    end

    # ── reference scaling lines, anchored at a mid-range point ───────────────────
    p_range  = range(minimum(error_rates_bell), maximum(error_rates_bell), length = 300)
    p_anchor = 1          # anchor point for the p and p² guide lines
    y_anchor = 1          # where the guide lines pass through at p_anchor

    # ── plot ──────────────────────────────────────────────────────────────────────
    fig = Figure(size = (640, 920), fontsize = 14)

    ax1  = Axis(fig[1, 1],
        xlabel             = "Physical noise rate",
        ylabel             = "Logical noise rate",
        title              = "Logical vs Physical Noise",
        xscale             = log10,
        yscale             = log10,
        xtickformat        = xs -> power_of_10_label.(xs),
        ytickformat        = ys -> power_of_10_label.(ys),
        xgridvisible       = true,
        ygridvisible       = true,
        xminorgridvisible  = true,
        yminorgridvisible  = true,
        xminorticksvisible = true,
        yminorticksvisible = true,
        xminorticks        = IntervalsBetween(9),
        yminorticks        = IntervalsBetween(9),
    )

    ax2 = Axis(fig[2, 1],
        xlabel = "Physical noise rate",
        ylabel = "Acceptance ratio",
        xscale = log10,
        xtickformat = xs -> power_of_10_label.(xs),
        xgridvisible = true,
        ygridvisible = true,
        xminorgridvisible = true,
        yminorgridvisible = true,
        xminorticksvisible = true,
        yminorticksvisible = true,
        xminorticks = IntervalsBetween(9),
        yminorticks = IntervalsBetween(9),
    )

    linkxaxes!(ax1, ax2)

    colors = Makie.wong_colors()

    # p_L = p  (slope 1 in log-log)
    lines!(ax1, p_range, y_anchor .* (p_range ./ p_anchor).^1,
        color     = :gray60,
        linestyle = :dash,
        linewidth = 1.5,
        label     = L"p_L \sim p",
    )

    
    # p_L = p²  (slope 2 in log-log)
    lines!(ax1, p_range, y_anchor .* (p_range ./ p_anchor).^2,
        color     = :gray40,
        linestyle = :dot,
        linewidth = 1.5,
        label     = L"p_L \sim p^2",
    )



    # for 5 telegate noises: filter by a certain p_bell noise and add data to te plot

        a, b = power_law_fit(error_rates_bell, logical_error_rates)


        # fitted line:  p_L = A * p^α
        lines!(ax1, p_range, a .* p_range .^ b,
            color     = colors[2],
            linestyle = :dashdot,
            linewidth = 2,
            #label     = L"fit: $p_L = %.2f\, p^{%.2f}$" % (a,b),   # see note below
        )

        # data
        scatterlines!(ax1, error_rates_bell, logical_error_rates,
            color      = colors[1],
            markersize = 10,
            linewidth  = 2,
            label      = "Logical noise",
        )
        scatterlines!(ax2, error_rates_bell, acceptance_ratio,
            color = colors[3], markersize = 10, linewidth = 2,
            label = "Acceptance ratio")

        for ax in (ax1, ax2)
            vlines!(ax, [pseudo_thresh],
                color = :red,
                linestyle = :dash,
                linewidth = 2,
                label = ax === ax1 ? "Pseudothreshold" : nothing)
        end
    

    axislegend(ax1, position = :rb)
    axislegend(ax2, position = :rb)

    #axislegend(ax, position = :rb)
    #display(fig)
    save(joinpath(data_path,"..", "qec_threshold.png"), fig)

end

function two_d_plot(data_path)#telegate_error_rates, logical_error_rate, acceptance_ratio)#, power_exp, power_interc)

    # read from file
    
    df = CSV.read(data_path, DataFrame)
    #print(df)
    local_error_rates = Float64.(df.p)
    error_rates_bell = Float64.(df.p_bell)#[3:60]#end]#50]
    logical_error_rates = Float64.(df.logical_error_rate)#[3:60]#end]#50]
    df.ratio_logical_phys= logical_error_rates ./ local_error_rates  
    
    ps      = sort(unique(local_error_rates))
    p_bells = sort(unique(error_rates_bell))

    # build Z so that Z[i,j] corresponds to ps[i], p_bells[j]
    Z = fill(NaN, length(ps), length(p_bells))
    p_to_i  = Dict(p => i for (i, p) in enumerate(ps))
    pb_to_j = Dict(pb => j for (j, pb) in enumerate(p_bells))

    for r in eachrow(df)
        Z[p_to_i[r.p], pb_to_j[r.p_bell]] = r.ratio_logical_phys# r.logical_error_rate
    end

  
    # ── plot ──────────────────────────────────────────────────────────────────────
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1];
        xlabel = "p",
        ylabel = "p_bell",        
        # xscale = log10,
        # yscale = log10,
        title  = "Logical error rates",
        #xtickformat        = xs -> power_of_10_label.(xs),
        #ytickformat        = ys -> power_of_10_label.(ys),
        xgridvisible       = true,
        ygridvisible       = true,
        xminorgridvisible  = true,
        yminorgridvisible  = true,
        xminorticksvisible = true,
        yminorticksvisible = true,
        xminorticks        = IntervalsBetween(9),
        yminorticks        = IntervalsBetween(9),
    )

    Zlog = log10.(Z) # for plotting, we use the log of the ratio
    maxabs = maximum(abs, Zlog)
    hm = heatmap!(ax, ps, p_bells, Zlog; colormap = cgrad([:blue, :white, :red]), colorrange = (-maxabs, maxabs))

    # maxabs = maximum(abs, Z)
    # hm = heatmap!(ax, ps, p_bells, Z; colormap = cgrad([:blue, :white, :red]), colorrange = (0, maxabs))
    Colorbar(fig[1, 2], hm; label = "ϵ_logical/ϵ_p")

    save(joinpath(data_path,"..", "2d_heatmap.png"), fig)

end

function power_of_10_label(val)
    n = round(Int, log10(val))
    return L"10^{%$n}"
end

#code =    # Set code and architecture here
data_path = joinpath(@__DIR__, "..", "..", "data", "Steane/[4, 3]/simulation/dqc_sim_data.csv")
#log_log_plot(data_path)
two_d_plot(data_path)