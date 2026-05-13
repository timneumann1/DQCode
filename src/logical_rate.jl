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

function log_log_plot()#telegate_error_rates, logical_error_rate, acceptance_ratio)#, power_exp, power_interc)


    # read from file
    path = joinpath(@__DIR__, "results", "TrivariateBicycle", "[6, 6]", "error_rates.csv")
    df = CSV.read(path, DataFrame)
    print(df)
    telegate_error_rates = Float64.(df.telegate_error_rates)[3:60]#end]#50]
    logical_error_rate = Float64.(df.logical_error_rates)[3:60]#end]#50]
    phys_v_log = telegate_error_rates - logical_error_rate
    println(phys_v_log)
    pseudo_thresh_ind = findfirst(i -> phys_v_log[i] * phys_v_log[i+1] < 0, eachindex(phys_v_log)[1:end-1])
    pseudo_thresh = telegate_error_rates[pseudo_thresh_ind]

    acceptance_ratio = Float64.(df.acceptance_ratios)[3:60]#end]#50]
    println(telegate_error_rates)
    a, b = power_law_fit(telegate_error_rates, logical_error_rate)


    # ── tick formatter ────────────────────────────────────────────────────────────
    function power_of_10_label(val)
        n = round(Int, log10(val))
        return L"10^{%$n}"
    end

    # ── reference scaling lines, anchored at a mid-range point ───────────────────
    p_range  = range(minimum(telegate_error_rates), maximum(telegate_error_rates), length = 300)
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

    # fitted line:  p_L = A * p^α
    lines!(ax1, p_range, a .* p_range .^ b,
        color     = colors[2],
        linestyle = :dashdot,
        linewidth = 2,
        #label     = L"fit: $p_L = %.2f\, p^{%.2f}$" % (a,b),   # see note below
    )

    # data
    scatterlines!(ax1, telegate_error_rates, logical_error_rate,
        color      = colors[1],
        markersize = 10,
        linewidth  = 2,
        label      = "Logical noise",
    )
    scatterlines!(ax2, telegate_error_rates, acceptance_ratio,
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
    save(joinpath(@__DIR__, "results", "TrivariateBicycle", "[6, 6]", "qec_threshold.png"), fig)

end