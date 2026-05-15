using CairoMakie, LaTeXStrings
using LsqFit
using CSV, DataFrames
using Statistics


function log_log_plot(df)#telegate_error_rates, logical_error_rate, acceptance_ratio)#, power_exp, power_interc)

    # ---------------------- Pre-Processing ----------------------
    ps_all = sort(unique(Float64.(df.p)))
    p_range = range(minimum(ps_all), maximum(ps_all), length=300)
    p_bells = sort(unique(Float64.(df.p_bell)))
    # for 3 telegate noises: filter by a certain p_bell noise and add data to te plot
    p_bells_plot = [minimum(p_bells), p_bells[Int(floor(end/2))], maximum(p_bells)]
   # print(p_bells_plot)
    
    #logical_error_rates = Float64.(df.logical_error_rate)
    #df.diff_logical_phys= logical_error_rates .- local_error_rates  
    #acceptance_ratios = Float64.(df.acceptance_ratios)

    #p_bells = sort(unique(error_rates_bell))
    
    #p_range  = range(minimum(error_rates_bell), maximum(error_rates_bell), length = 300)
    

    # ── plot ──────────────────────────────────────────────────────────────────────
    fig = Figure(size = (640, 920), fontsize = 14)

    ax1  = Axis(fig[1, 1],
        xlabel             = "Physical noise rate",
        ylabel             = "Logical noise rate",
        title              = "Logical vs Physical Noise",
        xscale             = log10,
        yscale             = log10,
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

    ax2 = Axis(fig[2, 1],
        xlabel = "Physical noise rate",
        ylabel = "Acceptance ratio",
        xscale = log10,
        #xtickformat = xs -> power_of_10_label.(xs),
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

    #p_anchor = 1          # anchor point for the p and p² guide lines
    #y_anchor = 1          # where the guide lines pass through at p_anchor

    # p_L = p  (slope 1 in log-log)
    lines!(ax1, p_range, p_range,# y_anchor .* (p_range ./ p_anchor).^1,
        color     = :gray60,
        linestyle = :dash,
        linewidth = 1.5,
        label     = L"p_L \sim p",
    )

    # p_L = p²  (slope 2 in log-log)
    lines!(ax1, p_range, p_range.^2, #y_anchor .* (p_range ./ p_anchor).^2,
        color     = :gray40,
        linestyle = :dot,
        linewidth = 1.5,
        label     = L"p_L \sim p^2",
    )

    colors = Makie.wong_colors()


    for (idx,p_bell) in enumerate(p_bells_plot)

        #print(idx,p_bell)
        df_ = copy(df[df.p_bell .== p_bell, :])
        sort!(df_, :p)
        
        df_.diff_logical_phys = df_.logical_error_rate .- df_.p  
        pseudo_thresh_ind = findfirst(i -> df_.diff_logical_phys[i] * df_.diff_logical_phys[i+1] < 0, 1:length(df_.diff_logical_phys)-1)
        pseudo_thresh = isnothing(pseudo_thresh_ind) ? nothing : df_.p[pseudo_thresh_ind]

        if idx ==1
            # We only extract the scaling line for the lowest noise
            @info "Retrieving power law scaling for p_bell = $p_bell, assuming power law scaling p_log = a p_phys^b"
            a, b = power_law_fit(df_.p, df_.logical_error_rate)
            fit_label = L"\text{fit for p_{bell}}=%$(round(p_bell, digits=2)): p_L=%$(round(a, digits=2))\, p^{%$(power_of_10_label( round(b, digits=10)))}"
            # fitted line:  p_L = A * p^α
            lines!(ax1, p_range, a .* p_range .^ b,
            color     = colors[idx],
            linestyle = :dashdot,
            linewidth = 2,
            label     = fit_label
        )
        end
        data_label = L"\text{Logical noise for p_{bell}}=%$( power_of_10_label( round(p_bell, digits=10) ) )"
        acc_label  = L"\text{Acceptance ratio for p_{bell}}=%$( power_of_10_label( round(p_bell, digits=10)) )"
       
        # data
        scatterlines!(ax1, df_.p, df_.logical_error_rate,
            color      = colors[idx],
            markersize = 10,
            linewidth  = 2,
            label      = data_label
        )

        scatterlines!(ax2, df_.p, df_.acceptance_ratio,
            color = colors[idx], markersize = 10, linewidth = 2,
            label = acc_label
        )

        if pseudo_thresh !== nothing
            vlines!(ax1, [pseudo_thresh], color= :green, linestyle=:dash)
            vlines!(ax2, [pseudo_thresh], color= :green, linestyle=:dash)
        end
    end


    axislegend(ax1, position = :rb)
    axislegend(ax2, position = :rb)

    #axislegend(ax, position = :rb)
    #display(fig)
    save(joinpath(data_path,"..", "qec_threshold.png"), fig)

end

function two_d_plot(df)#telegate_error_rates, logical_error_rate, acceptance_ratio)#, power_exp, power_interc)

    # read from file
    
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
    Colorbar(fig[1, 2], hm; label = "log(ϵ_logical/ϵ_p)")

    save(joinpath(data_path,"..", "2d_heatmap.png"), fig)

end


function power_law_fit(physical_noise, logical_noise)


    m(t,p) = p[1] .* t .^p[2]# exp.(p[2] * t)
    p0 = [1.0, 2.0]
    fit = curve_fit(m, physical_noise, logical_noise, p0)
    a, b = fit.param
    σ_a, σ_b = stderror(fit)

    @info "a = $(round(a, digits=4)) ± $(round(σ_a, digits=4)), b = $(round(b, digits=3)) ± $(round(σ_b, digits=3))"
    return a, b
end

function power_of_10_label(val)
    #print(val)
    n = round(Int, log10(val))
    return L"10^{%$n}"
end

#code =    # Set code and architecture here
data_path = joinpath(@__DIR__, "..", "..", "data", "Steane/[4, 3]/simulation/dqc_sim_data.csv")
df = CSV.read(data_path, DataFrame)
log_log_plot(df)
two_d_plot(df)