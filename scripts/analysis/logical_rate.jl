using CairoMakie, LaTeXStrings
using CairoMakie: Axis
using LsqFit
using CSV, DataFrames
using Statistics


# ---------------------------------------------------------------
# ------------------ 2D Heatmap ---------------------------------
# ---------------------------------------------------------------


function log_log_plot(df, data_path; monolithic = false)#telegate_error_rates, logical_error_rate, acceptance_ratio)#, power_exp, power_interc)

    # ---------------------- Pre-Processing ----------------------
    df.logical_error_rate = ifelse.(df.logical_error_rate .== 0.0, 1 / df.n_samples[1], df.logical_error_rate)

    ps_all = sort(unique(Float64.(df.p)))
    p_range = range(minimum(ps_all), maximum(ps_all), length=300)
    p_bells = sort(unique(Float64.(df.p_bell)))

    if monolithic
        p_bells_plot = [minimum(p_bells)]
    else
        # for 3 telegate noises: filter by a certain p_bell noise and add data to the plot

        p_bells_plot = [minimum(p_bells), p_bells[Int(floor(end/2))], maximum(p_bells)]
    end

   
    # ---------------------- Plot ----------------------
    fig = Figure(size = (640, 920), fontsize = 14)
    
    ys = filter(x -> x > 0, Float64.(df.logical_error_rate))
    minexp = floor(Int, log10(minimum(ys)))
    maxexp = ceil(Int,  log10(maximum(ys)))
    @info "min: $minexp"
    ytick_vals = 10.0 .^ (minexp:maxexp)
    ytick_labels = power_of_10_label.(ytick_vals)

    ax1  = Axis(fig[1, 1], ylabel = L"\text{Logical } |0 \rangle_L \text{ initialisation error}",  title = L"\text{Logical vs. Physical Initialisation Error}",
        xscale  = log10, yscale  = log10,
        yticks = (ytick_vals, ytick_labels), 
        ylabelsize = 22, titlesize = 28,
        xgridvisible= true, ygridvisible = true, xminorgridvisible  = true, yminorgridvisible  = true, xminorticksvisible = true, yminorticksvisible = true,
        xminorticks = IntervalsBetween(9), yminorticks = IntervalsBetween(9))

    ax2 = Axis(fig[2, 1],  xlabel = L"\text{Physical error rate } p", ylabel = L"\text{Acceptance ratio}",
        xscale = log10, xlabelsize = 22, ylabelsize = 22,
        xgridvisible = true, ygridvisible = true,  xminorgridvisible = true,  yminorgridvisible = true,   xminorticksvisible = true, yminorticksvisible = true,
        xminorticks = IntervalsBetween(9), yminorticks = IntervalsBetween(9))
    
    ylims!(ax2, min( minimum(df.acceptance_ratio), 0.75), 1.0)
    linkxaxes!(ax1, ax2)

    # p_L = p  (slope 1 in log-log)
    lines!(ax1, p_range, p_range, color = :gray60,  linestyle = :dash, linewidth = 1.5,  label = L"p_L \sim p")

    # p_L = p²  (slope 2 in log-log)
    lines!(ax1, p_range, p_range.^2, color = :gray40, linestyle = :dot, linewidth = 1.5, label = L"p_L \sim p^2")

    colors = Makie.wong_colors()


    for (idx,p_bell) in enumerate(p_bells_plot)

        df_bell = copy(df[(df.p_bell .== p_bell), :])
        @info minimum(df_bell.logical_error_rate)
        df_ = copy(df_bell[(df_bell.logical_error_rate.>0.0), :])
        @info "Deleted $(nrow(df_bell)-nrow(df_)) row(s) because of 0.0 logical error rate"
        @info minimum(df_.logical_error_rate)
        @info minimum(df_.acceptance_ratio)
        sort!(df_, :p)
        df_.diff_logical_phys = df_.logical_error_rate .- df_.p  
        df_.std_error = sqrt.(1/5e5 .* (df_.logical_error_rate) .* (1 .- df_.logical_error_rate))
        # assumes monotonous scaling
        pseudo_thresh_ind = findfirst(i -> df_.diff_logical_phys[i] * df_.diff_logical_phys[i+1] < 0, 1:length(df_.diff_logical_phys)-1)
        pseudo_thresh = isnothing(pseudo_thresh_ind) ? nothing : df_.p[pseudo_thresh_ind]
        
        if idx ==1
            a, b = power_law_fit(df_.p, df_.logical_error_rate)
            if monolithic
                fit_label = L"\text{Power law fit: p_L=%$(round(a, digits=2))\, p^{%$( round(b, digits=2))} }"
            else
                fit_label = L"\text{PL fit for p_{bell}}=%$(power_of_10_label( round(p_bell, digits=10))): p_L=%$(round(a, digits=2))\, p^{%$( round(b, digits=2))}"
            end
            # fitted line:  p_L = A * p^α
            lines!(ax1, p_range, a .* p_range .^ b, color = colors[idx], linestyle = :dashdot, linewidth = 2, label = fit_label)
        end

        if monolithic
            data_label = L"\text{Logical initialisation noise}"
            acc_label  = L"\text{Acceptance ratio}"
        else
            data_label = L"\text{Logical initialisation noise for p_{bell}}=%$( power_of_10_label( round(p_bell, digits=10) ) )"
            acc_label  = L"\text{Acceptance ratio for p_{bell}}=%$( power_of_10_label( round(p_bell, digits=10)) )"
        end
  
        # data
        #errorbars!(ax1, df_.p, df_.logical_error_rate, df_.std_error; color =colors[idx], label=L"\sigma^2(LER)")
        scatterlines!(ax1, df_.p, df_.logical_error_rate, color = colors[idx],markersize = 10, linewidth  = 2, label = data_label)

        scatterlines!(ax2, df_.p, df_.acceptance_ratio, color = colors[idx], markersize = 10, linewidth = 2,  label = acc_label )

        if pseudo_thresh !== nothing && idx ==1
            vlines!(ax1, [pseudo_thresh], color= colors[idx], linestyle=:dash)
            vlines!(ax2, [pseudo_thresh], color= colors[idx], linestyle=:dash)
        end
    end

    axislegend(ax1, position = :rb)
    axislegend(ax2, position = :rb)
    save(joinpath(data_path,"..", "qec_threshold.png"), fig)

end

# ---------------------------------------------------------------
# ------------------ Scaling & Threshold ------------------------
# ---------------------------------------------------------------



function two_d_plot(df, data_path)

    # ---------------------- Pre-Processing ----------------------

    local_error_rates = Float64.(df.p)
    error_rates_bell = Float64.(df.p_bell)
    logical_error_rates = Float64.(df.logical_error_rate)
    df.ratio_logical_phys = logical_error_rates ./ local_error_rates  
    ps      = sort(unique(local_error_rates))
    p_bells = sort(unique(error_rates_bell))

    # build Z so that Z[i,j] corresponds to ps[i], p_bells[j]
    ratio = fill(NaN, length(ps), length(p_bells))
    LER = fill(NaN, length(ps), length(p_bells))
    p_to_i  = Dict(p => i for (i, p) in enumerate(ps))
    pb_to_j = Dict(pb => j for (j, pb) in enumerate(p_bells))

    for r in eachrow(df)
        ratio[p_to_i[r.p], pb_to_j[r.p_bell]] = r.ratio_logical_phys
        LER[p_to_i[r.p], pb_to_j[r.p_bell]] = r.logical_error_rate
    end

  
    # ---------------------- Plot ----------------------
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1]; xlabel = L"\text{Physical error rate } p",  ylabel = L"p_{Bell}", title  = L"\text{Logical } |0 \rangle_L \text{ per physical initialisation error and p_{Bell}}",   
        xscale = log10, yscale = log10, 
        xgridvisible = true, ygridvisible  = true, xlabelsize = 22, ylabelsize = 22, titlesize = 26,
        xminorgridvisible  = true, yminorgridvisible  = true, xminorticksvisible = true, yminorticksvisible = true, 
        xminorticks = IntervalsBetween(9),  yminorticks = IntervalsBetween(9),
    )

    log_ratio = map(z -> (z > 0.0) ? log10(z) : NaN, ratio) # for plotting, we use the log of the ratio
    vals = filter(x -> !isnan(x), vec(log_ratio))
    maxabs = isempty(vals) ? 1.0 : maximum(abs, vals)
    
   
    hm = heatmap!(ax, ps, p_bells, log_ratio; colormap = cgrad([:blue, :white, :red]), colorrange = (-maxabs, maxabs))
    Colorbar(fig[1, 2], hm; label = L"\log_{10}(\text{LER}/p)", labelsize = 22)
   
    save(joinpath(data_path,"..", "2d_heatmap_ratio.png"), fig)

    # ---------------------- Plot ----------------------
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1]; xlabel = L"\text{Physical error rate } p",  ylabel = L"p_{Bell}", title  = L"\text{Logical } |0 \rangle_L \text{ per physical initialisation error and p_{Bell}}",   
        xscale = log10, yscale = log10,
        xgridvisible = true, ygridvisible  = true, xlabelsize = 22, ylabelsize = 22,titlesize = 26,
        xminorgridvisible  = true, yminorgridvisible  = true, xminorticksvisible = true, yminorticksvisible = true,
        xminorticks = IntervalsBetween(9),  yminorticks = IntervalsBetween(9),
    )

    log_LER = map(z -> (z > 0) ? log10(z) : NaN, LER) # for plotting, we use the log of the ratio
    vals = filter(x -> !isnan(x), vec(log_LER))
    maxabs = isempty(vals) ? 1.0 : maximum(abs, vals)
    hm = heatmap!(ax, ps, p_bells, log_LER; colormap = cgrad(:viridis, rev=true), colorrange = (-maxabs, 0))
    Colorbar(fig[1, 2], hm; label = L"\log_{10}(\text{LER})", labelsize = 22)
   
    save(joinpath(data_path,"..", "2d_heatmap.png"), fig)

end

function power_law_fit(physical_noise, logical_noise)
    m(t,p) = p[1] .* t .^p[2]# exp.(p[2] * t)
    p0 = [1.0, 2.0]
    fit = curve_fit(m, physical_noise, logical_noise, p0)
    a, b = fit.param
    σ_a, σ_b = stderror(fit)
    @info "Retrieving power law scaling, assuming p_log = a p_phys^b -> a = $(round(a, digits=4)) ± $(round(σ_a, digits=4)), b = $(round(b, digits=3)) ± $(round(σ_b, digits=3))"
    return a, b
end

function power_of_10_label(val)
    #print(val)
    val > 0 || return L""
    n = round(Int, log10(val))
    return L"10^{%$n}"
end


code = "Steane"
qpu_sizes = "[4, 3]"
monolithic = false

data_path = joinpath(@__DIR__, "..", "..", "data", "$code/$qpu_sizes", "simulation_FT/dqc_sim_data.csv") # stores to simulation_FT/ or simulation_non_FT/ folder, depending on the indicated path
df = CSV.read(data_path, DataFrame)

if monolithic
    log_log_plot(df, data_path, monolithic = monolithic)
else
    log_log_plot(df, data_path)
    two_d_plot(df,data_path)
end


