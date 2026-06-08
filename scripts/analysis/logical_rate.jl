# logical_rate.jl

using CairoMakie, LaTeXStrings
using CairoMakie: Axis
using LsqFit
using CSV, DataFrames
using Statistics

# ---------------------------------------------------------------
# ------------------ Scaling & Threshold ------------------------
# ---------------------------------------------------------------

function log_log_plot(df, data_path; monolithic = false)
    @info "Creating log-log plot..."
    # ---------------------- Pre-Processing ----------------------
    df = copy(df)
    ps_all = sort(unique(Float64.(df.p)))
    p_range = range(minimum(ps_all), maximum(ps_all), length=300)
    p_bells = sort(unique(Float64.(df.p_bell)))
    if monolithic
        p_bells_plot = [minimum(p_bells)]
    else
        # choose three telegate noises for log-log plot
        p_bells_plot = [minimum(p_bells), p_bells[Int(floor(end/2))], maximum(p_bells)]
    end
    # ---------------------- Plot ----------------------
    fig = Figure(size = (900, 1120), fontsize = 14)
    log_err = filter(x -> x > 0, Float64.(df.logical_error_rate))
    minexp = floor(Int, log10(minimum(log_err)))
    maxexp = ceil(Int,  log10(maximum(log_err)))
    minexp_x = floor(Int, log10(minimum(ps_all)))
    maxexp_x = ceil(Int,  log10(maximum(ps_all)))
    xtick_vals   = 10.0 .^ (minexp_x:maxexp_x)
    xtick_labels = power_of_10_label.(xtick_vals,4)
    ytick_vals = 10.0 .^ (minexp:maxexp)
    ytick_labels = power_of_10_label.(ytick_vals,4)
    _code_architecture_label = df.code_architecture_label[1]
    ax1  = Axis(fig[1, 1],
        ylabel = L"\text{Logical } |0\rangle_L \text{ initialisation error}",  
        title = L"\text{Logical vs. Physical Initialisation Error}",
        subtitle = latexstring(_code_architecture_label * "; " *L"p_{\text{two-qubit gate}}=p_{\text{meas}}=p_{\text{init}}=p;"*" "*L"p_{\text{single-qubit gate}}="*"$(df.p_single_ratio[1])"*L"p;"*" "*L"p_{\text{idle}}="*"$(df.p_idle_ratio[1])"*L"p;"*" "*L"\text{depth_{telegate}}="*"$(df.telegate_idle_depth[1])"),
        titlegap = 12,
        subtitlegap = 10,
        xscale  = log10, 
        yscale  = log10,
        xticks = (xtick_vals, xtick_labels),
        xticklabelsize = 20, 
        yticks = (ytick_vals, ytick_labels), 
        yticklabelsize = 20, 
        ylabelsize = 22, 
        titlesize = 28,
        subtitlesize = 16,
        xgridvisible= true, 
        ygridvisible = true, 
        xminorgridvisible  = true, 
        yminorgridvisible  = true, 
        xminorticksvisible = true, 
        yminorticksvisible = true,
        xminorticks = IntervalsBetween(10), 
        yminorticks = IntervalsBetween(10)
    )
    ax2 = Axis(fig[2, 1],  
        xlabel = L"\text{Physical initialisation error rate } p", 
        ylabel = L"\text{Acceptance rate}",
        xscale = log10, 
        xlabelsize = 22, 
        ylabelsize = 22,
        xgridvisible = true, 
        ygridvisible = true,  
        xminorgridvisible = true,  
        yminorgridvisible = true,   
        xminorticksvisible = true, 
        yminorticksvisible = true,
        xminorticks = IntervalsBetween(10), 
        yminorticks = IntervalsBetween(10))
    ylims!(ax2, min( minimum(df.acceptance_ratio), 0.75), 1.0)
    linkxaxes!(ax1, ax2)
    # p_L = p  (slope 1 in log-log)
    lines!(ax1, p_range, p_range, color = :gray60,  linestyle = :dash, linewidth = 1.5,  label = L"p_L \sim p")
    # p_L = p²  (slope 2 in log-log)
    lines!(ax1, p_range, p_range.^2, color = :gray40, linestyle = :dot, linewidth = 1.5, label = L"p_L \sim p^2")
    colors = [RGBf(255/255, 177/255,  42/255), RGBf( 0/255, 154/255, 207/255), RGBf(96/255, 184/255,  72/255) ]
    for (idx,p_bell) in enumerate(p_bells_plot)
        df_bell = copy(df[(df.p_bell .== p_bell), :])
        df_ = copy(df_bell[(df_bell.logical_error_rate.>0.0), :])
        @info "Deleted $(nrow(df_bell)-nrow(df_)) row(s) because of 0.0 logical error rate"
        sort!(df_, :p)
        df_.diff_logical_phys = df_.logical_error_rate .- df_.p  
        df_.var = 1 ./ (df_.n_samples-df_.discarded_runs) .* (df_.logical_error_rate) .* (1 .- df_.logical_error_rate)
        df_.std_error = sqrt.(max.(df_.var, 0.0))
        # assumes monotonous scaling
        pseudo_thresh_ind = findfirst(i -> df_.diff_logical_phys[i] * df_.diff_logical_phys[i+1] < 0, 1:length(df_.diff_logical_phys)-1)
        pseudo_thresh = isnothing(pseudo_thresh_ind) ? nothing : df_.p[pseudo_thresh_ind]
        if idx ==1
            a, b = power_law_fit(df_.p, df_.logical_error_rate) # only determine scaling exponent for the lowest bell pair noise
            if monolithic # bell pair noise is irrelevant in this case
                fit_label = L"\text{Power law fit: p_L=%$(round(a, digits=2))\, p^{%$( round(b, digits=2))} }"
            else
                fit_label = latexstring(L"\text{PL fit for }p_{\text{Bell}}="*"=$(power_of_10_label(p_bell,4))"*L": p_L="*"$(round(a, digits=2))p^{$( round(b, digits=2))}")
            end
            # fitted line:  p_L = a*p^b
            lines!(ax1, p_range, a .* p_range .^ b, color = colors[idx], linestyle = :dashdot, linewidth = 2, label = fit_label)
        end
        if monolithic
            data_label = L"\text{Logical initialisation error}"
            acc_label  = L"\text{Acceptance rate}"
        else
            data_label = latexstring(L"\text{Logical initialisation error for } p_{\text{Bell}}="*"$(power_of_10_label(p_bell,4))")
            acc_label  = latexstring(L"\text{Acceptance rate for } p_{\text{Bell}}="*"$(power_of_10_label(p_bell,4))")
        end
        errorbars!(ax1, df_.p, df_.logical_error_rate, df_.std_error; color=colors[idx], label= idx == 1 ? L"\sigma(LER)" : nothing)
        scatterlines!(ax1, df_.p, df_.logical_error_rate, color = colors[idx], markersize = 8, linewidth  = 2, label = data_label)
        scatterlines!(ax2, df_.p, df_.acceptance_ratio, color = colors[idx], markersize = 8, linewidth = 2,  label = acc_label )
        if pseudo_thresh !== nothing && idx ==1
            vlines!(ax1, [pseudo_thresh], color= colors[idx], linestyle=:dash)
            vlines!(ax2, [pseudo_thresh], color= colors[idx], linestyle=:dash)
        end
    end
    axislegend(ax1, position = :rb)
    axislegend(ax2, position = :rb)
    save(joinpath(data_path, "../..", "qec_threshold.png"), fig)
    @info "Saved qec_threshold.png"
end


# ---------------------------------------------------------------
# ------------------ 2D Heatmap ---------------------------------
# ---------------------------------------------------------------


function two_d_plot(df, data_path)
    @info "Creating 2d heatmap..."
    # ---------------------- Pre-Processing ----------------------
    df = copy(df)
    local_error_rates = Float64.(df.p)
    error_rates_bell = Float64.(df.p_bell)
    df.ratio_logical_phys = df.logical_error_rate ./ local_error_rates  
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
    _code_architecture_label = df.code_architecture_label[1]
    # ---------------------- Plot 2D Heatmap ----------------------
    fig = Figure(size = (800, 600), fontsize = 14)
    ax = Axis(fig[1, 1]; 
            xlabel = L"\text{Physical initialisation error rate } p",  
            ylabel = L"\text{Bell pair initialisation error rate }p_{\text{Bell}}",  
            title = L"\text{Logical } |0\rangle_L\text{ Initialisation Error Rate (LIER) vs. } (p,p_{\text{Bell}})",
            subtitle  =  latexstring(_code_architecture_label * "; " *L"p_{\text{two-qubit gate}}=p_{\text{meas}}=p_{\text{init}}=p;"*" "*L"p_{\text{single-qubit gate}}="*"$(df.p_single_ratio[1])"*L"p;"*" "*L"p_{\text{idle}}="*"$(df.p_idle_ratio[1])"*L"p;"*" "*L"\text{depth_{telegate}}="*"$(df.telegate_idle_depth[1])"),
            titlegap = 12,
            subtitlegap = 10,
            xscale = log10, 
            yscale = log10, 
            xgridvisible = true, 
            ygridvisible  = true, 
            xlabelsize = 22, 
            ylabelsize = 22,
            xticklabelsize = 20,  
            yticklabelsize = 20, 
            titlesize = 28,
            subtitlesize = 16,
            xminorgridvisible  = true, 
            yminorgridvisible  = true, 
            xminorticksvisible = true, 
            yminorticksvisible = true, 
            xminorticks = IntervalsBetween(10),  
            yminorticks = IntervalsBetween(10),
    )
    log_ratio = map(z -> (z > 0.0) ? log10(z) : NaN, ratio) # for plotting, we use the log of the ratio
    vals = filter(x -> !isnan(x), vec(log_ratio))
    maxabs = isempty(vals) ? 1.0 : maximum(abs, vals)
    hm = heatmap!(ax, ps, p_bells, log_ratio; colormap = cgrad([:blue, :white, :red]), colorrange = (-maxabs, maxabs))
    Colorbar(fig[1, 2], hm; label = L"\log_{10}(\text{LIER}/p)", labelsize = 22)
    # Pseudothreshold frontier
    ε = 0.01
    iso_xs = Float64[]
    iso_ys = Float64[]
    for i in eachindex(ps)
        for j in eachindex(p_bells)
            if !isnan(log_ratio[i, j]) && abs(log_ratio[i, j]) < ε
                push!(iso_xs, ps[i])
                push!(iso_ys, p_bells[j])
            end
        end
    end
    if !isempty(iso_xs)
        # Order pseudothreshold locations and filter one data point per x-coordinate (lowest y value)
        order = sortperm(iso_xs)
        iso_xs = iso_xs[order]
        iso_ys = iso_ys[order]
        unique_indices = unique(i -> iso_xs[i], eachindex(iso_xs))
        iso_xs = iso_xs[unique_indices]
        iso_ys = iso_ys[unique_indices]
        _color = RGBf( 74/255, 163/255, 209/255)
        lines!(ax, iso_xs, iso_ys;
            color = _color, linewidth = 2, linestyle = :solid)
        scatter!(ax, iso_xs, iso_ys;
            color = _color, markersize = 10, marker = :circle,
            label = L"\log_{10}(\text{LIER}/p) \approx 0")
        axislegend(ax, position = :lt)
    end
    save(joinpath(data_path,"../..", "2d_heatmap_ratio.png"), fig)
    @info "2d_heatmap_ratio.png"
    # ---------------------- Plot 2D Ratio Heatmap ----------------------
    fig = Figure(size = (900, 600))
    ax = Axis(fig[1, 1]; 
            title = L"\text{Logical } |0\rangle_L\text{ Initialisation Error Rate (LIER) vs. } (p,p_{\text{Bell}})",
            subtitle  =  latexstring(_code_architecture_label * "; " *L"p_{\text{two-qubit gate}}=p_{\text{meas}}=p_{\text{init}}=p;"*" "*L"p_{\text{single-qubit gate}}="*"$(df.p_single_ratio[1])"*L"p;"*" "*L"p_{\text{idle}}="*"$(df.p_idle_ratio[1])"*L"p;"*" "*L"\text{depth_{telegate}}="*"$(df.telegate_idle_depth[1])"),
            xlabel = L"\text{Physical initialisation error rate } p",  
            ylabel = L"\text{Bell pair initialisation error rate }p_{\text{Bell}}",  
            xscale = log10, 
            yscale = log10,
            titlegap = 12,
            subtitlegap = 10,
            xgridvisible = true, 
            ygridvisible  = true, 
            xlabelsize = 22, 
            ylabelsize = 22,
            xticklabelsize = 20,  
            yticklabelsize = 20, 
            titlesize = 28,
            subtitlesize = 16,
            xminorgridvisible  = true, 
            yminorgridvisible  = true, 
            xminorticksvisible = true, 
            yminorticksvisible = true,
            xminorticks = IntervalsBetween(10),  
            yminorticks = IntervalsBetween(10),
    )
    log_LER = map(z -> (z > 0) ? log10(z) : NaN, LER) # for plotting, we use the log of the ratio
    vals = filter(x -> !isnan(x), vec(log_LER))
    maxabs = isempty(vals) ? 1.0 : maximum(abs, vals)
    hm = heatmap!(ax, ps, p_bells, log_LER; colormap = cgrad(:viridis, rev=true), colorrange = (-maxabs, 0))
    Colorbar(fig[1, 2], hm; label = L"\log_{10}(\text{LIER})", labelsize = 22)
    save(joinpath(data_path,"../..", "2d_heatmap.png"), fig)
    @info "2d_heatmap.png"
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


# ---------------------------------------------------------------
# ------------------ Pre-Decoding vs. logical  ------------------
# ---------------------------------------------------------------

# Plots the probability of having some physical X error (determined by anti-commutation with
# logical Z or some Z type stabiliser) vs. the probability of having a logical X error after decoding

function pre_decoding_vs_logical_plot(df, data_path)
    @info "Creating Pre-decoding error vs. logical error plot..."
    # ---------------------- Pre-Processing ----------------------
    df = copy(df)
    p_bells = sort(unique(Float64.(df.p_bell)))
    ps_all  = sort(unique(Float64.(df.p)))
    p_bells_plot = [minimum(p_bells), maximum(p_bells)]
    # ---------------------- Plotting ----------------------
    all_ys = filter(x -> x > 0, vcat(Float64.(df.logical_error_rate), Float64.(df.x_error_pre_decoding_rate), Float64.(df.z_error_pre_decoding_rate)))
    minexp = floor(Int, log10(minimum(all_ys)))
    maxexp = ceil(Int,  log10(maximum(all_ys)))
    ytick_vals   = 10.0 .^ (minexp:maxexp)
    ytick_labels = power_of_10_label.(ytick_vals,4)
    minexp_x = floor(Int, log10(minimum(ps_all)))
    maxexp_x = ceil(Int,  log10(maximum(ps_all)))
    xtick_vals   = 10.0 .^ (minexp_x:maxexp_x)
    xtick_labels = power_of_10_label.(xtick_vals,4)
    fig = Figure(size = (1600, 640), fontsize = 14)
    ax_kwargs = (
        xlabel     = L"\text{Physical initialisation error rate } p",
        ylabel     = L"\text{Error rate}",
        xscale     = log10, yscale = log10,
        xticks     = (xtick_vals, xtick_labels),
        yticks     = (ytick_vals, ytick_labels),
        titlegap = 12,
        subtitlegap = 10,
        xlabelsize = 22, ylabelsize = 22,
        titlesize = 20,
        xgridvisible = true, ygridvisible = true,
        xminorgridvisible = true, yminorgridvisible = true,
        xminorticksvisible = true, yminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        yminorticks = IntervalsBetween(10),
    )
    _code_architecture_label = df.code_architecture_label[1]
    Label(fig[0, 1:2], L"\text{Pre-Decoding and Logical Error Rates vs. Physical Initialisation Error}"; fontsize = 28, tellwidth = false)
    Label(fig[1, 1:2], latexstring(_code_architecture_label * "; " *L"p_{\text{two-qubit gate}}=p_{\text{meas}}=p_{\text{init}}=p;"*" "*L"p_{\text{single-qubit gate}}="*"$(df.p_single_ratio[1])"*L"p;"*" "*L"p_{\text{idle}}="*"$(df.p_idle_ratio[1])"*L"p;"*" "*L"\text{depth_{telegate}}="*"$(df.telegate_idle_depth[1])"); fontsize = 16, tellwidth = false)
    ax1 = Axis(fig[2, 1]; title = latexstring(L"p_{\text{Bell}}="*"$(power_of_10_label(p_bells_plot[1],4))"), ax_kwargs...)  
    ax2 = Axis(fig[2, 2]; title = latexstring(L"p_{\text{Bell}}="*"$(power_of_10_label(p_bells_plot[2],4))"), ax_kwargs...)
    col_pre_x   = RGBf(255/255, 177/255,  42/255) 
    col_pre_z   = RGBf( 0/255, 154/255, 207/255)
    col_logical = RGBf(255/255, 177/255,  42/255) 
    linestyles = [:solid, :solid, :dot]
    for (idx, p_bell) in enumerate(p_bells_plot)
        ax = idx == 1 ? ax1 : ax2
        df_b = copy(df[df.p_bell .== p_bell, :])
        df_ = df_b[df_b.logical_error_rate .> 0.0 .&& df_b.x_error_pre_decoding_rate .> 0.0 .&& df_b.z_error_pre_decoding_rate .> 0.0, :]
        @info "Deleted $(nrow(df_b)-nrow(df_)) row(s) because of 0.0 logical or pre-decoding error rate"
        sort!(df_, :p)
        lines!(ax, df_.p, df_.x_error_pre_decoding_rate; color = col_pre_x, linestyle=linestyles[1], linewidth = 2, label=L"\text{Pre-decoding physical X error rate}")
        scatter!(ax, df_.p, df_.x_error_pre_decoding_rate; color = col_pre_x, marker = :circle, markersize = 12)
        lines!(ax, df_.p, df_.z_error_pre_decoding_rate; color = col_pre_z, linestyle=linestyles[2], linewidth = 2, label=L"\text{Pre-decoding physical Z error rate}")
        scatter!(ax, df_.p, df_.z_error_pre_decoding_rate; color = col_pre_z, marker = :circle, markersize = 12)
        lines!(ax, df_.p, df_.logical_error_rate; color = col_logical,linestyle=linestyles[3], linewidth = 2, label=L"\text{Logical } |0\rangle_L \text{ initialisation error rate}"), 
        scatter!(ax, df_.p, df_.logical_error_rate; color = col_logical, marker = :diamond, markersize = 12)
    end
    linkyaxes!(ax1, ax2)
    hideydecorations!(ax2, grid = false, minorgrid = false)
    axislegend(ax1, position = :lt)
    axislegend(ax2, position = :rb)
    save(joinpath(data_path, "../..", "pre_decoding_x_vs_logical.png"), fig)
    @info "Saved pre_decoding_x_vs_logical.png"
end


# ---------------------------------------------------------------
# ------------------ Logical error vs. fidelity  ------------------
# ---------------------------------------------------------------

# Plots 1-LER for each qubit vs. the fidelity (no logical errors) <- only meaningful for k>1
# (the logical error rate for other plots is the mean of the k logical error rates per logical qubit)

function per_qubit_fidelity_bar_plot(df, data_path)
    # ---------------------- Pre-Processing ----------------------
    df = copy(df)
    p_bell = minimum(Float64.(df.p_bell)) # we plot the comparison only for the lowest bell pair noise
    df = df[Float64.(df.p_bell) .== p_bell, :]
    sort!(df, :p)
    # parsing logical_failures to compute LER per logical qubit
    parsed   = [parse.(Int, split(strip(s, ['[', ']']), ',')) for s in df.logical_failures]
    k_log_qubits = maximum(length.(parsed))
    accepted = Float64.(df.n_samples .- df.discarded_runs)
    for k in 1:k_log_qubits
        df[!, "ler_qubit_$k"] = [
            length(row) >= k ? row[k] / accepted[i] : NaN
            for (i, row) in enumerate(parsed)
        ]
    end
    ps     = Float64.(df.p)
    n_p    = length(ps)
    # ---------------------- Plotting ----------------------
    qubit_colors = [RGBf(255/255, 177/255,  42/255), RGBf( 0/255, 154/255, 207/255), RGBf(96/255, 184/255,  72/255) , RGBf(180/255, 120/255, 210/255)]
    log_ps     = log10.(ps)
    bar_width  = 0.7 * minimum(diff(log_ps)) 
    slot_width = bar_width / k_log_qubits
    offsets = [(i - (k_log_qubits + 1) / 2) * slot_width for i in 1:k_log_qubits]
    _code_architecture_label = df.code_architecture_label[1]#_code_architecture_label = replace(df.code_architecture_label[1], "_" => "\\_")
    fig = Figure(size = (max(600, n_p * 55), 520), fontsize = 13)
    ax = Axis(fig[1, 1],
        xlabel     = L"\text{Physical error rate } p",
        ylabel     = L"1 - \text{LIER or Fidelity}",
        title      = latexstring(L"\text{Logical qubit } (1 - \text{LIER}) \text{ and Fidelity vs. } p \text{ for }p_{\text{Bell}}="*"$(power_of_10_label(p_bell,4))"), 
        subtitle  =  latexstring(_code_architecture_label * "; " *L"p_{\text{two-qubit gate}}=p_{\text{meas}}=p_{\text{init}}=p;"*" "*L"p_{\text{single-qubit gate}}="*"$(df.p_single_ratio[1])"*L"p;"*" "*L"p_{\text{idle}}="*"$(df.p_idle_ratio[1])"*L"p;"*" "*L"\text{depth_{telegate}}="*"$(df.telegate_idle_depth[1])"),
        xlabelsize = 24, ylabelsize = 24,
        titlesize = 28, subtitlesize = 18,
        titlegap = 12,
        subtitlegap = 10,
        xgridvisible = false, ygridvisible = true,
        yminorgridvisible = true, yminorticksvisible = true,
        yminorticks = IntervalsBetween(10),
    )
    ax.xticks = (log_ps, power_of_10_label.(ps,4))
    ax.xticklabelrotation = π / 4
    for k in 1:k_log_qubits
        col   = "ler_qubit_$k"
        vals  = 1.0 .- Float64.(df[!, col])
        xs    = log_ps .+ offsets[k]
        barplot!(ax, xs, vals;
            width = slot_width * 0.9,
            color = qubit_colors[k],
            label = L"\text{Logical qubit } %$k",
        )
    end
     # fidelity as a dashed line across each bar group
    for (i, p) in enumerate(log_ps)
        x_lo = p + offsets[1]   - slot_width * 0.45
        x_hi = p + offsets[end] + slot_width * 0.45
        lines!(ax, [x_lo, x_hi], [df.avg_fidelity[i], df.avg_fidelity[i]];
            color     = :red,
            linestyle = :dash,
            linewidth = 2,
        )
    end
    lines!(ax, [NaN, NaN], [NaN, NaN];
        color     = :red,
        linestyle = :dash,
        linewidth = 2,
        label     = L"\text{Avg. fidelity}",
    )
    all_vals = vcat((1.0 .- df[!, "ler_qubit_$k"] for k in 1:k_log_qubits)..., df.avg_fidelity)
    all_vals = filter(x -> !isnan(x), Float64.(all_vals))
    ylo = floor(minimum(all_vals), digits = 4)
    ylims!(ax, max(0.0, ylo - 1/25*(1.0-ylo)), 1.0)
    axislegend(ax, position = :lb)
    save(joinpath(data_path, "../..", "per_qubit_ler_fidelity_bars.png"), fig)
    @info "Saved per_qubit_ler_fidelity_bars.png"
end


# ----------------------------------------------------------------------
# ------------------ Logical error amongst cores/codes  ----------------
# ----------------------------------------------------------------------


function compare_logical_error_rate(df)
    # ---------------------- Plotting ----------------------
    ps_all = sort(unique(Float64.(df.p)))
    p_range = range(minimum(ps_all), maximum(ps_all), length=300)
    labels  = unique(df.code_architecture_label) 
    colors = [RGBf(255/255, 177/255,  42/255), RGBf( 0/255, 154/255, 207/255), RGBf(96/255, 184/255,  72/255) , RGBf(180/255, 120/255, 210/255)]
    all_ys = filter(x -> x > 0, Float64.(df.logical_error_rate))
    minexp = floor(Int, log10(minimum(all_ys)))
    maxexp = ceil(Int,  log10(maximum(all_ys)))
    minexp_x = floor(Int, log10(minimum(ps_all)))
    maxexp_x = ceil(Int,  log10(maximum(ps_all)))
    xtick_vals   = 10.0 .^ (minexp_x:maxexp_x)
    xtick_labels = power_of_10_label.(xtick_vals,4)
    ytick_vals = 10.0 .^ (minexp:maxexp)
    ytick_labels = power_of_10_label.(ytick_vals,4)
    fig = Figure(size = (720, 540), fontsize = 14)
    ax  = Axis(fig[1,1],
               xscale = log10, 
               yscale = log10,
               ylabel = L"\text{Logical } |0\rangle_L \text{ initialisation error}",  
               title = L"\text{Logical vs. Physical Initialisation Error}",
               subtitle = latexstring(L"p_{\text{two-qubit gate}}=p_{\text{meas}}=p_{\text{init}}=p;"*" "*L"p_{\text{single-qubit gate}}="*"$(df.p_single_ratio[1])"*L"p;"*" "*L"p_{\text{idle}}="*"$(df.p_idle_ratio[1])"*L"p;"*" "*L"\text{depth_{telegate}}="*"$(df.telegate_idle_depth[1])"),
               xlabel = L"\text{Physical initialisation error rate } p", 
               xticks = (xtick_vals, xtick_labels), 
               yticks = (ytick_vals, ytick_labels), 
               ylabelsize = 22, xlabelsize = 22, titlesize = 28,
               subtitlesize=16,
               xgridvisible= true, 
               ygridvisible = true, 
               xminorgridvisible  = true, 
               yminorgridvisible  = true, 
               xminorticksvisible = true, 
               yminorticksvisible = true,
               xminorticks = IntervalsBetween(10), 
               yminorticks = IntervalsBetween(10),
               titlegap = 12,
               subtitlegap = 10)
    # p_L = p  (slope 1 in log-log)
    lines!(ax, p_range, p_range, color = :gray60,  linestyle = :dash, linewidth = 1.5,  label = L"p_L \sim p")
    # p_L = p²  (slope 2 in log-log)
    lines!(ax, p_range, p_range.^2, color = :gray40, linestyle = :dot, linewidth = 1.5, label = L"p_L \sim p^2")
    for (idx, config_label) in enumerate(labels)
        df_ = df[df.code_architecture_label .== config_label, :]
        p_bell_min = minimum(df_.p_bell)
        df_ = df_[df_.p_bell .== p_bell_min, :]
        sort!(df_, :p)
        df_.var = 1 ./ (df_.n_samples-df_.discarded_runs) .* (df_.logical_error_rate) .* (1 .- df_.logical_error_rate)
        df_.std_error = sqrt.(max.(df_.var, 0.0))
        data_label = L"\text{%$(config_label)}"
        errorbars!(ax, df_.p, df_.logical_error_rate, df_.std_error; color =colors[idx], label=idx == 1 ? L"\sigma(LER)" : nothing)
        scatterlines!(ax, df_.p, df_.logical_error_rate;
                        color = colors[idx], linewidth = 2, markersize = 8, label = data_label)
    end
    axislegend(ax, position = :rb)
    @info @__DIR__
    save(joinpath(@__DIR__, "../../data/comparison_ler.png"), fig)
    @info "Saved comparison_ler.png"
end


# ------------------------------------------
# ------------------ Helper ----------------
# ------------------------------------------


function power_of_10_label(val::Float64, round_to_digits::Int)
    val > 0 || return L""
    n = round(log10(val),digits=round_to_digits)
    return L"10^{%$n}"
end


struct Config
    code::String          
    architecture::String  
    code_architecture_label::String    # for plot legend
    data_path::String     # for data storage
end

function load_configs(configs::Vector{Config})::DataFrame
    dfs = map(configs) do cfg
        df = CSV.read(cfg.data_path, DataFrame)
        df.code         = fill(cfg.code,         nrow(df))
        df.architecture = fill(cfg.architecture, nrow(df))
        df.code_architecture_label = fill(cfg.code_architecture_label, nrow(df))
        df
    end
    return vcat(dfs...; cols = :union) 
end


# ------------------------------------------
# --------------- Execution ----------------
# ------------------------------------------


function analysis_dqc_sim(configs)
    df = load_configs(configs)
    if length(configs) == 1
        log_log_plot(df, configs[1].data_path) # can add optional keyword monolithic=true
        two_d_plot(df, configs[1].data_path)
        pre_decoding_vs_logical_plot(df, configs[1].data_path)
        per_qubit_fidelity_bar_plot(df, configs[1].data_path)
    else
        compare_logical_error_rate(df)
    end
end


configs_dqc_sim = [
    Config("Steane", "[4,3]", L"\text{Steane }[[7,1,3]], \text{ 2 cores}",  "data/Steane/[4, 3]/simulation_FT/dqc_sim_data.csv"),
    Config("Steane", "[4,3]", L"\text{Steane }[[7,1,3]], \text{ 2 cores}",  "data/Steane/[4, 3]/simulation_FT/dqc_sim_data.csv"),
    # Config("Shor", "[3,3,3]", L"\text{Shor }[[9,1,3]], \text{ 3 cores}",  "data/Shor/[3, 3, 3]/simulation_FT/dqc_sim_data.csv"),
    # Config("TrivariateBicycle", "[6,6]", L"[[12,2,3]], \text{ 2 cores}",  "data/TrivariateBicycle/[6, 6]/simulation_FT/dqc_sim_data.csv"),
    # Config("TrivariateBicycle", "[4,4,4]", L"[[12,2,3]], \text{ 3 cores}",  "data/TrivariateBicycle/[4, 4, 4]/simulation_FT/dqc_sim_data.csv"),
    # Config("TrivariateBicycle", "[3,3,3,3]", L"[[12,2,3]], \text{ 4 cores}",  "data/TrivariateBicycle/[3, 3, 3, 3]/simulation_FT/dqc_sim_data.csv"),
    
]

analysis_dqc_sim(configs_dqc_sim)



