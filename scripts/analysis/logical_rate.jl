# logical_rate.jl

using CairoMakie, LaTeXStrings
using CairoMakie: Axis
using LsqFit
using CSV, DataFrames
using Statistics

# ---------------------------------------------------------------
# ------------------ Scaling & Threshold ------------------------
# ---------------------------------------------------------------

function log_log_plot(df, data_path; ft = true, monolithic = false)
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
    fig = Figure(size = (830, 1120), fontsize = 14)
    ax1  = Axis(fig[1, 1],
        ylabel = L"\text{Logical } |0\rangle_L \text{ Initialisation Error Rate}",  
        xscale  = log10, 
        yscale  = log10,
        ytickformat = ys -> [L"%$(power_of_10_label(y,1))" for y in ys],
        xticklabelsize = 32, 
        yticklabelsize = 32,
        xlabelsize = 32,  
        ylabelsize = 30, 
        subtitlesize = 20,
        xgridvisible= true, 
        ygridvisible = true, 
        xminorticksvisible = true, 
        yminorticksvisible = true,
        xminorticks = IntervalsBetween(10), 
        yminorticks = IntervalsBetween(10)
    )
    ax2 = Axis(fig[2, 1],  
        xlabel = L"\text{Physical Initialisation Error Rate } p", 
        ylabel = L"\text{Acceptance Rate}",
        xscale = log10, 
        xtickformat = xs -> [L"%$(power_of_10_label(x,1))" for x in xs],
        xlabelsize = 32, 
        xticklabelsize = 32, 
        yticklabelsize = 32, 
        ylabelsize = 32,
    )
    ax2.ytickformat = values -> [L"%$(round(v, digits=3))" for v in values]
    lower_y_lim = ft ? min( minimum(df.acceptance_ratio), 0.85) : 0.95 
    upper_y_lim = ft ? 1.01 : 1.01
    ylims!(ax2, lower_y_lim , upper_y_lim)
    if maximum(df.distance) < 5
        # p_L = p  (slope 1 in log-log)
        lines!(ax1, p_range, p_range, color = :gray20,  linestyle = :dash, linewidth = 2.5)#,  label = L"p_L = p")
        text!(ax1, p_range[end], p_range[end], text = L"\text{LIER} = p", color = :gray20, align = (:right, :bottom), offset = (13, -50), fontsize = 24)

        # p_L = p²  (slope 2 in log-log)
        lines!(ax1, p_range, 200*p_range.^2, color = :gray40, linestyle = :dot, linewidth = 2.5)#, label = L"p_L = 100 p^2")
        text!(ax1, p_range[end], 200*p_range[end]^2, text = L"\text{LIER} \sim p^2", color = :gray40, align = (:right, :bottom), offset = (20, -65), fontsize = 24)
    else
        # p_L = p  (slope 1 in log-log for pseudothreshold)
        lines!(ax1, p_range, p_range, color = :gray20,  linestyle = :dash, linewidth = 2.5)#,  label = L"p_L = p")
        text!(ax1, p_range[end], p_range[end], text = L"\text{LIER} = p", color = :gray20, align = (:right, :bottom), offset = (13, -45), fontsize = 24)

        # p_L = p^3  (slope 3 in log-log)
        lines!(ax1, p_range, 500*p_range.^3, color = :gray40, linestyle = :dot, linewidth = 2.5)#, label = L"p_L = 500 p^3")
        text!(ax1, p_range[end], 500*p_range[end]^3, text = L"\text{LIER} \sim p^3", color = :gray40, align = (:right, :bottom), offset = (19, -65), fontsize = 24)

    end
    colors = [RGBf( 35/255, 139/255,  69/255), RGBf(230/255, 159/255,  0/255), RGBf( 33/255, 113/255, 181/255)]  
    for (idx,p_bell) in enumerate(p_bells_plot)
        df_bell = copy(df[(df.p_bell .== p_bell), :])
        @info "Bell state error: $p_bell"
        df_ = copy(df_bell[(df_bell.logical_error_rate.>0.0), :])
        @info "Deleted $(nrow(df_bell)-nrow(df_)) row(s) because of 0.0 logical error rate"
        sort!(df_, :p)
        df_.var = 1 ./ ( df_.k .* (df_.n_samples-df_.discarded_runs)) .* (df_.logical_error_rate) .* (1 .- df_.logical_error_rate)
        df_.std_error = sqrt.(max.(df_.var, 0.0))
        df_.var_AR = 1 ./ ( df_.n_samples) .* (df_.acceptance_ratio) .* (1 .- df_.acceptance_ratio)
        df_.std_error_AR = sqrt.(max.(df_.var_AR, 0.0))
        if idx ==1
            a, b = power_law_fit(df_.p, df_.logical_error_rate) # a = C, b = \alpha in the main text. We only determine scaling exponent for the lowest bell pair noise
            if monolithic # bell pair noise is irrelevant in this case
                fit_label = L"\text{Power law fit: p_L=%$(round(a, digits=2))\, p^{%$( round(b, digits=2))} }"
            else
                fit_label = latexstring(L"\text{LIER}="*"$(Int(round(a, digits=0)))p^{$( round(b, digits=2))}")
            end
            # fitted line:  p_L = a*p^b
            lines!(ax1, p_range, a .* p_range .^ b, color = colors[idx], linestyle = :dashdot, linewidth = 3.5, label = fit_label)
            pseudo_thresh = ft ? a^(1/(1-b)) : nothing
            @info "Pseudothreshold at $pseudo_thresh."
            if pseudo_thresh !== nothing
                vlines!(ax1, [pseudo_thresh], color= colors[idx], linestyle=:dash, linewidth = 3.0)
                vlines!(ax2, [pseudo_thresh], color= colors[idx], linestyle=:dash, linewidth = 3.0)
            end
            @info "Acceptance Rate for smallest p_Bell: $( minimum(df_.acceptance_ratio) ) up to $(maximum(df_.acceptance_ratio) )"
        end
        if monolithic
            data_label = L"\text{Logical Initialisation Error Rate}"
            acc_label  = L"\text{Acceptance Rate}"
        else
            data_label = latexstring(L"p_{\text{Bell}}="*"$(power_of_10_label(p_bell,1))") 
            acc_label  = latexstring(L"p_{\text{Bell}}="*"$(power_of_10_label(p_bell,1))") 
        end
        scatterlines!(ax1, df_.p, df_.logical_error_rate, color = colors[idx], markersize = 8, linewidth  = 2, label = data_label)
        scatterlines!(ax2, df_.p, df_.acceptance_ratio, color = colors[idx], markersize = 8, linewidth = 2,  label = acc_label )
        errorbars!(ax1, df_.p, df_.logical_error_rate, df_.std_error; color=colors[idx], whiskerwidth = 5) 
        errorbars!(ax2, df_.p, df_.acceptance_ratio, df_.std_error_AR; color=colors[idx], whiskerwidth = 5) 
    end
    linkxaxes!(ax1, ax2)
    hidexdecorations!(ax1, grid = false, minorgrid = false)
    axislegend(ax1, position = :rb, framevisible=true, labelsize=22)
    if ft
        save(joinpath(data_path, "../..", "qec_threshold_FT.png"), fig)
        save(joinpath(data_path, "../..", "qec_threshold_FT.pdf"), fig)
    else
        save(joinpath(data_path, "../..", "qec_threshold_non_FT.png"), fig)
        save(joinpath(data_path, "../..", "qec_threshold_non_FT.pdf"), fig)
    end
    @info "Saved qec_threshold"
end


# ---------------------------------------------------------------
# ------------------ 2D Heatmap ---------------------------------
# ---------------------------------------------------------------


function two_d_plot(df, data_path; ft = true)
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
    # ---------------------- Plot 2D Ratio Heatmap ----------------------
    fig = Figure(size = (800, 620), fontsize = 14)
    ax = Axis(fig[1, 1]; 
            xlabel = L"\text{Physical initialisation error rate } p",  
            ylabel = L"\text{Bell pair initialisation error rate }p_{\text{Bell}}",  
            xscale = log10, 
            yscale = log10, 
            xtickformat = xs -> [L"%$(power_of_10_label(x,3))" for x in xs],
            ytickformat = ys -> [L"%$(power_of_10_label(x,3))" for x in ys],
            xgridvisible = true, 
            ygridvisible  = true, 
            xlabelsize = 32, 
            ylabelsize = 32,
            xticklabelsize = 32,  
            yticklabelsize = 32, 
            xminorticks = IntervalsBetween(10),  
            yminorticks = IntervalsBetween(10),
    )
    log_ratio = map(z -> (z > 0.0) ? log10(z) : NaN, ratio) # for plotting, we use the log of the ratio
    vals = filter(x -> !isnan(x), vec(log_ratio))
    if isempty(vals)
        @error "filtered values are empty"
    else
        maxabs = maximum(abs, vals)
    end
    @info "Not showing $(length(log_ratio)-length(vals)) cells because of 0.0 logical error rate"
    p_edges = exp10.(range(log10(ps[1]), log10(ps[end]), length=length(ps)+1))
    pbell_edges = exp10.(range(log10(p_bells[1]), log10(p_bells[end]), length=length(p_bells)+1))
    hm = heatmap!(ax, p_edges, pbell_edges, log_ratio; colormap = cgrad([:blue, :white, :red]), colorrange = (-maxabs, maxabs))
    Colorbar(fig[1, 2], hm; label = L"\log_{10}(\text{LIER}/p)", labelsize = 28, ticklabelsize=20)
    pseudo_threshold_marker = PolyElement(color = :white, strokecolor = :gray, strokewidth = 0.5)
    axislegend(ax, [pseudo_threshold_marker], [L"\text{LIER} = p"], position = :rt, labelsize = 24)
    if ft
        save(joinpath(data_path,"../..", "2d_heatmap_ratio_FT.png"), fig)
        save(joinpath(data_path,"../..", "2d_heatmap_ratio_FT.pdf"), fig)
    else
        save(joinpath(data_path,"../..", "2d_heatmap_ratio_non_FT.png"), fig)
        save(joinpath(data_path,"../..", "2d_heatmap_ratio_non_FT.pdf"), fig)
    end
    @info "Saved 2d_heatmap_ratio"
    # ---------------------- Plot 2D Heatmap ----------------------
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1]; 
            xlabel = L"\text{Physical initialisation error rate } p",  
            ylabel = L"\text{Bell pair initialisation error rate }p_{\text{Bell}}",  
            xscale = log10, 
            yscale = log10,
            xtickformat = xs -> [L"%$(power_of_10_label(x,3))" for x in xs],
            ytickformat = ys -> [L"%$(power_of_10_label(x,3))" for x in ys],
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
            xminorticks = IntervalsBetween(10),  
            yminorticks = IntervalsBetween(10),
    )
    log_LER = map(z -> (z > 0) ? log10(z) : NaN, LER) # for plotting, we use the log of the ratio
    vals = filter(x -> !isnan(x), vec(log_LER))
    if isempty(vals)
        @error "filtered values are empty"
    else
        maxabs = maximum(abs, vals)
    end
    @info "Not showing $(length(log_LER)-length(vals)) cells because of 0.0 logical error rate"
    p_edges = exp10.(range(log10(ps[1]), log10(ps[end]), length=length(ps)+1))
    pbell_edges = exp10.(range(log10(p_bells[1]), log10(p_bells[end]), length=length(p_bells)+1))
    hm = heatmap!(ax, p_edges, pbell_edges, log_LER; colormap = cgrad(:viridis, rev=true), colorrange = (-maxabs,0.0))
    Colorbar(fig[1, 2], hm; label = L"\log_{10}(\text{LIER})", labelsize = 22)
    if ft
        save(joinpath(data_path,"../..", "2d_heatmap_FT.png"), fig)
        save(joinpath(data_path,"../..", "2d_heatmap_FT.pdf"), fig)
    else
        save(joinpath(data_path,"../..", "2d_heatmap_non_FT.png"), fig)
        save(joinpath(data_path,"../..", "2d_heatmap_non_FT.pdf"), fig)
    end
    @info "Saved 2d_heatmap"
end

function power_law_fit(physical_noise, logical_noise)
    m(t,p) = p[1] .* t .^p[2]# exp.(p[2] * t)
    p0 = [1.0, 1.0] 
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

function pre_decoding_vs_logical_plot(df, data_path; ft = true)
    @info "Creating Pre-decoding error vs. logical error plot..."
    # ---------------------- Pre-Processing ----------------------
    df = copy(df)
    p_bells = sort(unique(Float64.(df.p_bell)))
    p_bell_plot = minimum(p_bells)
    # ---------------------- Plotting ----------------------
    fig = Figure(size = (960, 640), fontsize = 14)
    ax_kwargs = (
        xlabel     = L"\text{Physical Initialisation Error Rate } p",
        ylabel     = L"\text{Error Rate}",
        xscale     = log10, 
        yscale     = log10,
        xtickformat = xs -> [L"%$(power_of_10_label(x,1))" for x in xs],
        ytickformat = ys -> [L"%$(power_of_10_label(y,1))" for y in ys],
        xticklabelsize = 32, yticklabelsize = 32, 
        xlabelsize = 32, ylabelsize = 32,
        xgridvisible = true, ygridvisible = true,
        xminorticksvisible = true, yminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        yminorticks = IntervalsBetween(10),
    )
    ax = Axis(fig[1,1]; ax_kwargs...)  
    col_pre_x   = RGBf( 33/255, 113/255, 181/255)#RGBf(255/255, 177/255,  42/255) 
    col_pre_z   =  RGBf(230/255, 159/255,  0/255)# RGBf( 0/255, 154/255, 207/255)
    col_logical = RGBf( 33/255, 113/255, 181/255)# RGBf(255/255, 177/255,  42/255) 
    linestyles = [:solid, :solid, :dot]
    df_b = copy(df[df.p_bell .== p_bell_plot, :])
    df_ = df_b[df_b.logical_error_rate .> 0.0 .&& df_b.x_error_pre_decoding_rate .> 0.0 .&& df_b.z_error_pre_decoding_rate .> 0.0, :]
    @info "Deleted $(nrow(df_b)-nrow(df_)) row(s) because of 0.0 logical or pre-decoding error rate"
    sort!(df_, :p)
    lines!(ax, df_.p, df_.x_error_pre_decoding_rate; color = col_pre_x, linestyle=linestyles[1], linewidth = 2, label=L"\text{Pre-Decoding Physical X-Error Rate}")
    scatter!(ax, df_.p, df_.x_error_pre_decoding_rate; color = col_pre_x, marker = :circle, markersize = 12)
    lines!(ax, df_.p, df_.z_error_pre_decoding_rate; color = col_pre_z, linestyle=linestyles[2], linewidth = 2, label=L"\text{Pre-Decoding Physical Z-Error Rate}")
    scatter!(ax, df_.p, df_.z_error_pre_decoding_rate; color = col_pre_z, marker = :circle, markersize = 12)
    lines!(ax, df_.p, df_.logical_error_rate; color = col_logical,linestyle=linestyles[3], linewidth = 2, label=L"\text{LIER}"),   #L"\text{Logical } |0\rangle_L^{\otimes k} \text{ Initialisation Error Rate}"), 
    scatter!(ax, df_.p, df_.logical_error_rate; color = col_logical, marker = :diamond, markersize = 12)
    axislegend(ax, position = :lt, framevisible=false, labelsize=24)
    if ft
        save(joinpath(data_path, "../..", "pre_decoding_x_vs_logical_FT.png"), fig)
        save(joinpath(data_path, "../..", "pre_decoding_x_vs_logical_FT.pdf"), fig)
    else 
        save(joinpath(data_path, "../..", "pre_decoding_x_vs_logical_non_FT.png"), fig)
        save(joinpath(data_path, "../..", "pre_decoding_x_vs_logical_non_FT.pdf"), fig)
    end
    @info "Saved pre_decoding_x_vs_logical"
end


# ---------------------------------------------------------------
# ------------------ Logical error vs. fidelity  ------------------
# ---------------------------------------------------------------

# Plots 1-LER for each qubit vs. the fidelity (no logical errors) <- only meaningful for k>1
# (the logical error rate for other plots is the mean of the k logical error rates per logical qubit)

function per_qubit_fidelity_bar_plot(df, data_path; ft = true)
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
    qubit_colors = [ RGBf( 33/255, 113/255, 181/255), RGBf(230/255, 159/255,  0/255), RGBf( 35/255, 139/255,  69/255) , RGBf(180/255, 120/255, 210/255)]
    log_ps     = log10.(ps)
    bar_width  = 0.7 * minimum(diff(log_ps)) 
    slot_width = bar_width / k_log_qubits
    offsets = [(i - (k_log_qubits + 1) / 2) * slot_width for i in 1:k_log_qubits]
    fig = Figure(size = (max(1100, n_p * 55), 650), fontsize = 13)
    ax = Axis(fig[1, 1],
        xlabel     = L"\text{Physical error rate } p",
        ylabel     = L"1 - \text{LIER or Fidelity (F)}",
        xlabelsize = 46, ylabelsize = 46,
        xticklabelsize = 36, yticklabelsize = 36, 
        xgridvisible = false, ygridvisible = false,
        yminorgridvisible = false, yminorticksvisible = false,
        yminorticks = IntervalsBetween(10),
    )
    ax.xticks = (log_ps, power_of_10_label.(ps,2))
    ax.ytickformat = values -> [L"%$(round(v, digits=5))" for v in values]
    ax.xticklabelrotation = π / 4
    for k in 1:k_log_qubits
        col   = "ler_qubit_$k"
        vals  = 1.0 .- Float64.(df[!, col])
        xs    = log_ps .+ offsets[k]
        barplot!(ax, xs, vals;
            width = slot_width * 0.9,
            color = qubit_colors[k],
            label = L"\text{LIER}_%$k",
        )
    end
     # fidelity as a dashed line across each bar group
    for (i, p) in enumerate(log_ps)
        x_lo = p + offsets[1]   - slot_width * 0.45
        x_hi = p + offsets[end] + slot_width * 0.45
        lines!(ax, [x_lo, x_hi], [df.avg_fidelity[i], df.avg_fidelity[i]];
            color     = :red,
            linestyle = :solid,
            linewidth = 3,
        )
    end
    lines!(ax, [NaN, NaN], [NaN, NaN];
        color     = :red,
        linestyle = :solid,
        linewidth = 5,
        label     = L"\text{Avg. F}",
    )
    all_vals = vcat((1.0 .- df[!, "ler_qubit_$k"] for k in 1:k_log_qubits)..., df.avg_fidelity)
    all_vals = filter(x -> !isnan(x), Float64.(all_vals))
    ylo = floor(minimum(all_vals), digits = 4)
    ylims!(ax, max(0.0, ylo - 1/25*(1.0-ylo)), 1.0)
    axislegend(ax, position = :rt, framevisible=false, labelsize = 32)
    if ft
        save(joinpath(data_path, "../..", "per_qubit_ler_fidelity_bars_FT.png"), fig)
        save(joinpath(data_path, "../..", "per_qubit_ler_fidelity_bars_FT.pdf"), fig)
    else
        save(joinpath(data_path, "../..", "per_qubit_ler_fidelity_bars_non_FT.png"), fig)
        save(joinpath(data_path, "../..", "per_qubit_ler_fidelity_bars_non_FT.pdf"), fig)
    end
    @info "Saved per_qubit_ler_fidelity_bars"
end


# ----------------------------------------------------------------------
# ------------------ Logical error amongst cores/codes  ----------------
# ----------------------------------------------------------------------

function compare_logical_error_rate(df)
    # ---------------------- Plotting ----------------------
    ps_all = sort(unique(Float64.(df.p)))
    labels  = unique(df.code_architecture_label) 
    colors = [RGBf(230/255, 159/255,   0/255),   # [[7,1,3]], 2 cores
          RGBf(221/255,  28/255, 119/255),   # [[9,1,3]], 3 cores
          RGBf(102/255, 194/255, 165/255),   # [[12,2,3]], 2 cores
          RGBf( 35/255, 139/255,  69/255),   # [[12,2,3]], 3 cores
          RGBf(  0/255,  68/255,  27/255),   # [[12,2,3]], 4 cores
          RGBf(166/255,  86/255,  40/255),   # [[17,1,5]], 2 cores
          RGBf(116/255, 169/255, 207/255),   # [[18,4,4]], 2 cores
          RGBf( 33/255, 113/255, 181/255),   # [[18,4,4]], 3 cores
          RGBf(  8/255,  48/255, 107/255),   # [[18,4,4]], 6 cores
          RGBf(166/255,  86/255,  40/255)]
    fig = Figure(size = (600, 540), fontsize = 14)
    ax  = Axis(fig[1,1],
        xscale = log10, 
        yscale = log10,
        ylabel = L"\text{Logical } |0\rangle_L \text{ Initialisation Error Rate}",  
        xlabel = L"\text{Physical Initialisation Error Rate } p", 
        xtickformat = xs -> [L"%$(power_of_10_label(x,1))" for x in xs],
        ytickformat = ys -> [L"%$(power_of_10_label(y,1))" for y in ys],
        ylabelsize = 24, xlabelsize = 24,
        xticklabelsize = 24, 
        yticklabelsize = 24, 
        xgridvisible= true, 
        ygridvisible = true, 
        xminorticksvisible = true, 
        yminorticksvisible = true,
        xminorticks = IntervalsBetween(10), 
        yminorticks = IntervalsBetween(10),
        titlegap = 12,
        subtitlegap = 10)
    p_range_ext = range(minimum(ps_all), maximum(ps_all)*1.45, length=300)
    # p_L = p  (slope 1 in log-log)
    lines!(ax, p_range_ext, p_range_ext, color = :gray20,  linestyle = :dash, linewidth = 2.5)
    text!(ax, p_range_ext[end], p_range_ext[end], text = L"\text{LIER} = p", color = :gray20, align = (:right, :bottom), offset = (13, -40), fontsize = 16)
    # p_L = p²  (slope 2 in log-log)
    lines!(ax, p_range_ext, 15*(p_range_ext).^2, color = :gray40, linestyle = :dot, linewidth = 2.5)#, label = L"\text{LIER} \sim p^2")
    text!(ax, p_range_ext[end], 15*p_range_ext[end]^2, text = L"\text{LIER} \sim p^2", color = :gray40, align = (:right, :top), offset = (17, -50), fontsize = 16)
    for (idx, config_label) in enumerate(labels)
        df_ = df[df.code_architecture_label .== config_label, :]
        p_bell_ = minimum(df_.p_bell)  # min p_bell
        df_ = df_[df_.p_bell .== p_bell_, :]
        sort!(df_, :p)
        df_.var = 1 ./ (  df_.k .* (  df_.n_samples-df_.discarded_runs) ) .* (df_.logical_error_rate) .* (1 .- df_.logical_error_rate)
        df_.std_error = sqrt.(max.(df_.var, 0.0))
        data_label = latexstring( config_label)
        errorbars!(ax, df_.p, df_.logical_error_rate, df_.std_error; color =colors[idx], whiskerwidth = 5) 
        scatterlines!(ax, df_.p, df_.logical_error_rate;
                        color = colors[idx], linewidth = 2, markersize = 8, label = data_label)
    end
    axislegend(ax, position = :rb, framevisible=false, labelsize = 18)
    save(joinpath(@__DIR__, "../../data/comparison_lier.png"), fig)
    save(joinpath(@__DIR__, "../../data/comparison_lier.pdf"), fig)
    @info "Saved comparison_ler"
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
    logical_q::Int
    distance::Int
    data_path::String     # for data storage
end

function load_configs(configs::Vector{Config})::DataFrame
    dfs = map(configs) do cfg
        df = CSV.read(cfg.data_path, DataFrame)
        df.code         = fill(cfg.code,         nrow(df))
        df.architecture = fill(cfg.architecture, nrow(df))
        df.k = fill(cfg.logical_q, nrow(df))
        df.distance = fill(cfg.distance, nrow(df))
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
        log_log_plot(df, configs[1].data_path; ft=true, monolithic=false) # can add optional keyword monolithic=true, ft_encoding = true
        two_d_plot(df, configs[1].data_path; ft=true)
        pre_decoding_vs_logical_plot(df, configs[1].data_path; ft=true)
        per_qubit_fidelity_bar_plot(df, configs[1].data_path; ft=true)
    else
        compare_logical_error_rate(df)
    end
end

# select one or more code-architecture pair(s) for plotting
configs_dqc_sim = [
    Config("Steane", "[4,3]", L"[[7,1,3]], \text{ 2 cores}", 1, 3,  "data/Steane/[4, 3]/simulation_FT/dqc_sim_data.csv"),
    Config("Shor", "[3,3,3]", L"[[9,1,3]], \text{ 3 cores}", 1, 3, "data/Shor/[3, 3, 3]/simulation_FT/dqc_sim_data.csv"),
    Config("TrivariateBicycle", "[6,6]", L"[[12,2,3]], \text{ 2 cores}", 2, 3,  "data/TrivariateBicycle/[6, 6]/simulation_FT/dqc_sim_data.csv"),
    Config("TrivariateBicycle", "[4,4,4]", L"[[12,2,3]], \text{ 3 cores}", 2,3, "data/TrivariateBicycle/[4, 4, 4]/simulation_FT/dqc_sim_data.csv"),
    #Config("TrivariateBicycle", "[4,4,4]", L"[[12,2,3]], \text{ 3 cores}", 2,3, "data/TrivariateBicycle/[4, 4, 4]/simulation_non_FT/dqc_sim_data.csv"),  # need to set ft=false above
    Config("TrivariateBicycle", "[3,3,3,3]", L"[[12,2,3]], \text{ 4 cores}", 2,3, "data/TrivariateBicycle/[3, 3, 3, 3]/simulation_FT/dqc_sim_data.csv"),
    Config("Color", "[8,9]", L"[[17,1,5]], \text{ 2 cores}",  4, 5,"data/Triangular/[8, 9]/simulation_FT/dqc_sim_data.csv"),
    Config("BivariateBicycle", "[9,9]", L"[[18,4,4]], \text{ 2 cores}",  4,4, "data/BivariateBicycle/[9, 9]/simulation_FT/dqc_sim_data.csv"),
    Config("BivariateBicycle", "[6,6,6]", L"[[18,4,4]], \text{ 3 cores}",  4,4, "data/BivariateBicycle/[6, 6, 6]/simulation_FT/dqc_sim_data.csv"),
    Config("BivariateBicycle", "[3,3,3,3,3,3]", L"[[18,4,4]], \text{ 6 cores}", 4,4,  "data/BivariateBicycle/[3, 3, 3, 3, 3, 3]/simulation_FT/dqc_sim_data.csv"),
]

analysis_dqc_sim(configs_dqc_sim)

