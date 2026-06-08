# logical_rate.jl

using CairoMakie, LaTeXStrings
using CairoMakie: Axis
using LsqFit
using CSV, DataFrames
using Statistics


# ------------------------------------------
# --------- Gate count comparison ----------
# ------------------------------------------

function compare_telegate_counts(df)
    @info "Creating gate count plot..."
    # ---------------------- Pre-Processing ----------------------
    df.telegate_counts   = [parse.(Int, split(strip(s, ['[', ']']), ','))[3] for s in df.gate_counts]
    # ---------------------- Plotting ----------------------
    opt_methods = ["qiskit_encoding", "mqt_encoding", "Gottesman", "dqc_compiled", "GA", "MCTS"]
    colors = [
        RGBf(  0/255, 154/255, 207/255),     
        RGBf(200/255, 173/255, 127/255),
        RGBf(255/255, 20/255, 64/255),
        RGBf(50/255, 50/255, 55/255),
        RGBf(255/255, 177/255,  42/255), 
        RGBf( 96/255, 184/255,  72/255),
    ]
    method_labels = ["Qiskit A-G", "MQT-QECC", "Gottesman (standard)", "Gottesman (compiled)", "warmstart GA", "MCTS"]
    method_to_color_labels = Dict(m => (colors[i], method_labels[i]) for (i, m) in enumerate(opt_methods))
    code_architecture_labels = unique(df.code_architecture_label)   
    num_code_architectures  = length(code_architecture_labels)
    n_methods   = length(opt_methods)
    group_width = 0.7
    bar_w       = group_width / n_methods
    offsets     = [(i - (n_methods + 1) / 2) * bar_w for i in 1:n_methods]
    code_architecture_positions = collect(1:num_code_architectures)  
    fig = Figure(size = (max(640, num_code_architectures * 160), 560), fontsize = 14)
    ax = Axis(fig[1, 1];
        title          = L"\text{Telegate Count per Optimisation Method and Code–Architecture Pair}",
        ylabel         = L"\text{Telegate count}",
        titlesize      = 24,
        titlegap       = 12,
        xlabelsize     = 20,
        ylabelsize     = 20,
        xticks         = (Float64.(code_architecture_positions), latexstring.(code_architecture_labels)),
        xticklabelsize = 13,
        xgridvisible   = false,
        ygridvisible   = true,
        yminorgridvisible  = true,
        yminorticksvisible = true,
        yminorticks    = IntervalsBetween(10),
    )
    for (mi, method) in enumerate(opt_methods)
        xs   = Float64[]
        vals = Float64[]
        for (ai, ca_lbl) in enumerate(code_architecture_labels)
            row = df[(df.code_architecture_label .== ca_lbl) .& (df.method .== method), :]
            isempty(row) && continue
            push!(xs,   code_architecture_positions[ai] + offsets[mi])
            push!(vals, Float64(row.telegate_counts[1]))
        end
        isempty(xs) && continue
        barplot!(ax, xs, vals;
            width = bar_w * 0.92,
            color = method_to_color_labels[method][1],
            label = L"\text{%$(method_to_color_labels[method][2])}",
        )
    end
    xlims!(ax, 0.5, num_code_architectures + 0.5)
    ylims!(ax, 0, nothing)
    axislegend(ax, position = :rt, framevisible = true)
    save(joinpath(@__DIR__, "../../data/gate_count_comparison.png"), fig)
    @info "Saved gate_count_comparison.png"
end


# ----------------------------------------------------------------------
# ------------------ Resource comparison ----------------
# ----------------------------------------------------------------------

function compare_resources(df, data_path)
    @info "Creating circuit vs. measurement resource comparison plot..."
    # ---------------------- Plotting ----------------------
    encoding_methods = ["encoding_circ", "stabiliser_measurement"]
    colors = [
        RGBf(  0/255, 154/255, 207/255),     
        RGBf(255/255, 177/255,  42/255), 
    ]
    encoding_labels = [latexstring(L"\text{FT encoding circuit (acceptance rate=}"*"$( round( df.acceptance_ratio[1], digits=3))"*L"\text{ for } p="*"$(power_of_10_label(df.p[1], 3))"*L"\text{, }p_{\text{Bell}}="*"$( power_of_10_label(df.p_bell[1], 3) ))"), 
                        latexstring(L"\text{Distributed Stabiliser Measurements (per round)}")]
    function get_val(method, col)
        rows = df[df.method .== method, :]
        isempty(rows) && return NaN
        v = rows[1, col]
        ismissing(v) && return NaN
        return Float64(v)
    end
    row1_metrics = [
        ("single_qubit_gates", L"\text{Single-qubit gates}", L"\text{Count}"),
        ("cx_gates",           L"\text{CX gates}", L"\text{Count}"),
        ("telegates",          L"\text{Telegates}", L"\text{Count}")
    ]
    row2_metrics = [
        ("measurements",       L"\text{Measurements}", L"\text{# of measurements}"),
        ("depth_cx_layers",       L"\text{Depth: CX Layers}", L"\text{# of layers}"),
        ("depth_telegate_layers", L"\text{Depth: Telegate Layers}", L"\text{# of layers}"),
    ]
    bar_xs = [1.0, 2.0] 
    common_ax_kwargs = (
        xgridvisible       = false,
        ygridvisible       = true,
        yminorgridvisible  = true,
        yminorticksvisible = true,
        yminorticks        = IntervalsBetween(10),
        ylabelsize         = 16,
        titlesize          = 16,
    ) 
    _code_architecture_label = df.code_architecture_label[1]
    fig = Figure(size = (1200, 760), fontsize = 14)
    Label(fig[1, 1:3], L"\text{Resource Comparison: Postselected FT encoding circuit vs. distributed stabiliser measurement initialisation}",
            fontsize = 24, tellwidth = false)
    Label(fig[2, 1:3], latexstring("$(_code_architecture_label)"),
            fontsize = 18, tellwidth = false)
    function draw_panel!(fig_pos, col_key, title_str, ylabel_str; show_legend = false)
        ax = Axis(fig_pos;
            title = title_str,
            ylabel = ylabel_str,
            common_ax_kwargs...,
        )
        xlims!(ax, 0.5, 2.5)
        for (mi, method) in enumerate(encoding_methods)
            v = get_val(method, col_key)
            isnan(v) && continue
            barplot!(ax, [bar_xs[mi]], [v];
                width = 0.5,
                color = colors[mi],
                label = encoding_labels[mi]
            )
        end
        ylims!(ax, 0, nothing)
        show_legend && Legend(fig[5, 1:3], ax,  framevisible = true, labelsize = 16)
        return ax
    end
    for (ci, (col_key, title_str, y_label)) in enumerate(row1_metrics)
        draw_panel!(fig[3, ci], col_key, title_str, y_label)
    end
    for (ci, (col_key, title_str, y_label)) in enumerate(row2_metrics)
        draw_panel!(fig[4, ci], col_key, title_str, y_label, ; show_legend = (ci == 1))
    end
    #olsize!(fig.layout, 5, Fixed(1))
    rowgap!(fig.layout, 4,5) 
    save(joinpath(data_path, "../../resource_comparison.png"), fig)
    @info "Saved resource_comparison.png"
end



# ------------------------------------------
# ------------------ Helper ----------------
# ------------------------------------------

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

function power_of_10_label(val::Float64, round_to_digits::Int)
    val > 0 || return L""
    n = round(log10(val),digits=round_to_digits)
    return L"10^{%$n}"
end


# ------------------------------------------
# --------------- Execution ----------------
# ------------------------------------------


function analysis_optimiser(configs)
    # Compare the telegate counts per optimisation method for a specific code-architecture configuration
    df = load_configs(configs)
    @info df
    compare_telegate_counts(df)
end


function analysis_resources(configs)
    # Compare the resources needed for FT circuit vs. measurement for a specific code-architecture configuration
    df = load_configs(configs)
    @info df
    compare_resources(df, configs[1].data_path)
end


configs_optimiser = [
    Config("Steane", "[4,3]", L"\text{Steane }[[7,1,3]], \text{ 2 cores}",  "data/Steane/[4, 3]/qiskit_encoding/qiskit_encoding_stats.csv"),
    Config("Steane", "[4,3]", L"\text{Steane }[[7,1,3]], \text{ 2 cores}",  "data/Steane/[4, 3]/mqt_encoding/mqt_encoding_stats.csv"),
    Config("Steane", "[4,3]", L"\text{Steane }[[7,1,3]], \text{ 2 cores}",  "data/Steane/[4, 3]/warmstart_ga/warm_start_ga_stats.csv"),
    Config("Steane", "[4,3]", L"\text{Steane }[[7,1,3]], \text{ 2 cores}",  "data/Steane/[4, 3]/mcts/mcts_stats.csv"),

    Config("Shor", "[3,3,3]", L"\text{Shor }[[9,1,3]], \text{ 3 cores}",  "data/Shor/[3, 3, 3]/qiskit_encoding/qiskit_encoding_stats.csv"),
    Config("Shor", "[3,3,3]", L"\text{Shor }[[9,1,3]], \text{ 3 cores}",  "data/Shor/[3, 3, 3]/mqt_encoding/mqt_encoding_stats.csv"),
    Config("Shor", "[3,3,3]", L"\text{Shor }[[9,1,3]], \text{ 3 cores}",  "data/Shor/[3, 3, 3]/warmstart_ga/warm_start_ga_stats.csv"),
    Config("Shor", "[3,3,3]", L"\text{Shor }[[9,1,3]], \text{ 3 cores}",  "data/Shor/[3, 3, 3]/mcts/mcts_stats.csv"),

    Config("TrivariateBicycle", "[6,6]", L"[[12,2,3]], \text{ 2 cores}",  "data/TrivariateBicycle/[6, 6]/qiskit_encoding/qiskit_encoding_stats.csv"),
    Config("TrivariateBicycle", "[6,6]", L"[[12,2,3]], \text{ 2 cores}",  "data/TrivariateBicycle/[6, 6]/mqt_encoding/mqt_encoding_stats.csv"),
    Config("TrivariateBicycle", "[6,6]", L"[[12,2,3]], \text{ 2 cores}",  "data/TrivariateBicycle/[6, 6]/warmstart_ga/warm_start_ga_stats.csv"),
    Config("TrivariateBicycle", "[6,6]", L"[[12,2,3]], \text{ 2 cores}",  "data/TrivariateBicycle/[6, 6]/mcts/mcts_stats.csv"),

    Config("TrivariateBicycle", "[4,4,4]", L"[[12,2,3]], \text{ 3 cores}",  "data/TrivariateBicycle/[4, 4, 4]/qiskit_encoding/qiskit_encoding_stats.csv"),
    Config("TrivariateBicycle", "[4,4,4]", L"[[12,2,3]], \text{ 3 cores}",  "data/TrivariateBicycle/[4, 4, 4]/mqt_encoding/mqt_encoding_stats.csv"),
    Config("TrivariateBicycle", "[4,4,4]", L"[[12,2,3]], \text{ 3 cores}",  "data/TrivariateBicycle/[4, 4, 4]/warmstart_ga/warm_start_ga_stats.csv"),
    Config("TrivariateBicycle", "[4,4,4]", L"[[12,2,3]], \text{ 3 cores}",  "data/TrivariateBicycle/[4, 4, 4]/mcts/mcts_stats.csv"),

    Config("TrivariateBicycle", "[3,3,3,3]", L"[[12,2,3]], \text{ 4 cores}",  "data/TrivariateBicycle/[3, 3, 3, 3]/qiskit_encoding/qiskit_encoding_stats.csv"),
    Config("TrivariateBicycle", "[3,3,3,3]", L"[[12,2,3]], \text{ 4 cores}",  "data/TrivariateBicycle/[3, 3, 3, 3]/mqt_encoding/mqt_encoding_stats.csv"),
    Config("TrivariateBicycle", "[3,3,3,3]", L"[[12,2,3]], \text{ 4 cores}",  "data/TrivariateBicycle/[3, 3, 3, 3]/warmstart_ga/warm_start_ga_stats.csv"),
    Config("TrivariateBicycle", "[3,3,3,3]", L"[[12,2,3]], \text{ 4 cores}",  "data/TrivariateBicycle/[3, 3, 3, 3]/mcts/mcts_stats.csv"),
]

#analysis_optimiser(configs_optimiser)


configs_resources = [
    Config("Steane", "[4,3]", L"\text{Steane }[[7,1,3]], \text{ 2 cores}",  "data/Steane/[4, 3]/simulation_FT/resources_info_circ.csv"),
    Config("Steane", "[4,3]", L"\text{Steane }[[7,1,3]], \text{ 2 cores}",  "data/Steane/[4, 3]/simulation_FT/resources_info_meas.csv"),
]

analysis_resources(configs_resources)




