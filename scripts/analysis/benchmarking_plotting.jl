using CairoMakie, LaTeXStrings
using CairoMakie: Axis
using LsqFit
using CSV, DataFrames
using Statistics

function dqc_simulation_benchmarking(df)
    @info "Creating DQC simulation scaling plot..."
     # ---------------------- Pre-Processing ----------------------
    sort!(df, :num_samples)
    xs   = Float64.(df.num_samples)
    ys   = Float64.(df.mean_time_s)
    errs = Float64.(df.std_time_s)
    # ---------------------- Plotting ----------------------
    color_line  = RGBAf( 0/255, 154/255, 207/255)
    fig = Figure(size = (720, 480), fontsize = 14)
    ax = Axis(fig[1, 1];
        title          = L"\text{DQC Simulation Runtime vs. Number of Samples}",
        xlabel         = L"\text{Number of samples}",
        ylabel         = L"\text{Runtime (seconds)}",
        titlesize      = 28,
        subtitlesize   = 18,
        titlegap       = 10,
        subtitlegap = 5,
        xlabelsize     = 28,
        ylabelsize     = 28,
        xticklabelsize = 24,
        yticklabelsize = 24,
        ytickformat = ys -> [L"{%$y}" for y in ys],
        xgridvisible   = true,
        ygridvisible   = true,
        yminorticks    = IntervalsBetween(5),
    )
    lines!(ax, xs, ys;
        color     = color_line,
        linewidth = 2,
        label     = L"\text{Runtime}\pm 1\sigma",
    )
    scatter!(ax, xs, ys;
        color     = color_line,
        markersize = 8,
    )
    errorbars!(ax, xs, ys, errs;
        color     = color_line,
        linewidth = 1.5,
        whiskerwidth = 5
    )
    ylims!(ax, 0.1, nothing)
    axislegend(ax; position = :lt, framevisible = false, labelsize = 24)
    save("data/dqc_simulation_scaling.png", fig)
    save("data/dqc_simulation_scaling.pdf", fig)
    @info "Saved dqc_simulation_scaling to data/"
    return fig
end

function power_of_10_label(val::Float64, round_to_digits::Int)
    val > 0 || return L""
    n = round(log10(val),digits=round_to_digits)
    return L"10^{%$n}"
end

df = CSV.read("data/dqc_sim_scaling.csv", DataFrame)
dqc_simulation_benchmarking(df)
