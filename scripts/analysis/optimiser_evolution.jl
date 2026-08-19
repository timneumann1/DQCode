# optimiser_evolution.jl

using CSV, DataFrames
using CairoMakie


function plot_evolution(dir, optimiser_label::String, fitness_scores, fidelities, single_q_counts,two_q_counts, telegate_counts )#, genetic_params::GeneticParameters)
    fig = Figure(size = (800, 900), fontsize = 18)
    if optimiser_label == "Warm-Start Genetic Algorithm"
        ax_fit   = Axis(fig[1, 1], ylabel=L"\text{Fitness}", 
                            ylabelsize = 36, xlabelsize = 36,
                            xticklabelsize =28, yticklabelsize = 28)
    elseif optimiser_label == "Monte Carlo Tree Search"
        ax_fit   = Axis(fig[1, 1], ylabel=L"\text{Reward}", 
                            titlegap = 8,ylabelsize = 36, xlabelsize = 36,
                            xticklabelsize = 28, yticklabelsize = 28)
    end
    ax_gates = Axis(fig[2, 1], ylabel=L"\text{Gate Counts}",# titlesize = 26,
                            ylabelsize = 36, xlabelsize = 36,
                            xticklabelsize = 28, yticklabelsize = 28)
    ax_fid   = Axis(fig[3, 1], xlabel=L"\text{Generation}", ylabel=L"\text{Fidelity}",
                            ylabelsize = 36, xlabelsize = 36,
                            xticklabelsize = 28, yticklabelsize = 28)
    generations = 1:length(fitness_scores)
    lines!(ax_fit, generations, fitness_scores, color=:seagreen4, linewidth=3)
    lines!(ax_gates, generations, single_q_counts, label=L"\text{Hadamard}", color=:goldenrod, linewidth=3)
    lines!(ax_gates, generations, two_q_counts,    label=L"\text{Intra-Core}",    color=:steelblue,  linewidth=3)
    lines!(ax_gates, generations, telegate_counts, label=L"\text{Inter-Core}",    color=:crimson, linewidth=3)
    axislegend(ax_gates, position=:rc, labelsize=24, framevisible = false,
           foreground_color_legend = nothing) 
    lines!(ax_fid, generations, fidelities, color=:green, linewidth=3)
    ylims!(ax_fid, minimum(fidelities) - 0.05, 1.05)
    linkxaxes!(ax_fit, ax_gates, ax_fid)
    hidexdecorations!(ax_fit, grid=false)
    hidexdecorations!(ax_gates, grid=false)
    rowgap!(fig.layout, 10)
    outpath = joinpath(dir, "optimisation_evolution.pdf")
    save(outpath, fig)
end

code = "TrivariateBicycle"   # QEC Code
qpu_sizes = "[6, 6]"         # QPU Size

# ------ GA ------ 
ga_dir = joinpath(@__DIR__, "..", "..", "data", "$code/$qpu_sizes", "warmstart_ga")
ga_evol = CSV.read(joinpath(ga_dir,"genetic_evolution.csv"), DataFrame)
fitness_evolution = ga_evol.fitness_evolution
fidelity_evolution = ga_evol.fidelity_evolution
single_count = ga_evol.single_count
two_qubit_count = ga_evol.two_qubit_count
telegate_count = ga_evol.telegate_count
plot_evolution(ga_dir, "Warm-Start Genetic Algorithm", fitness_evolution, fidelity_evolution, single_count,two_qubit_count,telegate_count)
# ------ MCTS ------
mcts_dir = joinpath(@__DIR__, "..", "..", "data", "$code/$qpu_sizes", "mcts")
mcts_evol = CSV.read(joinpath(mcts_dir,"mcts_evolution.csv"), DataFrame)
fidelity_evolution = mcts_evol.fidelity_evolution
single_count = mcts_evol.single_count
two_qubit_count = mcts_evol.two_qubit_count
telegate_count = mcts_evol.telegate_count
reward_evolution = mcts_evol.reward_evolution
plot_evolution(mcts_dir, "Monte Carlo Tree Search", reward_evolution, fidelity_evolution, single_count,two_qubit_count,telegate_count)