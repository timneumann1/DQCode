# optimiser_evolution.jl

using CSV, DataFrames
using CairoMakie


function plot_evolution(dir, optimiser_label::String, fitness_scores, fidelities, single_q_counts,two_q_counts, telegate_counts )#, genetic_params::GeneticParameters)
    title_str = "Evolution of optimiser metrics for $optimiser_label"         
    fig = Figure(size = (800, 900))
    if optimiser_label == "Warm-Start Genetic Algorithm"
        ax_fit   = Axis(fig[1, 1], ylabel="Fitness", title=title_str)
    elseif optimiser_label == "Monte Carlo Tree Search"
        ax_fit   = Axis(fig[1, 1], ylabel="Reward", title=title_str)
    end
    ax_gates = Axis(fig[2, 1], ylabel="Gate Counts")
    ax_fid   = Axis(fig[3, 1], xlabel="Generation", ylabel="Fidelity")
    generations = 1:length(fitness_scores)
    lines!(ax_fit, generations, fitness_scores, color=:seagreen4, linewidth=2)
    lines!(ax_gates, generations, single_q_counts, label="Single-qubit", color=:goldenrod, linewidth=2)
    lines!(ax_gates, generations, two_q_counts,    label="Two-qubit",    color=:steelblue,  linewidth=2)
    lines!(ax_gates, generations, telegate_counts, label="Telegates",    color=:crimson, linewidth=2)
    axislegend(ax_gates, position=:rt) 
    lines!(ax_fid, generations, fidelities, color=:green, linewidth=2)
    ylims!(ax_fid, minimum(fidelities) - 0.05, 1.05)
    linkxaxes!(ax_fit, ax_gates, ax_fid)
    hidexdecorations!(ax_fit, grid=false)
    hidexdecorations!(ax_gates, grid=false)
    rowgap!(fig.layout, 10)
    outpath = joinpath(dir, "optimisation_evolution.png")
    save(outpath, fig)
end

code = "Steane"
qpu_sizes = "[4, 3]"
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