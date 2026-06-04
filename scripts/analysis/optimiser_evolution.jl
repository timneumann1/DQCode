using CSV, DataFrames
using CairoMakie


function plot_evolution(dir, optimiser_label::String, fitness_scores, fidelities, single_q_counts,two_q_counts, telegate_counts )#, genetic_params::GeneticParameters)
    title_str = "Evolution of optimiser metrics for $optimiser_label"     
    #subtitle_str = "$(genetic_params.num_individuals) individuals over $(genetic_params.num_generations) generations"
    
    fig = Figure(size = (800, 900))

    ax_fit   = Axis(fig[1, 1], ylabel="Fitness", title=title_str)#, subtitle = subtitle_str)
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

    outpath = joinpath(dir, "Optimisation_Evolution.png")
    save(outpath, fig)
end

function plot_evolution(dir, optimiser_label::String, fidelities, gate_counts)#, mcts_params::MCTSParameters)
    
    title_str = "Evolution of optimiser metrics for $optimiser_label"     
    subtitle_str = "$(mcts_params.n_iterations) Iterations over depth $(mcts_params.depth), exploration_constant $(mcts_params.exploration_constant)"
    
    fig = Figure(size = (800, 900))

    ax_gates = Axis(fig[1, 1], ylabel="Gate Counts",  title=title_str, subtitle = subtitle_str)
    ax_fid   = Axis(fig[2, 1], xlabel="Step", ylabel="Fidelity")

    generations = 1:length(fidelities)
    lines!(ax_fid, generations, fidelities, color=:seagreen4, linewidth=2)
    ylims!(ax_fid, -0.05, 1.05) 

    tick_labels = string.(generations)
    ax_gates.xticks = (generations, tick_labels)
    ax_fid.xticks = (generations, tick_labels)

    ax_gates.xticklabelrotation = pi/4
    ax_fid.xticklabelrotation = pi/4

    single_q_counts = [g[1] for g in gate_counts]
    two_q_counts    = [g[2] for g in gate_counts]
    telegate_counts = [g[3] for g in gate_counts]

    lines!(ax_gates, generations, single_q_counts, label="Single-qubit", color=:goldenrod, linewidth=2)
    lines!(ax_gates, generations, two_q_counts,    label="Two-qubit",    color=:steelblue,  linewidth=2)
    lines!(ax_gates, generations, telegate_counts, label="Telegates",    color=:crimson, linewidth=2)
    axislegend(ax_gates, position=:rt) 

    linkxaxes!(ax_gates, ax_fid)
    hidexdecorations!(ax_gates, grid=false)

    rowgap!(fig.layout, 10)

    outpath = joinpath(dir, "optimisation_evolution.png")
    save(outpath, fig)
end

code = "Steane"
qpu_sizes = "[4, 3]"
ga_dir = joinpath(@__DIR__, "..", "..", "data", "$code/$qpu_sizes", "warmstart_ga")
ga_evol = CSV.read(joinpath(ga_dir,"genetic_evolution.csv"), DataFrame)

fitness_evolution = ga_evol.fitness_evolution
fidelity_evolution = ga_evol.fidelity_evolution
single_count = ga_evol.single_count
two_qubit_count = ga_evol.two_qubit_count
telegate_count = ga_evol.telegate_count

plot_evolution(ga_dir, "Warm-Start Genetic Algorithm", fitness_evolution, fidelity_evolution, single_count,two_qubit_count,telegate_count)



#mcts_evol = CSV.read(joinpath(dir, "mcts/mcts_evolution.csv"), DataFrame)

#fidelity_evolution = mcts_evol.fidelity_evolution
#gate_count_evolution = mcts_evol.gate_count_evolution

#plot_evolution(mcts_dir, "MCTS", fidelity_evolution, gate_count_evolution, mcts_params)




