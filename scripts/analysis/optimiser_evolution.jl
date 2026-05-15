include(joinpath(@__DIR__, "..", "src", "DQCode.jl"))
include(joinpath(@__DIR__, "..", "src", "Helper.jl"))

using .DQCode
using .Helper
using CSV
# specify csv script of interest

dir = "/Users/tim/Tim/projects/DQCode/data/Shor/[3, 3, 3]"

ga_evol = CSV.read(joinpath(dir, "/warmstart_ga/genetic_evolution.csv"))
mcts_evol = CSV.read(joinpath(dir, "/mcts/mcts_evolution.csv"))

fitness_evolution = ga_evol.fitness_evolution
fidelity_evolution = ga_evol.fidelity_evolution
single_count = ga_evol.single_count
two_qubit_count = ga_evol.ga_evol.single_count
telegate_count = ga_evol.telegate_count

Helper.plot_evolution(ga_dir, "Warm-Start Genetic Algorithm", fitness_evol, fid_evol, gc_evol, gen_params)

fidelity_evolution = mcts_evol.fidelity_evolution
gate_count_evolution = mcts_evol.gate_count_evolution

Helper.plot_evolution(mcts_dir, "MCTS", fidelity_evolution, gate_count_evolution, mcts_params)
