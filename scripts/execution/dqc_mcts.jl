include(joinpath(@__DIR__, "..", "..", "src", "DQCode.jl"))
using .DQCode

exp_label = "steane_4_3" # available configurations are stored and can be adapted in src/experiment/config.jl
DQCode.circuit_search_mcts(exp_label) # perform Monte Carlo Tree Search (MCTS)

