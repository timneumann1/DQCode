include(joinpath(@__DIR__, "..", "..", "src", "DQCode.jl"))
using .DQCode

exp_label = "shor_3_3_3" # available configurations are stored and can be adapted in src/experiment/config.jl
DQCode.circuit_search_mcts(exp_label) # perform Monte Carlo Tree Search (MCTS)

