include(joinpath(@__DIR__, "..", "src", "DQCode.jl"))
using .DQCode

exp_label = "bivariate_9_9" # available configurations are stored and can be adapted in src/experiment/config.jl

# Initialise the Monte Carlo Tree Search
DQCode.circuit_search_mcts(exp_label)

