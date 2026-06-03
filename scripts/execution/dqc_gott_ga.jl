include(joinpath(@__DIR__, "..", "..", "src", "DQCode.jl"))
using .DQCode

exp_label = "trivariate_4_4_4" # available configurations are stored and can be adapted in src/experiment/config.jl
DQCode.circuit_search_gott_ga(exp_label) # initialise Gottesman Encoding > DQC Compilation > Genetic search pipeline


        


