include(joinpath(@__DIR__, "..", "src", "DQCode.jl"))
using .DQCode

exp_label = "shor_3_3_3" # available configurations are stored and can be adapted in src/experiment/config.jl

# Initialise the Gottesman Encoding + GA pipeline

DQCode.circuit_search_gott_ga(exp_label)

@info "Note: The overall best individual might differ from the best-fidelity individual of the last generation, since we are conditioning on fidelity-1.0 circuits here."

        


