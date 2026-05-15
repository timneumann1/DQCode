include(joinpath(@__DIR__, "..", "..", "src", "DQCode.jl"))
using .DQCode

exp_label = "trivariate_6_6" # available configurations are stored and can be adapted in src/experiment/config.jl

# Set up the code parameters and network specifications
DQCode.create_code_network_data(exp_label)



