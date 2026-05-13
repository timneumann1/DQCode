include(joinpath(@__DIR__, "..", "src", "DQCode.jl"))
using DQCode: create_code_network_data, baseline_encoding_qiskit, baseline_encoding_mqt, circuit_search_gott, circuit_search_MCTS

exp_label = "trivariate_6_6" # available configurations are stored and can be adapted in src/experiment/config.jl

# Create two baseline circuits
baseline_encoding_qiskit(exp_label)
baseline_encoding_mqt(exp_label)



