include(joinpath(@__DIR__, "..", "..", "src", "DQCode.jl"))
using .DQCode

exp_label = "steane_4_3" # available configurations are stored and can be adapted in src/experiment/config.jl

# Create Qiskit baseline circuit
DQCode.baseline_encoding_qiskit(exp_label)



