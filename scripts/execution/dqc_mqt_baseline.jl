include(joinpath(@__DIR__, "..", "..", "src", "DQCode.jl"))
using .DQCode

exp_label = "steane_4_3" # available configurations are stored and can be adapted in src/experiment/config.jl

# Create MQT baseline circuit

# Assumes that the MQT filepath has been set correctly in DQCode, and that PyCall was build correspondingly (see setup)
prep_method = "optimal" # "heuristic"
DQCode.baseline_encoding_mqt(exp_label, DQCode.MQT_PATH, prep_method)

