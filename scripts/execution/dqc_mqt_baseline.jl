include(joinpath(@__DIR__, "..", "..", "src", "DQCode.jl"))
using .DQCode

mqt_path = normpath(joinpath(@__DIR__, "..", "..", "..", "qecc")) # assumes that PyCall was build correspondingly (see setup instructions in README.md)
isdir(mqt_path) || error("Cannot find qecc at $mqt_path. Clone qecc next to DQCode, or set mqt_path manually in this script.")

exp_label = "steane_4_3" # available configurations are stored and can be adapted in src/experiment/config.jl
prep_method = "optimal" # "heuristic" or "optimal"
DQCode.baseline_encoding_mqt(exp_label, mqt_path, prep_method) # evaluate MQT-QECC baseline


