include(joinpath(@__DIR__, "..", "..", "src", "DQCode.jl"))
using .DQCode

exp_label = "bivariate_6_6_6" # available configurations are stored and can be adapted in src/experiment/config.jl
#circuit_path = "mqt_encoding/mqt_encoding_circuit.jls"
circuit_path = "warmstart_GA/GA_circuit.jls" # set the circuit path (within the {code}/{architecture}/ folder of interest) that you wish to simulate

# Initialise the DQC simulation
method = "heuristic" # "optimal", "heuristic" or "none"

num_samples = 1e5
ps = 10 .^ range(log10(5e-5),log10(1e-3),length=20) 
p_bells = 10 .^ range(log10(1e-3),log10(5e-2), length = 20) 
telegate_idle_depth = 8 # THIS SHOULD STAY FIXED (MOTIVATED BY TELEGATE EXECUTION TIME)

DQCode.dqc_simulation(exp_label, DQCode.MQT_PATH, circuit_path, num_samples, ps, p_bells, telegate_idle_depth, method)