include(joinpath(@__DIR__, "..", "..", "src", "DQCode.jl"))
using .DQCode

exp_label = "trivariate_4_4_4" # available configurations are stored and can be adapted in src/experiment/config.jl
circuit_path = "mqt_encoding/mqt_encoding_circuit.jls"
#circuit_path = "warmstart_GA/GA_circuit.jls" # set the circuit path (within the {code}/{architecture}/ folder of interest) that you wish to simulate

# Initialise the DQC simulation
method = "heuristic" # "optimal", "heuristic" or "none"

num_samples = 3e5
ps = 10 .^ range(log10(5e-5),log10(1e-3),length=30) 
p_bells = 10 .^ range(log10(1e-3),log10(5e-2), length = 30) 
telegate_idle_depth = 10 # THIS SHOULD STAY FIXED (MOTIVATED BY TELEGATE EXECUTION TIME)
# We assume classical communication to be noiseless and crosstalk noise to be negligible
DQCode.dqc_simulation(exp_label, DQCode.MQT_PATH, circuit_path, num_samples, ps, p_bells, telegate_idle_depth, method)