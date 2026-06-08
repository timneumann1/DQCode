include(joinpath(@__DIR__, "..", "..", "src", "DQCode.jl"))
using .DQCode

mqt_path = normpath(joinpath(@__DIR__, "..", "..", "..", "qecc")) # assumes that PyCall was build correspondingly (see setup instructions in README.md)
isdir(mqt_path) || error("Cannot find qecc at $mqt_path. Clone qecc next to DQCode, or set mqt_path manually in this script.")

exp_label = "steane_4_3" # available configurations are stored and can be adapted in src/experiment/config.jl

circuit_path = "warmstart_ga/GA_circuit.jls" # set circuit path (within {code}/{architecture}/ folder) pointing to the circit to simulate
method = "optimal" # "optimal", "heuristic" or "none" 

num_samples = 2.5e5
ps = 10 .^ range(log10(5e-4),log10(1e-2),length=35)
p_bells = 10 .^ range(log10(1e-3),log10(5e-2), length = 35)  
telegate_idle_depth = 25
p_single_ratio = 1/100
p_idle_ratio = 1/10
DQCode.dqc_simulation(exp_label, mqt_path, circuit_path, Int(num_samples), ps, p_bells,
                        telegate_idle_depth, p_single_ratio, p_idle_ratio, method) # perform the DQC simulation

