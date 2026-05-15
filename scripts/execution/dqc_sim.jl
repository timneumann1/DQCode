include(joinpath(@__DIR__, "..", "..", "src", "DQCode.jl"))
using .DQCode

exp_label = "steane_4_3" # available configurations are stored and can be adapted in src/experiment/config.jl
circuit_path = "warmstart_GA/GA_circuit.jls" # set the circuit path (within the {code}/{architecture}/ folder of interest) that you wish to simulate

# Initialise the Monte Carlo Tree Search
DQCode.dqc_simulation(exp_label, DQCode.MQT_PATH, circuit_path)