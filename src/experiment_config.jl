module ExperimentConfig

using ..Types
using QECCore
using QECCore: Steane7, Shor9, QuantumReedMuller, QuantumTannerGraphProduct, CyclicQuantumTannerGraphProduct, Triangular488, distance
using ..TrivariateBicycleCode

function steane()
    distributed_qec_code = Steane7()
    type_two_register_sizes = [4,3]
    opt_params = OptimisationParameters(
    "hamming" # tableau distance metric
    ) 

    genetic_params = GeneticParameters(
        1000, # individuals
        200, # generations
        100, # max length of (raw) circuit individual
        0.85,  # mutation rate
        5, # tournament size
        0.5, # selection_ratio
        1, # num_elite
        true, #standard_encoding
        false, # warm_start
        [2e4,1,10,100] # fitness weights (fidelity, single-, two- and telegates)
    )

    mcts_params = MCTSParameters(
        10, # max depth
        2500, #iterations
        20, # steps before termination
        [1e4,1,2,3],
        0.85, # discount factor
        1.5 # exploration constant
    )

    gate_set = GateSet(
        [HadamardGate],#HadamardGate],#, SqrtXGate, SGate],#PauliXGate, PauliYGate, PauliZGate, InvSGate, SqrtXGate, InvSqrtXGate],  # single-qubit gates
        [CX_Gate]#CZ_Gate]  # two-qubit gates
    )
    return distributed_qec_code, type_two_register_sizes, opt_params, genetic_params, mcts_params, gate_set
end


function trivariate()
    distributed_qec_code = TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]) # [12,2,3]
    type_two_register_sizes = [3,3,3,3]
    opt_params = OptimisationParameters(
    "hamming" # tableau distance metric
    ) 

    genetic_params = GeneticParameters(
        10000, # individuals
        3500, # generations
        100, # max length of (raw) circuit individual
        0.85,  # mutation rate
        5, # tournament size
        0.5, # selection_ratio
        1, # num_elite
        true, #standard_encoding
        false, # warm_start
        [2e4,1,10,100] # fitness weights (fidelity, single-, two- and telegates)
    )

    mcts_params = MCTSParameters(
        10, # max depth
        50000, #iterations
        30, # steps before termination
        [1e5,1,5,10],
        0.85, # discount factor
        2.5 # exploration constant
    )
#    [1e4,1,5,100],

    gate_set = GateSet(
        [HadamardGate],#HadamardGate],#, SqrtXGate, SGate],#PauliXGate, PauliYGate, PauliZGate, InvSGate, SqrtXGate, InvSqrtXGate],  # single-qubit gates
        [CX_Gate]#CZ_Gate]  # two-qubit gates
    )
    return distributed_qec_code, type_two_register_sizes, opt_params, genetic_params, mcts_params, gate_set
end


function quantum_reed_muller()
    distributed_qec_code = QuantumReedMuller(4) # [15,1,3]
    type_two_register_sizes = [3,3,3,3,3]
    opt_params = OptimisationParameters(
    "hamming" # tableau distance metric
    ) 

    genetic_params = GeneticParameters(
        5000, # individuals
        3500, # generations
        100, # max length of (raw) circuit individual
        0.85,  # mutation rate
        5, # tournament size
        0.5, # selection_ratio
        1, # num_elite
        true, # standard_encoding
        false, # warm_start
        [1.5e4,1,10,100] # fitness weights (fidelity, single-, two- and telegates)
    )

    mcts_params = MCTSParameters(
        12, # max depth
        1000, #iterations
        20, # steps before termination
        [1e4,1,5,100],
        0.85, # discount factor
        2.5 # exploration constant
    )

    gate_set = GateSet(
        [HadamardGate],#HadamardGate],#, SqrtXGate, SGate],#PauliXGate, PauliYGate, PauliZGate, InvSGate, SqrtXGate, InvSqrtXGate],  # single-qubit gates
        [CX_Gate]#CZ_Gate]  # two-qubit gates
    )
    return distributed_qec_code, type_two_register_sizes, opt_params, genetic_params, mcts_params, gate_set
end

#distributed_qec_code = Steane7()
#distributed_qec_code = Shor9()
#distributed_qec_code = TrivariateBicycleViaCirculantMat(2, 2, [(:x, 1), (:y, 0)],[(:x, 1), (:z, 1)]) #[[8,2]] code
# distributed_qec_code = TrivariateBicycleViaCirculantMat(2, 1, [(:x, 1), (:y, 0)],[(:x, 0), (:z, 1)])
#distributed_qec_code = TrivariateBicycleViaCirculantMat(3, 1, [(:x, 2), (:y, 0)],[(:x, 1), (:z, 2)]) 
#distributed_qec_code = BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]) # [18,4,4]
# distributed_qec_code = Triangular488(5)
#distributed_qec_code = BivariateBicycleViaCirculantMat(12, 6, [(:x, 3), (:y, 1), (:y, 2)], [(:y, 3), (:x, 1), (:x, 2)])

#type_two_register_sizes = [4,3]
#type_two_register_sizes = [3,3,3]
#type_two_register_sizes = [4,4]
#type_two_register_sizes = [6]
#type_two_register_sizes = [4]
#type_two_register_sizes = [3,3,3,3]
#type_two_register_sizes = [6,6,6]

##### Optimisation Parameters #####

# opt_params = OptimisationParameters(
#     "hamming" # tableau distance metric
# ) 

# genetic_params = GeneticParameters(
#     2500, # individuals
#     5000, # generations
#     100, # max length of (raw) circuit individual
#     0.85,  # mutation rate
#     5, # tournament size
#     0.5, # selection_ratio
#     1, # num_elite
#     true, # warm_start
#     [1000,1,5,100] # fitness weights (fidelity, single-, two- and telegates)
# )

# mcts_params = MCTSParameters(
#     5, # max depth
#     5000, #iterations
#     15, # steps before termination
#     [5e5,1,5,100],
#     0.85, # discount factor
#     2.5 # exploration constant
# )

# gate_set = GateSet(
#     [HadamardGate],#HadamardGate],#, SqrtXGate, SGate],#PauliXGate, PauliYGate, PauliZGate, InvSGate, SqrtXGate, InvSqrtXGate],  # single-qubit gates
#     [CX_Gate]#CZ_Gate]  # two-qubit gates
# )

# gate_set = GateSet(
#     [HadamardGate, SGate],  # single-qubit gates
#     [CX_Gate]  # two-qubit gates
# )

distributed_qec_code, type_two_register_sizes, opt_params, genetic_params, mcts_params, gate_set = trivariate()

end