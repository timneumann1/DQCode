module ExperimentConfig

using ..Types
using ..Helper
using QECCore
using QECCore: Steane7, Shor9, QuantumReedMuller, QuantumTannerGraphProduct, CyclicQuantumTannerGraphProduct, Triangular488, distance
using ..TrivariateBicycleCode

export experiment_configurations

function experiment_configurations()
    ConfigKeyType = String
    ConfigValueType = @NamedTuple{
        code::AbstractCSSCode,
        qpu_sizes::Vector{Int},
        opt_params::OptimisationParameters, 
        genetic_params::GeneticParameters, 
        mcts_params::MCTSParameters, 
        folder::String
    }
    configs = Dict{ConfigKeyType, ConfigValueType}()

    # ------- Steane7 - [4,3] ---------

    configs["steane_4_3"] = (
        Steane7(),
        [4,3],
        OptimisationParameters("jaccard", GateSet([HadamardGate], [CX_Gate])),
        GeneticParameters(1000, 1000, 100, 0.85, 5, 0.5, 1, [2e4, 1, 10, 100], 2),
        MCTSParameters( 5, 1e6, 15, [1e6, 1, 5, 1e3], 0.999, 10),
        joinpath(@__DIR__, "results", string(code_dirname(Steane7())), string([4, 3]))
        )
    
    # ------- Shor9 - [3,3,3] ---------

    configs["shor_3_3_3"] = (
        Shor9(),
        [3,3,3],
        OptimisationParameters("jaccard", GateSet([HadamardGate], [CX_Gate])),
        GeneticParameters(15000, 7500, 100, 0.85, 5, 0.5, 1, [2e4, 1, 10, 100], 2),
        MCTSParameters(6, 1e5, 15, [1e6, 1, 5, 1e2], 0.999, 3.5),
        joinpath(@__DIR__, "results", string(code_dirname(Shor9())), string([3,3,3]))
        )

    # ------- Trivariate Bicycle - [3,3,3,3] ---------

    configs["trivariate_3_3_3_3"] = (
        TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]),
        [3,3,3,3],
        OptimisationParameters("jaccard", GateSet([HadamardGate], [CX_Gate])),
        GeneticParameters(15000, 7500, 100, 0.85, 5, 0.5, 1, [2e4, 1, 10, 100], 2),
        MCTSParameters(2, 5e4, 30, [1e6, 1, 5, 1e4], 0.999, 5.0),
        joinpath(@__DIR__, "results", string(code_dirname(TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]))), string([3,3,3,3]))
        )
    
    # ------- Trivariate Bicycle - [4,4,4] ---------

    configs["trivariate_4_4_4"] = (
        TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]),
        [4,4,4],
        OptimisationParameters("jaccard", GateSet([HadamardGate], [CX_Gate])),
        GeneticParameters(15000, 7500, 100, 0.85, 5, 0.5, 1, [2e4, 1, 10, 100], 2),
        MCTSParameters(2, 5e4, 30, [1e6, 1, 5, 1e4], 0.999, 5.0),
        joinpath(@__DIR__, "results", string(code_dirname(TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]))), string([4,4,4]))
        )

    # ------- Trivariate Bicycle - [6,6] ---------

    configs["trivariate_6_6"] = (
        TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]),
        [6,6],
        OptimisationParameters("jaccard", GateSet([HadamardGate], [CX_Gate])),
        GeneticParameters(15000, 7500, 100, 0.85, 5, 0.5, 1, [2e4, 1, 10, 100], 2),
        MCTSParameters(2, 5e4, 30, [1e6, 1, 5, 1e4], 0.999, 5.0),
        joinpath(@__DIR__, "results", string(code_dirname(TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]))), string([6,6]))
        )

    # ------- Quantum Reed-Muller - [3,3,3,3,3] ---------

    configs["quantum_reed_muller_3_3_3_3_3"] = (
        QuantumReedMuller(4),
        [3,3,3,3,3],
        OptimisationParameters("jaccard", GateSet([HadamardGate], [CX_Gate])),
        GeneticParameters(5000, 3500, 100, 0.85, 5, 0.5, 1, [1.5e4, 1, 10, 100], 2),
        MCTSParameters(3, 5e4, 30, [1e4, 1, 5, 2e2], 0.999, 5.0),
        joinpath(@__DIR__, "results", string(code_dirname(QuantumReedMuller(4))), string([3,3,3,3,3]))
        )
        
    # ------- Quantum Reed-Muller - [5,5,5] ---------

    configs["quantum_reed_muller_5_5_5"] = (
        QuantumReedMuller(4),
        [5,5,5],
        OptimisationParameters("jaccard", GateSet([HadamardGate], [CX_Gate])),
        GeneticParameters(5000, 3500, 100, 0.85, 5, 0.5, 1, [1.5e4, 1, 10, 100], 2),
        MCTSParameters(3, 5e4, 30, [1e4, 1, 5, 2e2], 0.999, 5.0),
        joinpath(@__DIR__, "results", string(code_dirname(QuantumReedMuller(4))), string([5,5,5]))
        )

    # ------- Bivariate Bicycle - [3,3,3,3,3,3] ---------

    configs["bivariate_3_3_3_3_3_3"] = (
        BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]),
        [3,3,3,3,3,3],
        OptimisationParameters("jaccard", GateSet([HadamardGate], [CX_Gate])),
        GeneticParameters(10000, 3500, 100, 0.85, 5, 0.5, 1, [1.5e4, 1, 10, 100], 2),
        MCTSParameters(3, 1e6, 50, [1e5, 1, 5, 1e2], 0.999, 10.0),
        joinpath(@__DIR__, "results", string(code_dirname(BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]))), string([3,3,3,3,3,3]))
        )
    
    # ------- Bivariate Bicycle - [6,6,6] ---------

    configs["bivariate_6_6_6"] = (
        BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]),
        [6,6,6],
        OptimisationParameters("jaccard", GateSet([HadamardGate], [CX_Gate])),
        GeneticParameters(10000, 3500, 100, 0.85, 5, 0.5, 1, [1.5e4, 1, 10, 100], 2),
        MCTSParameters(3, 1e6, 50, [1e5, 1, 5, 1e2], 0.999, 10.0),
        joinpath(@__DIR__, "results", string(code_dirname(BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]))), string([6,6,6]))
        )

    # ------- Bivariate Bicycle - [9,9] ---------

    configs["bivariate_9_9"] = (
        BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]),
        [9,9],
        OptimisationParameters("jaccard", GateSet([HadamardGate], [CX_Gate])),
        GeneticParameters(10000, 3500, 100, 0.85, 5, 0.5, 1, [1.5e4, 1, 10, 100], 2),
        MCTSParameters(3, 1e6, 50, [1e5, 1, 5, 1e2], 0.999, 10.0),
        joinpath(@__DIR__, "results", string(code_dirname(BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]))), string([9,9]))
        )
    return configs

end


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

#distributed_qec_code, type_two_register_sizes, opt_params, genetic_params, mcts_params, gate_set = bivariate_3_3_3_3_3_3()


# init_noise::Float64     # Initialisation noise
#     idle_depolarising_noise::Float64  # idling depolarising probability
#     idle_depolarising_noise_tele::Float64 # idle depolarising probability under telegate
#     single_q_gate_noise::Float64      # single qubit gate noise probability
#     two_q_gate_noise::Float64         # two-qubit gate noise probability
#     measurement_noise::Float64        # Measurement noise
#     two_q_gate_noise_diff_species::Float64         # two-qubit gate noise probability between communication and memory qubit
#     comm_qubit_init_noise::Float64 # Communication qubit init noise, de facto two qubit depolarising noise to mimic the imperfect creation of Bell pairs
#     comm_idle_depolarising_noise::Float64  # idling depolarising probability
#     single_comm_q_gate_noise::Float64 
#     comm_qubit_measurement_noise::Float64
#     classical_comm_noise::Float64


#end