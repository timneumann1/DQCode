module ExperimentConfig

using ..Types: GeneticParameters, MCTSParameters
using QECCore: Steane7, Shor9, AbstractCSSCode, BivariateBicycleViaCirculantMat#, distance
using ..TrivariateBicycleCode

# This file indicates the code-network configurations we are investigating, with code and gpu_sizes indicating the differnet setups,
# and the parameters being used as defaults for the circuit search (can be changed in scripts)

export experiment_configurations

function experiment_configurations()
    ConfigKeyType = String
    ConfigValueType = @NamedTuple{
        code::AbstractCSSCode,
        qpu_sizes::Vector{Int},
        genetic_params::GeneticParameters, 
        mcts_params::MCTSParameters, 
    }
    configs = Dict{ConfigKeyType, ConfigValueType}()

    # ------- Steane7 - [4,3] ---------

    configs["steane_4_3"] = (
        Steane7(),
        [4,3],
        GeneticParameters(7500, 1500, 100, 0.85, 5, 0.5, 1, [1e4, 1, 10, 1e2], "jaccard"),
        #MCTSParameters(15, [1e6,1,5,2e3], 0.999, "jaccard", false, 5, 5e5, 5.0),
        MCTSParameters(15, [1e6,1, 5, 2.5e4], 0.999, "jaccard", true, 3, 5e4, 10.0), #treereuse  #[1e6, 1, 5, 1e3] #1e6,1,5,1e5], 0.85, "jaccard", 10, 5e4, 1.5),
        )
    
    # ------- Shor9 - [3,3,3] ---------

    configs["shor_3_3_3"] = (
        Shor9(),
        [3,3,3],
        GeneticParameters(7500, 1500, 100, 0.85, 5, 0.5, 1, [1e4, 1, 10, 1e2], "jaccard"),
        MCTSParameters(15, [1e6,1,5,1e2], 0.9, "jaccard", true, 5, 1e5, 3.5),
        #MCTSParameters(15, [1e6,1,5,1e2], 0.999, "jaccard", 4, 5e4, 10.0), #[1e6, 1, 5, 1e2]
        )

    # ------- Trivariate Bicycle - [3,3,3,3] ---------

    configs["trivariate_3_3_3_3"] = (
        TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]),
        [3,3,3,3],
        GeneticParameters(25000, 2500, 100, 0.85, 5, 0.5, 1, [1e4, 1, 10, 5e2], "jaccard"), # changed from 15 to 25000 individuals, and from 1000 to 2500 iterations
        MCTSParameters(30, [1e6,1,5,1e4], 0.99, "jaccard", false, 3, 5e4, 10.0),#, 5, 1e5, 10), 
        #MCTSParameters(30, [1e6,1,5,1e4], 0.99, "jaccard", 2, 5e4, 5.0),#, 5, 1e5, 10), # tree reuse
        #joinpath(@__DIR__, "results", string(code_dirname(TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]))), string([3,3,3,3]))
        )
    
    # # ------- Trivariate Bicycle - [4,4,4] ---------

    configs["trivariate_4_4_4"] = (
        TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]),
        [4,4,4],
        GeneticParameters(25000, 2500, 100, 0.85, 5, 0.5, 1, [1e4, 1, 10, 5e2], "jaccard"), #working version for 596 had 15000 ind and 2500 iter
        #GeneticParameters(15000, 1000, 100, 0.85, 5, 0.5, 1, [1e4, 1, 10, 8e2], "jaccard"),
        MCTSParameters(30, [1e6,1,5,2.5e4], 0.999,"jaccard", false, 3, 5e4, 10.0),
        #MCTSParameters(30, [1e6,1,5,2e4], 0.99,"jaccard", true, 2, 5e4, 5.0),#tree reuse
        #joinpath(@__DIR__, "results", string(code_dirname(TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]))), string([4,4,4]))
        )

    # # ------- Trivariate Bicycle - [6,6] ---------

    configs["trivariate_6_6"] = (
        TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]),
        [6,6],
        GeneticParameters(25000, 2500, 100, 0.85, 5, 0.5, 1, [1e4, 1, 10, 5e2], "jaccard"),
        #MCTSParameters(30, [1e6,1,5,1e4], 0.999,"jaccard", false, 4, 5e4, 10.0),
        MCTSParameters(30, [1e6,1,5,4e4], 0.99,"jaccard", true, 2, 5e4, 5.0), # reuse tree
        #joinpath(@__DIR__, "results", string(code_dirname(TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]))), string([6,6]))
        )

    # ------- Quantum Reed-Muller - [3,3,3,3,3] ---------

    # configs["quantum_reed_muller_3_3_3_3_3"] = (
    #     QuantumReedMuller(4),
    #     [3,3,3,3,3],
    #     GeneticParameters(GateSet([HadamardGate], [CX_Gate]), 15000, 2000, 100, 0.85, 5, 0.5, 1, [1e4, 1, 10, 1e2], "hamming"),
    #     MCTSParameters(GateSet([HadamardGate], [CX_Gate]), 3, 5e5, 30, [1e6, 1, 5, 5e2], 9.999, 7.5, "hamming"),
    #     joinpath(@__DIR__, "results", string(code_dirname(QuantumReedMuller(4))), string([3,3,3,3,3]))
    #     )
        
    # # ------- Quantum Reed-Muller - [5,5,5] ---------

    # configs["quantum_reed_muller_5_5_5"] = (
    #     QuantumReedMuller(4),
    #     [5,5,5],
    #     OptimisationParameters("jaccard", GateSet([HadamardGate], [CX_Gate])),
    #     GeneticParameters(5000, 3500, 100, 0.85, 5, 0.5, 1, [1.5e4, 1, 10, 100], 2),
    #     MCTSParameters(3, 5e4, 40, [1e4, 1, 5, 2e2], 0.999, 5.0),
    #     joinpath(@__DIR__, "results", string(code_dirname(QuantumReedMuller(4))), string([5,5,5]))
    #     )

    # ------- Bivariate Bicycle [[18,4,4]] - [3,3,3,3,3,3] ---------

    configs["bivariate_3_3_3_3_3_3"] = (
        BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]),
        [3,3,3,3,3,3],
        GeneticParameters(25000, 5000, 100, 0.85, 5, 0.5, 1, [5e4, 1, 10, 1e3], "jaccard"), # changed from 15 to 25000 and from 3e4 to 5e4
        MCTSParameters(70, [1e6,1,5,5e3], 0.999, "jaccard", false, 3, 5e4, 10.0), # 3, 5e6, 10 # don't reuse tree
        #MCTSParameters(70, [1e6,1,5,6e3], 0.999, "jaccard", true, 3, 5e4, 10.0),
        #joinpath(@__DIR__, "results", string(code_dirname(BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]))), string([3,3,3,3,3,3]))
        )
    
    # # # ------- Bivariate Bicycle [[18,4,4]] - [6,6,6] ---------

    configs["bivariate_6_6_6"] = (
        BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]),
        [6,6,6],
        GeneticParameters(25000, 5000, 100, 0.85, 5, 0.5, 1, [5e4, 1, 10, 1e3], "jaccard"), # changed from 15 to 25000 and from 4e4 to 5e4
        MCTSParameters(70, [1e6,1,5,2e3], 0.999, "jaccard", false, 3, 5e4, 10.0), # don't reuse tree version
        #MCTSParameters(70, [1e6,1,5,1e4], 0.999, "jaccard", 2, 5e4, 10.0), # reuse tree version
        #joinpath(@__DIR__, "results", string(code_dirname(BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]))), string([6,6,6]))
        )

    # # # ------- Bivariate Bicycle [[18,4,4]] - [9,9] ---------

    configs["bivariate_9_9"] = (
        BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]),
        [9,9],
        GeneticParameters(25000, 5000, 100, 0.85, 5, 0.5, 1, [5e4, 1, 10, 1e3], "jaccard"),
        MCTSParameters(85, [1e6,1,5,1e3], 0.999, "jaccard", false, 3, 5e4, 10.0), # dob;t reuse
        #MCTSParameters(85, [1e6,1,5,6e3], 0.999, "jaccard", 2, 5e4, 10.0), # reuse
        #joinpath(@__DIR__, "results", string(code_dirname(BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]))), string([9,9]))
        )
    
    # # # ------- Bivariate Bicycle [[36,4,6]] - [6,6,6,6,6,6] ---------

    configs["bivariate_6_6_6_6_6_6"] = (
        BivariateBicycleViaCirculantMat(3, 6, [(:x, 1), (:y, 2), (:y, 3)], [(:x, 0), (:y, 1), (:x, 2)]),
        [6,6,6,6,6,6],
        GeneticParameters(5000, 5000, 100, 0.85, 5, 0.5, 1, [5e4, 1, 10, 1e3], "jaccard"),
        MCTSParameters(70, [1e5, 1, 5, 1e2], 0.999, "jaccard", false, 3, 1e6, 10.0),
        #joinpath(@__DIR__, "results", string(code_dirname( BivariateBicycleViaCirculantMat(3, 6, [(:x, 1), (:y, 2), (:y, 3)], [(:x, 0), (:y, 1), (:x, 2)]) )), string([6,6,6,6,6,6]))
        )

    # # ------- Bivariate Bicycle - [[144,12,12]] ---------

    # configs["bivariate_12"] = (
    #     BivariateBicycleViaCirculantMat(12, 6, [(:x, 3), (:y, 1), (:y, 2)], [(:y, 3), (:x, 1), (:x, 2)]),
    #     [12,12,12,12,12,12,12,12,12,12,12,12],
    #     GeneticParameters(15000, 5000, 100, 0.85, 5, 0.5, 1, [5e4, 1, 10, 1e3], "jaccard"),
    #     MCTSParameters(2000, [1e6, 1, 5, 1e2], 0.999, "jaccard", 2, 1e4, 10.0),
    #     joinpath(@__DIR__, "results", string(code_dirname(BivariateBicycleViaCirculantMat(12, 6, [(:x, 3), (:y, 1), (:y, 2)], [(:y, 3), (:x, 1), (:x, 2)]))), string([12,12,12,12,12,12,12,12,12,12,12,12]))
    #     )

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