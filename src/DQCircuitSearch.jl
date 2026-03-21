module DQCircuitSearch

include("types.jl")
include("trivariate_bicycle_code.jl")
include("experiment_config.jl")
include("helper.jl")
include("circsim.jl")
include("logical_encoding.jl")
#include("plots.jl")
#include("dtsimulation.jl")
include("genetic.jl")
#include("mcts.jl")
#include("deep_q.jl")
#include("parameters.jl")

using QECCore
using QECCore: distance
using .TrivariateBicycleCode

using QuantumClifford: MixedDestabilizer, Stabilizer, stabilizerview, logicalzview, canonicalize_rref!, tab
using QuantumClifford.ECC: DistanceMIPAlgorithm

using .Helper: tableau_to_bitmatrix, data_qubit_partitioning, perm_to_transpositions, create_lookup_array
using .Types

using .ExperimentConfig: distributed_qec_code, type_two_register_sizes, opt_params, genetic_params, mcts_params

using HiGHS
using JuMP

# TODO: Retrieve parameters here and build the partitoning etc. centrally before calling optimisers

function experiment_setup(qec_code, register_sizes)
    
    ###### QEC code #####
    code = MixedDestabilizer(qec_code)
    target_state = vcat(stabilizerview(code), logicalzview(code))
    target_canon = canonicalize_rref!(target_state)
    target_tableau = tab(target_canon[1])
    target_bit_matrix = tableau_to_bitmatrix(target_tableau)
    code_distance = distance(qec_code, DistanceMIPAlgorithm(solver=HiGHS))

    print("CODE STABILISERS LOW WEIGHT: \n$(Stabilizer(qec_code))\n")
    println("Logical operators are \n$(logicalzview(code))\n")
    println("\nTarget state:$target_state\n")
    println("\n$(qec_code)[$(code_n(qec_code)), $(code_k(qec_code)), $code_distance]]-code.\n\n")

    code_params = CodeParameters(
        qec_code,
        Stabilizer(qec_code),
        target_state,
        target_bit_matrix,
        code_distance
    )

    ###### Quantum Network #####

    @assert sum(register_sizes) == code_n(code_params.qec_code)
    permutation = data_qubit_partitioning(register_sizes, Stabilizer(code_params.qec_code))
    # Note: Whereas for partitioning, we use the original stabiliser formalism (used in QEC cycles), 
    # for optimisation, we use the canonical form again for commensurability -- the target state is identical since 
    # it is fixed by the entire stabiliser subgroup
    inv_perm = invperm(permutation)
    mapping = perm_to_transpositions(deepcopy(permutation))    
    #num_registers = length(networking_params.register_sizes)
    register_lookup_array, data_qubits, num_data_qubits, num_comm_qubits_per_register = create_lookup_array(register_sizes)
    @assert num_data_qubits == code_n(code_params.qec_code)
    num_qubits = num_data_qubits + num_comm_qubits_per_register * length(register_sizes)
    all_qubits = collect(1:num_qubits)
    comm_qubits = setdiff(all_qubits, data_qubits)
    comm_idx = [ ( index + num_comm_qubits_per_register * (register_lookup_array[index]-1) ) for index in 1:num_data_qubits]
    @assert comm_idx == data_qubits # both arrays are capturing the same mapping
    comm_inv_perm_idx = [ (inv_perm[index] + num_comm_qubits_per_register * (register_lookup_array[inv_perm[index]]-1) ) for index in 1:num_data_qubits]

    println("Number of qubits is $num_qubits, number of comm. qubits per register is $num_comm_qubits_per_register")
    println("Lookup Array: $register_lookup_array")
    println("Data qubits: $data_qubits")

    network_specs = NetworkSpecifications(
        register_sizes, #register sizes of Type-II architecture (here: only fewn memory qubits per core), CircuitSim automatically adds comm. qubits (ancillas are only added in DTS)
        length(register_sizes), # number of registers
        permutation, # permutation for initial mapping
        mapping,     # corresponding transpositions for initial mapping
        inv_perm,    # inverse permutation (to be used in circuit execution)
        register_lookup_array,  # lookup array for core membership 
        data_qubits, # array of data qubit indices
        comm_qubits, # array of communication qubut indices
        num_data_qubits, # number of data qubits
        num_comm_qubits_per_register, # number of communication qubits per register
        num_qubits,     # total number of qubits
        comm_idx,       # array containing indices under communication qubit reindexing
        comm_inv_perm_idx, # array containing indices under communication qubit reindexing after inverse permutation
        0.0, # depolarising noise
        0.0, # gate noise
        0.0, # telegate noise            
        1 # number of shots
    )

    return code_params, network_specs
end


function run_baseline_encoding()
    code_params, network_specs = experiment_setup(distributed_qec_code, type_two_register_sizes)
    return LogicalEnc.baseline_encoding(code_params, network_specs)
end
export run_baseline_encoding

#using .Genetic: run_genetic_search
function run_genetic_search()
    code_params, network_specs = experiment_setup(distributed_qec_code, type_two_register_sizes)
    return Genetic.genetic_search(code_params, network_specs, opt_params, genetic_params)
end
export run_genetic_search

#using .MonteCarloTreeSearch: run_MCTS
# function run_MCTS()
#     code_params, network_specs, opt_params, _ , mcts_params = experiment_setup()
#     return MonteCarloTreeSearch.monte_carlo_tree_search(code_params, network_specs, opt_params, mcts_params)
# end
# export run_MCTS


# using .Parameters: run_parameter_sweep
# export run_parameter_sweep

# using .Circuit_Plots: plot_gate_teleportation
# export plot_gate_teleportation

using .LogicalEnc: run_tests
export run_tests

#using .LogicalEnc: baseline_encoding


# using .MonteCarloTreeSearch: run_MCTS
# export run_MCTS

# using .PPO: run_PPO, run_PPO_curriculum
# export run_PPO, run_PPO_curriculum


end