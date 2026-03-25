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

using QuantumClifford
using QuantumClifford: MixedDestabilizer, Stabilizer, Tableau, stabilizerview, logicalxview, logicalzview, canonicalize_rref!, tab, sX, AbstractOperation
using QuantumClifford.ECC: DistanceMIPAlgorithm

using .Helper: tableau_to_bitmatrix, data_qubit_partitioning, perm_to_transpositions, create_lookup_array, verify_success
using .Types

using .LogicalEnc: baseline_encoding_circuit

using .CircuitSimulator: execute_circuit

using .ExperimentConfig: distributed_qec_code, type_two_register_sizes, opt_params, genetic_params, mcts_params


using HiGHS
using JuMP

function experiment_setup_logical_CNOT(qec_code, register_sizes)
    # Experiment performing a logical CNOT between the first and second logical qubit encoded in the qLDPC code
    
    code = MixedDestabilizer(qec_code)
    code_distance = distance(qec_code, DistanceMIPAlgorithm(solver=HiGHS))

    ###### Quantum Network #####

    @assert sum(register_sizes) == code_n(qec_code)
    permutation = data_qubit_partitioning(register_sizes, Stabilizer(qec_code))
    # Note: Whereas for partitioning, we use the original stabiliser formalism (used in QEC cycles), 
    # for optimisation, we use the canonical form again for commensurability -- the target state is identical since 
    # it is fixed by the entire stabiliser subgroup
    inv_perm = invperm(permutation)
    mapping = perm_to_transpositions(deepcopy(permutation))    
    #num_registers = length(networking_params.register_sizes)
    register_lookup_array, data_qubits, num_data_qubits, num_comm_qubits_per_register = create_lookup_array(register_sizes)
    @assert num_data_qubits == code_n(qec_code)
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

    ###### QEC code #####

    #stabilizers = collect(stabilizerview(code))
    #logical_Zs = collect(logicalzview(code))
    stabilizers = [copy(p) for p in stabilizerview(code)]
    logical_Zs  = [copy(p) for p in logicalzview(code)]
    logical_Xs  = [copy(p) for p in logicalxview(code)]
    @assert length(logical_Zs) == 2 "Currently searching for logical CNOT gate in [[n,2,d]] codes"

    # the inital state needs to be provided to the MC trajectories function, hence it should be a MixedDestabilizer
    # in contrast, the target states only serve for verification, and can thus be Stabilizers
    initial_states = MixedDestabilizer{QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}}}[] 

    target_states = Stabilizer{QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}}}[]
    target_bit_matrices = Matrix{Int}[]

    for b1 in 0:1, b2 in 0:1
        println("\n\n\n\n $b1,$b2 \n\n\n\n")

        _ , baseline_exec_circuit, _ = baseline_encoding_circuit(qec_code, network_specs, logical_Xs = true)
        if b1 == 1
            pushfirst!(baseline_exec_circuit, sX(comm_idx[code_n(qec_code)-1])) #code_n(qec_code)-code_k(qec_code)+1
        end
        if b2 == 1
            pushfirst!(baseline_exec_circuit, sX(comm_idx[code_n(qec_code)])) #code_n(qec_code)-code_k(qec_code)+1
        end
        # for the initial state, we want to encode the physical 00, 01, 10 and 11 into their logical counterparts add the +Z logical for |0>, and the -Z logical for |1>
        #x1_init = (b1 == 0) ? nothing : logical_Xs[1] 
        #x2_init = (b2 == 0) ? nothing : logical_Xs[2]
        #logicalX_init = [x1_init, x2_init]
        #init_generators = copy(stabilizers)
        #push!(init_generators, z1_init)
        #push!(init_generators, z2_init)
        #for i in code_n(qec_code)-code_k(qec_code)+1:code_n(qec_code)
        #_ , baseline_exec_circuit, _ = baseline_encoding_circuit(qec_code, network_specs, logical_Xs = true)
        println("HEEEEEEERE: $baseline_exec_circuit")
        #println(logical_Xs[1])
        # Need to traverse each symbol in logical X and add the respective operation (for CSS code, only X?) to the exec_circuit, with gate index comm_idx(index)
        
        # initial state is a mixed_destabilizer of the corresponding logical state
        mc_result = execute_circuit(baseline_exec_circuit, network_specs.num_qubits, network_specs.num_registers, num_traj=network_specs.num_shots)
        println(mc_result)
        initial_state = only(mc_result)
        # we need to append the communication qubits in the |0> state here
        #initial_state = MixedDestabilizer(Stabilizer(Tableau(init_generators)))

        # for the target state, we want to add the logical +Z/-Z, based on both the target and the control qubit being in states |0> or |1>
        z1_target = (b1 == 0) ? logical_Zs[1] : -logical_Zs[1] 
        z2_target = (b1 == 0) ? ((b2 == 0) ? logical_Zs[2] : -logical_Zs[2] ) : ((b2 == 0) ? -logical_Zs[2] : logical_Zs[2] )
        #target_state = Stabilizer(vcat(stabilizers, [z1_target, z2_target]))
        targ_generators = copy(stabilizers)
        push!(targ_generators, z1_target)
        push!(targ_generators, z2_target)
        target_state = Stabilizer(Tableau(targ_generators))
            
        target_canon = canonicalize_rref!(copy(target_state))  
        target_tableau = tab(target_canon[1])
        target_bit_matrix = tableau_to_bitmatrix(target_tableau)
        push!(initial_states, initial_state)
        push!(target_states, target_state)
        push!(target_bit_matrices, target_bit_matrix)
    end

    print("CODE STABILISERS LOW WEIGHT: \n$(Stabilizer(qec_code))\n")
    println("Logical Z operators are \n$(logicalzview(code))\n")

    test_init = tab(canonicalize_rref!(traceout!(copy(stabilizerview(initial_states[1])), network_specs.comm_qubits))[1])
    test_targ = tab(canonicalize_rref!(copy(target_states[1]))[1])
    println("\nFirst Initial state:$(test_init)\n")
    println("\nFirst target state:$(test_targ)\n")
    # println("\nSecond Initial state:$(initial_states[2])\n")
    # println("\nSecond target state:$(target_states[2])\n")
    # println("\nThird Initial state:$(initial_states[3])\n")
    # println("\nThird target state:$(target_states[3])\n")
    # println("\nFourth Initial state:$(initial_states[4])\n")
    # println("\nFourth target state:$(target_states[4])\n")
    println("\n$(qec_code)[[$(code_n(qec_code)), $(code_k(qec_code)), $code_distance]]-code.\n\n")

    init_tab = tableau_to_bitmatrix(test_init)
    targ_tab = tableau_to_bitmatrix(test_targ)
    test_tab_dist = Helper.tableau_distance(init_tab, targ_tab, network_specs.data_qubits, network_specs.comm_qubits, "hamming")
    println(test_tab_dist)
    success = verify_success(AbstractOperation[], initial_states[1], target_states[1], network_specs) 

    println("success:$success")
    code_params = CodeParametersLog(
        qec_code,
        Stabilizer(qec_code),
        initial_states,
        target_states,
        target_bit_matrices,
        code_n(qec_code),
        code_k(qec_code),
        code_distance
    )

    return code_params, network_specs
end

function experiment_setup_logical_zero_state(qec_code, register_sizes)
    
    ###### QEC code #####
    code = MixedDestabilizer(qec_code)

    target_state = vcat(stabilizerview(code), logicalzview(code))
    target_canon = canonicalize_rref!(copy(target_state))
    target_tableau = tab(target_canon[1])
    target_bit_matrix = tableau_to_bitmatrix(target_tableau)
    code_distance = distance(qec_code, DistanceMIPAlgorithm(solver=HiGHS))

    print("CODE STABILISERS LOW WEIGHT: \n$(Stabilizer(qec_code))\n")
    println("Logical Z operators are \n$(logicalzview(code))\n")
    println("Logical X operators are \n$(logicalxview(code))\n")
    println("\nTarget state:$target_state\n")
    println("\n$(qec_code)[$(code_n(qec_code)), $(code_k(qec_code)), $code_distance]]-code.\n\n")

    code_params = CodeParameters(
        qec_code,
        Stabilizer(qec_code),
        target_state,
        target_bit_matrix,
        code_n(qec_code),
        code_k(qec_code),
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
    code_params, network_specs = experiment_setup_logical_zero_state(distributed_qec_code, type_two_register_sizes)
    return LogicalEnc.baseline_encoding(code_params, network_specs)
end
export run_baseline_encoding

#using .Genetic: run_genetic_search
function run_genetic_search()
    code_params, network_specs = experiment_setup_logical_zero_state(distributed_qec_code, type_two_register_sizes)
    return Genetic.genetic_search(code_params, network_specs, opt_params, genetic_params)
end
export run_genetic_search

function run_genetic_search_logical_CNOT()
    code_params, network_specs = experiment_setup_logical_CNOT(distributed_qec_code, type_two_register_sizes)
    return Genetic.genetic_search(code_params, network_specs, opt_params, genetic_params)
end
export run_genetic_search_logical_CNOT

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