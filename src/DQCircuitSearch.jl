module DQCircuitSearch

include("types.jl")
include("trivariate_bicycle_code.jl")
include("experiment_config.jl")
include("helper.jl")
#include("circsim.jl")
include("logical_encoding.jl")
#include("plots.jl")
#include("dtsimulation.jl")
include("genetic.jl")
include("mcts.jl")
include("qec_tools.jl")
include("dqc_state_prep_sim.jl")
#include("deep_q.jl")
#include("parameters.jl")

using QECCore
using QECCore: distance
using .TrivariateBicycleCode

using QuantumClifford
using QuantumClifford: MixedDestabilizer, Stabilizer, Tableau, stabilizerview, logicalxview, logicalzview, canonicalize_rref!, tab, sX, AbstractOperation
using QuantumClifford.ECC: DistanceMIPAlgorithm
using Serialization


using .Helper: tableau_to_bitmatrix, data_qubit_partitioning, perm_to_transpositions, create_lookup_array, verify_success, execute_circuit
using .Types

using .LogicalEnc: baseline_encoding_circuit

#using .CircuitSimulator: execute_circuit

using .ExperimentConfig: distributed_qec_code, type_two_register_sizes, opt_params, genetic_params, mcts_params, gate_set, noise_model, n_shots

using .QECTools
using .DQCLogicalStatePrepSimulator

using HiGHS
using JuMP

function network_setup(qec_code, register_sizes)

    @assert sum(register_sizes) == code_n(qec_code) "$(code_n(qec_code))"
    println("Stabilizers of $qec_code: $(Stabilizer(qec_code))")
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
        n_shots # number of shots
    )
    return network_specs

end

function code_setup(qec_code)

    # setup for encoding circuit search

     ###### QEC code #####
    code = MixedDestabilizer(qec_code)

    target_state = vcat(stabilizerview(code), logicalzview(code))
    target_canon = canonicalize!(copy(target_state))
    num_X_checks = count(any( tableau_to_bitmatrix(tab(target_canon)) .== 1, dims = 2))# count how many rows contain X stabiliser (clean for CSS codes) #count(i -> tableau_to_bitmatrix(tab(target_canon))[i,i]==1, 1:size(target_canon,1))
    println("Number of X checks is $num_X_checks")
    target_canon_rref = canonicalize_rref!(copy(target_state))
    target_tableau = tab(target_canon_rref[1])
    target_bit_matrix = tableau_to_bitmatrix(target_tableau)
    code_distance = distance(qec_code, DistanceMIPAlgorithm(solver=HiGHS))

    print("CODE STABILISERS LOW WEIGHT: \n$(Stabilizer(qec_code))\n")
    println("Logical Z operators are \n$(logicalzview(code))\n")
    println("Logical X operators are \n$(logicalxview(code))\n")
    println("\nTarget state:$target_state\n")
    println("\n$(qec_code)[$(code_n(qec_code)), $(code_k(qec_code)), $code_distance]]-code.\n\n")


    stabilizer_generators = [copy(p) for p in stabilizerview(code)]
    logical_Zs  = logicalzview(code)
    stabilizer_group = [stabilizer_generators[1] * stabilizer_generators[1]] # start with the identity element of the stabilizer group

    for gen in stabilizer_generators
        new_elements = [s * gen for s in stabilizer_group]
        append!(stabilizer_group, new_elements)
    end

    @assert length(Set(stabilizer_group)) == 2^(code_n(qec_code)-code_k(qec_code)) "Stabilizer group does not have the correct size of 2^|stabilizer_generators|"
    
    code_params = CodeParameters(
        qec_code,
        Stabilizer(qec_code),
        stabilizer_group,
        num_X_checks,
        logical_Zs,
        target_state,
        target_bit_matrix,
        code_n(qec_code),
        code_k(qec_code),
        code_distance
    )
    return code_params

end


# function code_setup(qec_code, op )

#     # Experiment performing a logical CNOT between the first and second logical qubit encoded in the qLDPC code
    
#     code = MixedDestabilizer(qec_code)
#     code_distance = distance(qec_code, DistanceMIPAlgorithm(solver=HiGHS))

#     #stabilizers = collect(stabilizerview(code))
#     #logical_Zs = collect(logicalzview(code))
#     stabilizer_generators = [copy(p) for p in stabilizerview(code)]
#     logical_Zs  = [copy(p) for p in logicalzview(code)]
#     logical_Xs  = [copy(p) for p in logicalxview(code)]
#     logical_Ys = [(-im)*logical_Zs[p]*logical_Xs[p] for p in eachindex(logical_Xs)]
#     println("Logical Ys:$logical_Ys")
#     println("\n$(qec_code)[$(code_n(qec_code)), $(code_k(qec_code)), d]]-code.\n\n")
#     @assert length(logical_Zs) == 2 "Currently searching for logical CNOT gate in [[n,2,d]] codes"

#     # the inital state needs to be provided to the MC trajectories function, hence it should be a MixedDestabilizer
#     # in contrast, the target states only serve for verification, and can thus be Stabilizers

#     # build the entire stabiliser group, in general O(2^n)
#     stabilizer_group = [stabilizer_generators[1] * stabilizer_generators[1]] # start with the identity element of the stabilizer group

#     for gen in stabilizer_generators
#         new_elements = [s * gen for s in stabilizer_group]
#         append!(stabilizer_group, new_elements)
#     end

#     @assert length(Set(stabilizer_group)) == 2^(code_n(qec_code)-code_k(qec_code)) "Stabilizer group does not have the correct size of 2^|stabilizer_generators|"
    
#     if op == "H1"
#         # a logical Hadamard maps logical Xs to logical Zs and vice versa
#         target_logical_X1 = [logical_Zs[1]*s for s in stabilizer_group]
#         target_logical_X2 = [logical_Xs[2]*s for s in stabilizer_group]
#         target_logical_Z1 = [logical_Xs[1]*s for s in stabilizer_group]
#         target_logical_Z2 = [logical_Zs[2]*s for s in stabilizer_group]
#     end

   
    # if op == "CX" ### THIS IS CURRENTLY WRONG
    #     for b1 in 0:1, b2 in 0:1
    #         println("\n\n\n\n $b1,$b2 \n\n\n\n")

    #         _ , baseline_exec_circuit, _ = baseline_encoding_circuit(qec_code, network_specs, gate_set, logical_Xs = true)
    #         if b1 == 1
    #             pushfirst!(baseline_exec_circuit, sX(network_specs.comm_idx[code_n(qec_code)-1])) #code_n(qec_code)-code_k(qec_code)+1
    #         end
    #         if b2 == 1
    #             pushfirst!(baseline_exec_circuit, sX(network_specs.comm_idx[code_n(qec_code)])) #code_n(qec_code)-code_k(qec_code)+1
    #         end
    #         # for the initial state, we want to encode the physical 00, 01, 10 and 11 into their logical counterparts add the +Z logical for |0>, and the -Z logical for |1>
    #         #x1_init = (b1 == 0) ? nothing : logical_Xs[1] 
    #         #x2_init = (b2 == 0) ? nothing : logical_Xs[2]
    #         #logicalX_init = [x1_init, x2_init]
    #         #init_generators = copy(stabilizers)
    #         #push!(init_generators, z1_init)
    #         #push!(init_generators, z2_init)
    #         #for i in code_n(qec_code)-code_k(qec_code)+1:code_n(qec_code)
    #         #_ , baseline_exec_circuit, _ = baseline_encoding_circuit(qec_code, network_specs, logical_Xs = true)
    #         println("HEEEEEEERE: $baseline_exec_circuit")
    #         #println(logical_Xs[1])
    #         # Need to traverse each symbol in logical X and add the respective operation (for CSS code, only X?) to the exec_circuit, with gate index comm_idx(index)
            
    #         # initial state is a mixed_destabilizer of the corresponding logical state
    #         mc_result = execute_circuit(baseline_exec_circuit, network_specs.num_qubits, network_specs.num_registers, num_traj=network_specs.num_shots)
    #         println(mc_result)
    #         initial_state = only(mc_result)
    #         # we need to append the communication qubits in the |0> state here
    #         #initial_state = MixedDestabilizer(Stabilizer(Tableau(init_generators)))

    #         # for the target state, we want to add the logical +Z/-Z, based on both the target and the control qubit being in states |0> or |1>
    #         z1_target = (b1 == 0) ? logical_Zs[1] : -logical_Zs[1] 
    #         z2_target = (b1 == 0) ? ((b2 == 0) ? logical_Zs[2] : -logical_Zs[2] ) : ((b2 == 0) ? -logical_Zs[2] : logical_Zs[2] )
            
    #         #target_state = Stabilizer(vcat(stabilizers, [z1_target, z2_target]))
    #         targ_generators = copy(stabilizers)
    #         push!(targ_generators, z1_target)
    #         push!(targ_generators, z2_target)
    #         target_state = Stabilizer(Tableau(targ_generators))
                
    #         target_canon = canonicalize_rref!(copy(target_state))  
    #         target_tableau = tab(target_canon[1])
    #         target_bit_matrix = tableau_to_bitmatrix(target_tableau)
    #         push!(initial_states, initial_state)
    #         push!(target_states, target_state)
    #         push!(target_bit_matrices, target_bit_matrix)
    #     end

    # elseif op == "S"

    #     # TODO: Need to use encoding circuits for the baseline, for which we have XX: ++, so push two Hadamards, ZZ: 00, so don't push anything (see above)
    #     initial_stab1 = copy(stabilizers)
    #     push!(initial_stab1, logical_Xs[1])
    #     push!(initial_stab1, logical_Xs[2])
    #     push!(initial_states, MixedDestabilizer(initial_stab1))

    #     initial_stab2 = copy(stabilizers)
    #     push!(initial_stab2, logical_Zs[1])
    #     push!(initial_stab2, logical_Zs[2])
    #     push!(initial_states, MixedDestabilizer(initial_stab2))


    #     target_stab1 = copy(stabilizers)
    #     push!(target_stab1, logical_Ys[1])
    #     push!(target_stab1, logical_Xs[2])
    #     target_stab1 = Stabilizer(Tableau(target_stab1))
    #     target_canon1 = canonicalize_rref!(copy(target_stab1))  
    #     target_tableau1 = tab(target_canon1[1])
    #     target_bit_matrix1 = tableau_to_bitmatrix(target_tableau1)
    #     push!(target_bit_matrices, target_bit_matrix1)
    #     push!(target_states, target_stab1)

    #     target_stab2 = copy(stabilizers)
    #     push!(target_stab2, logical_Zs[1])
    #     push!(target_stab2, logical_Zs[2])
    #     target_stab2 = Stabilizer(Tableau(target_stab2))
    #     target_canon2 = canonicalize_rref!(copy(target_stab2))  
    #     target_tableau2 = tab(target_canon2[1])
    #     target_bit_matrix2 = tableau_to_bitmatrix(target_tableau2)
    #     push!(target_bit_matrices, target_bit_matrix2)
    #     push!(target_states, target_stab2)
    # end

    # print("CODE STABILISERS LOW WEIGHT: \n$(Stabilizer(qec_code))\n")
    # println("Logical Z operators are \n$(logicalzview(code))\n")

    # test_init = tab(canonicalize_rref!(traceout!(copy(stabilizerview(initial_states[1])), network_specs.comm_qubits))[1])
    # test_targ = tab(canonicalize_rref!(copy(target_states[1]))[1])
    # println("\nFirst Initial state:$(test_init)\n")
    # println("\nFirst target state:$(test_targ)\n")
    # # println("\nSecond Initial state:$(initial_states[2])\n")
    # # println("\nSecond target state:$(target_states[2])\n")
    # # println("\nThird Initial state:$(initial_states[3])\n")
    # # println("\nThird target state:$(target_states[3])\n")
    # # println("\nFourth Initial state:$(initial_states[4])\n")
    # # println("\nFourth target state:$(target_states[4])\n")
    # println("\n$(qec_code)[[$(code_n(qec_code)), $(code_k(qec_code)), $code_distance]]-code.\n\n")

    # init_tab = tableau_to_bitmatrix(test_init)
    # targ_tab = tableau_to_bitmatrix(test_targ)
    # test_tab_dist = Helper.tableau_distance(init_tab, targ_tab, network_specs.data_qubits, network_specs.comm_qubits, "hamming")
    # println(test_tab_dist)
    # success = verify_success(AbstractOperation[], initial_states[1], target_states[1], network_specs) 
    # println("success:$success")
    
    # code_params = CodeParametersLog(
    #     qec_code,
    #     Stabilizer(qec_code),
    #     initial_states,
    #     target_states,
    #     target_bit_matrices,
    #     code_n(qec_code),
    #     code_k(qec_code),
    #     code_distance
    # )

#     code_params = CodeParametersLogical(
#         qec_code,
#         stabilizer_generators,
#         stabilizer_group,
#         logical_Xs,
#         logical_Zs,
#         [target_logical_X1, target_logical_X2],
#         [target_logical_Z1, target_logical_Z2],
#         code_n(qec_code),
#         code_k(qec_code),
#         code_distance
#     )

#     return code_params
# end


function run_baseline_encoding()
    network_specs = network_setup(distributed_qec_code, type_two_register_sizes)
    code_params = code_setup(distributed_qec_code) 
    return LogicalEnc.baseline_encoding(code_params, network_specs, gate_set)
end
export run_baseline_encoding

#using .Genetic: run_genetic_search
function run_genetic_search()
    network_specs = network_setup(distributed_qec_code, type_two_register_sizes)
    code_params = code_setup(distributed_qec_code) 
    return Genetic.genetic_search(code_params, network_specs, opt_params, genetic_params, gate_set)
end
export run_genetic_search

# function run_genetic_search_logical_op()
#     network_specs = network_setup(distributed_qec_code, type_two_register_sizes)
#     code_params = code_setup(distributed_qec_code)#, network_specs, logical_op) 
#     return Genetic.genetic_search(code_params, network_specs, opt_params, genetic_params, gate_set)
# end
# export run_genetic_search_logical_op

#using .MonteCarloTreeSearch: run_MCTS
function run_MCTS()
    network_specs = network_setup(distributed_qec_code, type_two_register_sizes)
    code_params = code_setup(distributed_qec_code)#, logical_op) 
    return MonteCarloTreeSearch.monte_carlo_tree_search(code_params, network_specs, opt_params, mcts_params, gate_set)
end
export run_MCTS

function run_circuit_search()
    network_specs = network_setup(distributed_qec_code, type_two_register_sizes)
    code_params = code_setup(distributed_qec_code)#, logical_op) 
    mcts_gates = MonteCarloTreeSearch.monte_carlo_tree_search(code_params, network_specs, opt_params, mcts_params, gate_set)
    @assert genetic_params.warm_start == true "Set warm start to true in order to use MCTS warm start"
    return Genetic.genetic_search(code_params, network_specs, opt_params, genetic_params, gate_set, warm_start_gates_mcts = mcts_gates)
end
export run_circuit_search

function run_dqc_state_prep()
    network_specs = network_setup(distributed_qec_code, type_two_register_sizes)
    code_params = code_setup(distributed_qec_code)
    filepath = "/Users/tim/Tim/projects/thesis/src/results/TrivariateBicycle/GA/90/raw_circuit.jls"
    gates = deserialize(filepath)
    return DQCLogicalStatePrepSimulator.dqc_state_prep(gates, code_params, network_specs, noise_model)
end
export run_dqc_state_prep

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