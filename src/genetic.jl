module Genetic

using ..Types
using ..CircuitSimulator
using ..Helper
using ..LogicalEnc

using Random
using QuantumClifford
using QuantumClifford: MixedDestabilizer, @S_str, true_success_stat, false_success_stat, continue_stat, failure_stat, PauliMeasurement, VerifyOp
using BenchmarkTools
using CairoMakie
using Serialization

export genetic_search


function fitness_function(fidelities, circuit_sizes, gen, g)
    penalties = map(cs -> sum(g.fitness_weights[2:4] .* cs) , circuit_sizes) #w[1]*cs[1] + w[2]*cs[2] + w[3]*cs[3]
    return g.fitness_weights[1] .* fidelities .- 2*exp10(-gen/g.num_generations)*penalties  # fitness can decrease over time since weighting is time-dependent if (gen/genetic_params.num_generations)*
end
#(1.2-gen/g.num_generations)
# the prefactor allows for some variability in the first third of the generations, if the parameters are tuned well (keep also some non-1 fidelity in the beginning)
# in the end, we are only optimising while maintaining one fidelity

function initialise_population(code_params, network_specs, genetic_params, gate_set; standard_encoding = false, warm_start = false, warm_start_gates = [], min_len=10)
    
    println("\nInitialising Population... \n")
    population = Vector{CircuitIndividual}(undef, genetic_params.num_individuals) # Vector{Circuit}(undef, num_individuals)
    
    #println("Naive encoding circuit: $( standard_logical_zero_encoding_circuit(qec_code))) of size $(length(standard_logical_zero_encoding_circuit(qec_code)[2]))")
    @assert code_params.qec_code !== nothing 
    #standard_circuit = standard_logical_zero_encoding_circuit(qec_code)
    #standard_circuit_length = length(standard_circuit[2]) 
    #permutation = transpositions_to_perm(reverse(standard_circuit[3]), num_data_qubits)
    #if standard_encoding || warm_start
    #    baseline_gates = baseline_encoding(code_params, network_specs, gate_set, data_storage = false)
    #end
    print(warm_start_gates)
    
    for i in eachindex(population)
        
        gates = Gate[]#Circuit(num_data_qubits, depth)  
        if standard_encoding 
            
            gates = copy(baseline_encoding(code_params, network_specs, gate_set, data_storage = false))
            hadamard_indices = [g.index for g in gates if typeof(g) in gate_set.single_qubit_gates]
            filter!(g -> !(typeof(g) in gate_set.single_qubit_gates), gates)
            for index in hadamard_indices
                pushfirst!(gates, HadamardGate(index))
            end
            
            # TODO: convert to my own type but add indices there!
            
            #permutation = transpositions_to_perm(reverse(standard_circuit[3]), num_data_qubits)
            #circ.gates = cyclic_tanner_encoding_1(circ)#steane_encoding_circuit(circ) # warm start
            
            # for op in standard_circuit[2]
            #     if op isa QuantumClifford.sHadamard
            #         push!(gates, HadamardGate(permutation[op.q]))
            #     elseif op isa QuantumClifford.sX
            #         push!(gates, PauliXGate(permutation[op.q]))
            #     elseif op isa QuantumClifford.sZ
            #         push!(gates, PauliZGate(permutation[op.q]))
            #     elseif op isa QuantumClifford.sZCX
            #         control, target = Tuple(affectedqubits(op))
            #         push!(gates, CNOT_Gate(permutation[control], permutation[target]))
            #     # elseif op isa QuantumClifford.sZCY
            #     #     control, target = Tuple(affectedqubits(op))
            #     #     push!(gates, CNOT_Gate(permutation[control], permutation[target]))
            #     # elseif op isa QuantumClifford.sZCZ
            #     #     control, target = Tuple(affectedqubits(op))
            #     #     push!(gates, CNOT_Gate(permutation[control], permutation[target]))
            #     elseif op isa QuantumClifford.sSWAP
            #         continue
            #     else
            #         error("Unsupported warm-start gate type: $(typeof(op))")
            #     end  
               
            # end
        else
            if warm_start
                # baseline_size = length(baseline_gates)
                # @assert baseline_size > 10
                # start_idx = rand(1:(baseline_size ÷ 2))
                # end_idx = start_idx + (baseline_size - (baseline_size ÷ 2))
                # gates = baseline_gates[start_idx:end_idx]
                gates = copy(warm_start_gates)
                # for op in random_initialisation
                #     if op isa QuantumClifford.sHadamard
                #         push!(gates, HadamardGate(permutation[op.q]))
                #     elseif op isa QuantumClifford.sX
                #         push!(gates, PauliXGate(permutation[op.q]))
                #     elseif op isa QuantumClifford.sZ
                #         push!(gates, PauliZGate(permutation[op.q]))
                #     elseif op isa QuantumClifford.sZCX
                #         control, target = Tuple(affectedqubits(op))
                #         push!(gates, CNOT_Gate(permutation[control], permutation[target]))
                #     # elseif op isa QuantumClifford.sZCY
                #     #     control, target = Tuple(affectedqubits(op))
                #     #     push!(gates, CNOT_Gate(permutation[control], permutation[target]))
                #     # elseif op isa QuantumClifford.sZCZ
                #     #     control, target = Tuple(affectedqubits(op))
                #     #     push!(gates, CNOT_Gate(permutation[control], permutation[target]))
                #     elseif op isa QuantumClifford.sSWAP
                #         continue
                #     else
                #         error("Unsupported warm-start gate type: $(typeof(op))")
                #     end                  
                # end

            else
                ind_length = rand(min_len:2*network_specs.num_data_qubits)
                for _ in 1:ind_length
                    push!(gates,_random_gate(network_specs.num_data_qubits, gate_set))
                end
            end
        end
        population[i] = CircuitIndividual(gates)  # Default constructor of CircuitIndividual

        #println("\n$gates")
    end
    
    println("Individual: $(population[1])")
    return population
end


function evaluate_population(population, code_params, network_specs, opt_params, gate_set)

    fidelities = Vector{Float64}(undef, length(population))
    circuit_sizes = Vector{Vector{Int64} }(undef, length(population))

    #if code_params isa CodeParameters # logical zero state genetic search
    for (idx, circ_individual) in enumerate(population)
        gate_counts = [0,0,0]
        #quantum_clifford_circuit, num_single_qubit_gates, num_two_qubit_gates, num_telegates = construct_executable_circuit(circ_individual.gates, gate_set, network_specs)
        @assert network_specs.num_shots == 1
        circuit = Vector{QuantumClifford.AbstractOperation}()  
        for gate in circ_individual.gates
            T = typeof(gate)
            #println(gate)
            if T in gate_set.single_qubit_gates# isa Union{PauliXGate, PauliYGate, PauliZGate, HadamardGate, SGate} 
                #reward -= mdp.mcts_params.fitness_weights[2]
                gate_counts += [1,0,0]
                #new_quantum_state = execute_circuit([gate_to_apply(T, mdp.network_specs.inv_perm[qubit]) ], initial_quantum_state, num_traj = 1)
                push!(circuit, gate_to_apply(T, gate.index) ) 

            elseif T in gate_set.two_qubit_gates 
                control = gate.control #network_specs.inv_perm[gate.control]
                target = gate.target #network_specs.inv_perm[gate.target]
                control_register = network_specs.register_lookup_array[network_specs.inv_perm[control]] 
                target_register = network_specs.register_lookup_array[network_specs.inv_perm[target]]
                
                if control_register == target_register # the lookup array does not account for the communication qubits
                    #reward -= mdp.mcts_params.fitness_weights[3]
                    gate_counts += [0,1,0]
                else
                    #reward -= mdp.mcts_params.fitness_weights[4]
                    gate_counts += [0,0,1]
                
                end
                push!(circuit, gate_to_apply(T, control, target ))
                
                
                #new_quantum_state = execute_circuit([gate_to_apply(T, mdp.network_specs.inv_perm[control], mdp.network_specs.inv_perm[target]) ], initial_quantum_state, num_traj = 1)
            end
        end


        new_quantum_state = execute_circuit(circuit, network_specs.num_data_qubits, 0; num_traj=1)
        new_quantum_state = only(new_quantum_state)
        new_quantum_state_tab = tab(canonicalize_rref!( stabilizerview(new_quantum_state) )[1])
    #print(new_quantum_state)
    #stab_view = traceout!(copy(stab_view), mdp.network_specs.comm_qubits) # TODO: This can be refactored to ptrace upon stable QS release
    # NOTE: if we swtich to ptrace, then also tableau_distance in the helper.jl needs to be adapted!
    #stab_canon = canonicalize_rref!( new_quantum_state )
    #tableau = tab(stab_canon[1])
    #println("Tableau: $tableau")
        new_quantum_state_bit_matrix = tableau_to_bitmatrix(new_quantum_state_tab) # extract the stabiliser tableau from MixedDestabilizer object
        tab_distance = tableau_distance(new_quantum_state_bit_matrix, code_params.target_bit_matrix, metric = opt_params.tableau_metric)#, mdp.network_specs.data_qubits, mdp.network_specs.comm_qubits, mdp.opt_params.tableau_metric)
        #println(tab_distance)



        # # for perturbative expansion, only the leading order is kept, so probabilities can be smaller than 1, 
        # # also, PauliMeasurement don't work with pert. expansion currently
        # mc_result = execute_circuit(quantum_clifford_circuit, network_specs.num_qubits, network_specs.num_registers; num_traj=network_specs.num_shots)#, keepstates = true) # if specifying num_traj, we use MC sampling, otherwise perturbation.
        # stab_view = stabilizerview(only(mc_result))
        # stab_view = traceout!(copy(stab_view), network_specs.comm_qubits) # TODO: This can be refactored to ptrace upon stable QS release
        # # NOTE: if we swtich to ptrace, then also tableau_distance in the helper.jl needs to be adapted!
        # stab_canon = canonicalize_rref!( stab_view )
        # tableau = tab(stab_canon[1])
        # current_bit_matrix = tableau_to_bitmatrix(tableau) # extract the stabiliser tableau from MixedDestabilizer object
        # tab_distance = tableau_distance(current_bit_matrix, code_params.target_bit_matrix, network_specs.data_qubits, network_specs.comm_qubits, opt_params.tableau_metric)
        fidelities[idx] = 1 - tab_distance # 1 is perfect alignment, here we are in the noiseless setting (one shot)
        circuit_sizes[idx] =  gate_counts#(num_single_qubit_gates, num_two_qubit_gates, num_telegates) #circuit_size(quantum_clifford_circuit) #  length(quantum_clifford_circuit)
        @assert gate_counts[1] == code_params.num_X_checks
        
    end

    # elseif code_params isa CodeParametersLog  # logical CNOT search

    #     for (idx, circ_individual) in enumerate(population)
    #         quantum_clifford_circuit, num_single_qubit_gates, num_two_qubit_gates, num_telegates = construct_executable_circuit(circ_individual.gates, gate_set, network_specs)
    #         tableau_distances = Float64[]
    #         verification_results = Float64[]
    #         # for each of the four logical basis states |00>, |01>, |10> and |11>, we determine the result of applying the circuit to the given initial state
    #         for state in eachindex(code_params.initial_states)
    #             mc_result = execute_circuit(quantum_clifford_circuit, network_specs.num_qubits, network_specs.num_registers, code_params.initial_states[state] ; num_traj=network_specs.num_shots)#, keepstates = true) # if specifying num_traj, we use MC sampling, otherwise perturbation.            
    #             stab_view = stabilizerview(only(mc_result))
    #             stab_view = traceout!(copy(stab_view), network_specs.comm_qubits) # TODO: This can be refactored to ptrace upon stable QS release
    #             # NOTE: if we swtich to ptrace, then also tableau_distance in the helper.jl needs to be adapted!
    #             stab_canon = canonicalize_rref!( stab_view )
    #             tableau = tab(stab_canon[1])
    #             current_bit_matrix = tableau_to_bitmatrix(tableau) # extract the stabiliser tableau from MixedDestabilizer object
    #             #println(tableau)
    #             #println(code_params.target_bit_matrices[idx])
    #             #println(current_bit_matrix)
    #             push!(tableau_distances, tableau_distance(current_bit_matrix, code_params.target_bit_matrices[state], network_specs.data_qubits, network_specs.comm_qubits, opt_params.tableau_metric))  
    #             #verification_logical_state = verify_success(quantum_clifford_circuit, code_params.initial_states[state], code_params.target_states[state], network_specs)  
    #             #push!(verification_results, verification_logical_state )            
    #         end
    #         rand() < 0.00025 ? println("Distance for individual $idx: 00:$(tableau_distances[1]),01:$(tableau_distances[2]), 10:$(tableau_distances[3]), 11:$(tableau_distances[4])" ) : ""
    #         fidelities[idx] = 1 - sum(tableau_distances)/length(tableau_distances) # 1 is perfect alignment, average over four different executions (for four logical basis states)
    #         circuit_sizes[idx] =  (num_single_qubit_gates, num_two_qubit_gates, num_telegates) #circuit_size(quantum_clifford_circuit) #  length(quantum_clifford_circuit)
        
    #     end
    #     println()
    
    # end
        
    return fidelities, circuit_sizes
end

function selection(generation, fitness_scores; tournament_size::Int=5, selection_ratio::Float64=1.0, num_elite = 1)
    length_generation = length(generation) 
    @assert length_generation == length(fitness_scores)
    num_selected = Int(floor(length_generation * selection_ratio))

    #num_elite = max(1, Int(round(num_selected * elite_fraction)))
    elite_idx = sortperm(fitness_scores, rev=true)[1:num_elite]
    
    #elite_idx = argmax(fitness_scores)
    best_individuals = generation[elite_idx]
    
#    best_individuals = Vector{eltype(generation)}()
    remaining = setdiff(collect(eachindex(generation)), elite_idx)

    for _ in 1:(num_selected-num_elite)
        tsize = min(tournament_size, length(remaining))
        tournament = remaining[randperm(length(remaining))[1:tsize]]
        # pick best fitness (max)
        best_idx = tournament[argmax(fitness_scores[tournament])]
        push!(best_individuals, generation[best_idx])
        deleteat!(remaining, findfirst(==(best_idx), remaining))
    end
    return best_individuals
end

function _ensure_min_size!(ind::CircuitIndividual, min_len::Int, num_data_qubits::Int, num_hadamards, gate_set)
    @assert length(ind.gates) >= num_hadamards "Hadamard layer inflicted, something went wrong in the crossover"
    while length(ind.gates) < min_len
        push!(ind.gates, _random_two_qubit_gate(num_data_qubits, gate_set))
    end
    return ind
end

function _cap_individual_size(ind::CircuitIndividual, max_len::Int)
    if length(ind.gates) > max_len
        ind.gates = ind.gates[1:max_len]
    end
    return ind
end

function crossover(best_individuals, mutation_rate, num_data_qubits, max_len, gate_set, num_hadamards)
    new_generation = copy(best_individuals) # TODO: verify there is no shadowing
    parents = copy(best_individuals)
    # could shuffle parents here

    # ADD CONDITION FOR ODD NUMBER: Keep size of generatios constnat 
    i = 1
    #println(length(parents))
    while i < length(parents)
        # first parts of respective partens are intentionally maintained in order to preserve H-CNOT structure
        p1 = parents[i]
        p2 = parents[i+1]

        p1_size = length(p1.gates)
        p2_size = length(p2.gates)
        #print("num hadamard is $num_hadamards")
        #println("p1 is $p1 of size $p1_size,\n p2 is $p2 of size $p2_size")
        cp_1 = rand(num_hadamards+1:p1_size-1)  # crossover point (in vector)
        cp_2 = rand(num_hadamards+1:p2_size-1) 

        ##
        #child = CircuitIndividual(length(parents[i].gates))
        #child.gates = copy(parents[i].gates)
        ##

        child1 = CircuitIndividual(cp_1 + p2_size-cp_2)
        child2 = CircuitIndividual(cp_2 + p1_size-cp_1)

        child1.gates[1:num_hadamards] = p1.gates[1:num_hadamards]
        child1.gates[num_hadamards+1:cp_1] = p1.gates[num_hadamards+1:cp_1]
        child1.gates[cp_1+1:end] = p2.gates[cp_2+1:end]

        child2.gates[1:num_hadamards] = p2.gates[1:num_hadamards]
        child2.gates[num_hadamards+1:cp_2] = p2.gates[num_hadamards+1:cp_2]
        child2.gates[cp_2+1:end] = p1.gates[cp_1+1:end]

        child1 = mutation(child1, mutation_rate, num_data_qubits, gate_set)
        child2 = mutation(child2, mutation_rate, num_data_qubits, gate_set)
        
        child1 = _clean_circuit(child1)
        child2 = _clean_circuit(child2)

        _ensure_min_size!(child1, num_hadamards+2, num_data_qubits, num_hadamards, gate_set)
        _ensure_min_size!(child2, num_hadamards+2, num_data_qubits, num_hadamards, gate_set)

        _cap_individual_size(child1, max_len)
        _cap_individual_size(child2, max_len)

        push!(new_generation, child1, child2)

        i += 2
    end

    if length(best_individuals)%2 != 0
        child = CircuitIndividual(copy(parents[end].gates))
        child = mutation(child, mutation_rate, num_data_qubits, gate_set)
        child = _clean_circuit(child)
        _ensure_min_size!(child, num_hadamards+2, num_data_qubits, num_hadamards, gate_set)
        _cap_individual_size(child, max_len)
        push!(new_generation, child)
    end

    return new_generation
end

function _clean_circuit(circuit)
    # Removes gate duplicates
    new_gates = empty(circuit.gates) # will contain the clean gate sequence
    for gate in circuit.gates
        if !isempty(new_gates) && new_gates[end] == gate
            pop!(new_gates)
        else
            push!(new_gates, gate)
        end
    end
    circuit.gates = new_gates
    return circuit
end

function _random_single_qubit_gate(num_data_qubits, gate_set)
    # choose from 1‑qubit gates you already define
    #gates = (HadamardGate, SGate)# IdentityGate, PauliXGate, PauliYGate, PauliZGate)
    #return gates[rand(1:length(gates))](index)
    gate_type = gate_set.single_qubit_gates[rand(1:length(gate_set.single_qubit_gates))]
    index = rand(1:num_data_qubits)
    return gate_type(index)
end

function _random_two_qubit_gate(num_data_qubits, gate_set)
    gate_type = gate_set.two_qubit_gates[rand(1:length(gate_set.two_qubit_gates))]
    control = rand(1:num_data_qubits)
    target = rand(1:num_data_qubits)
    while target == control
        target = rand(1:num_data_qubits)
    end
    #return CNOT_Gate(control, target)
    return gate_type(control, target)
end

function _random_two_qubit_gate(gate_set, control, target)
    gate_type = gate_set.two_qubit_gates[rand(1:length(gate_set.two_qubit_gates))]    
    return gate_type(control, target)
end

function _random_gate(num_data_qubits, gate_set::GateSet; p_two_qubit=0.7)
    if rand() < p_two_qubit
        return _random_two_qubit_gate(num_data_qubits, gate_set)
    else
        return _random_single_qubit_gate(rand(1:num_data_qubits), gate_set)
    end
end

function mutation(individual, mutation_rate, num_data_qubits, gate_set)
    
    if rand() < mutation_rate
        # Pick one matrix elemenet randomly
        circuit_length = length(individual.gates)
        index = rand(1:circuit_length)
        
        T = typeof(individual.gates[index])
        if T in gate_set.two_qubit_gates# mutate existing CNOT gates
            control = individual.gates[index].control
            target = individual.gates[index].target
            r = rand()
            if r > 0.8 # SWAP control and target
                individual.gates[index] = CX_Gate(target, control)
            elseif r > 0.6
                individual.gates[index] = _random_two_qubit_gate(num_data_qubits, gate_set)
                #individual.gates[index] = _random_single_qubit_gate(rand(1:num_data_qubits), gate_set)
            elseif r > 0.4
                push!(individual.gates, _random_two_qubit_gate(num_data_qubits, gate_set))
            else
                deleteat!(individual.gates, index)
            end


        # elseif T in gate_set.single_qubit_gates
        #     r = rand() 
        #     if r>0.9 # mutate single qubit gates into single-qubit gates on random qubit
        #         individual.gates[index] = _random_single_qubit_gate(rand(1:num_data_qubits), gate_set)
        #     elseif r>0.5 && index != circuit_length
        #         individual.gates[index], individual.gates[index+1] = individual.gates[index+1], individual.gates[index]
        #     elseif r>0.25   # ... or multi-qubit gates
        #         control_index = rand(1:num_data_qubits)
        #         target_index = rand(1:num_data_qubits)
        #         tries = 0
        #         #while ( (ind.gates[target_index,c] isa CNOT_Gate) || (target_index == r) ) && tries < 10
        #         while target_index == control_index
        #             target_index = rand(1:num_data_qubits)
        #             tries +=1
        #         end
        #         if tries >= 10
        #             return individual
        #         end    
        #         individual.gates[index] = _random_two_qubit_gate(gate_set, control_index, target_index)
        #     else
        #         deleteat!(individual.gates, index)
        #     end
        end
    end

    return individual
end


function genetic_search(code_params, network_specs, opt_params, genetic_params, gate_set; warm_start_gates_mcts = [])

    ############## Define environment for GA ##########################

    # TODO: Define noise for the telegates (MOST IMPORTANT)
    # TODO: Outsource this to a function
    #networking_params, genetic_params = define_parameters()                             # retrieve parameters
    # TODO: Mapping stage -> use dictionary to map indices to one another
    # As extracted from Hypergraph Partitoning
    # For partitioning, we should not be using the target state canonical tableau, but the low-weight(!) stabilisers originally obtained from the tableau
    
    
    ##########################
    
    # permutation = data_qubit_partitioning(networking_params, Stabilizer(genetic_params.qec_code))

    #  # For optimisation, we can use the canonical form again
    
    # qec_code_required_qubits = code_n(genetic_params.qec_code)
    # #permutation = collect(1:qec_code_required_qubits)#[1,2,3,4,5,6,7]
    # inv_perm = invperm(permutation)
    # mapping = perm_to_transpositions(deepcopy(permutation)) # careful: without deepcopy, this does in-place substitution of permutation    
    # # NOTE: When generating the infromation for hypergraph part., we need to consult the naive encoding function in the logical encoding script to obtain the logical oeprators.
    # # For the inversion of the circuit, we have a custoim function in circsim.jl since this requries applicaiton of correct indices, accounting for communication qubits.
    # data_qubit_capacities = networking_params.register_sizes
    # @assert sum(data_qubit_capacities) === qec_code_required_qubits
    # num_registers  = length(data_qubit_capacities)
    # register_lookup_array, data_qubits, num_data_qubits = create_lookup_array_cliff(data_qubit_capacities)      # create lookup array

    # num_comm_qubits_per_register = num_registers-1
    # num_qubits = num_data_qubits + num_comm_qubits_per_register*(num_registers) # one verification qubit
    # all_qubits = collect(1:num_qubits)
    # comm_qubits = setdiff(all_qubits, data_qubits)
    
    

    #######################

    ### Product codes
    # H1 = Bool[1 1 1] #Bool[1 0 1 0; 0 1 0 1; 1 1 0 0]
    # H2 = Bool[1 0 0; 1 1 1] #Bool[1 1 0; 0 1 1]
    # #H1 = H2 = parity_matrix(RepCode(3))
    
    # ### (Cyclic) Tanner,  from QS source code: https://github.com/QuantumSavory/QuantumClifford.jl/blob/master/lib/QECCore/src/codes/quantum/quantumtannergraphproduct.jl
    # m = 1
    # #qec_code = QuantumTannerGraphProduct(H1, H2)# Steane7()
    # qec_code = CyclicQuantumTannerGraphProduct(m)

    # Bivariate Bicycle codes, from QS source code: https://github.com/QuantumSavory/QuantumClifford.jl/blob/master/lib/QECCore/src/codes/quantum/generalized_circulant_bivariate_bicycle.jl

    # l, m = 3, 3;
    # A = [(:x, 0), (:x, 1), (:y, 1)];
    # B = [(:y, 0), (:x, 2), (:y, 2)];
    # code = MixedDestabilizer(genetic_params.qec_code)#S"XIXIXIX IXXIIXX IIIXXXX ZIZZIZI ZZIIZZI ZZIZIIZ IZIZIZI"
    # code_stabilizer = stabilizerview(code)
    # print("CODE STABILISERS LOW WEIGHT: $(Stabilizer(genetic_params.qec_code))")
    # logical_z = logicalzview(code)
    # println("Logical operators are $(logical_z)")
    # target_state = vcat(code_stabilizer, logical_z)
    # println("\nTarget state:$target_state and qubit size: $(qec_code_required_qubits) as well as logical qubit size: $(code_k(genetic_params.qec_code))")
    # code_distance = distance(genetic_params.qec_code, DistanceMIPAlgorithm(solver=HiGHS))
    # println("\nCode distance is $code_distance.\n\n")
    
    # instead, can also do MixedDestabiliser(Steane7()) and then extract the stabiliser tableau

    # We want to compare tableaus, so we canonicalize
    # target_canon = canonicalize_rref!(target_state)
    #     #println("target canon: $target_canon")
    #     #println(typeof(target_canon))
    # target_tableau = tab(target_canon[1])
    # target_bit_matrix = tableau_to_bitmatrix(target_tableau)
        #println()
        #println("Target tableau: $target_tableau")
        #println("Bit Matrix: $target_bit_matrix")
    
    #println("Standard encoding circuit: $( standard_logical_zero_encoding_circuit(genetic_params.qec_code))) of size $(length(standard_logical_zero_encoding_circuit(genetic_params.qec_code)[2]))")
    #baseline_raw_circuit, baseline_exec_circuit = baseline_comparison(genetic_params.qec_code, num_data_qubits, networking_params, mapping, inv_perm, register_lookup_array, data_qubits, num_comm_qubits_per_register, num_qubits, data_qubit_capacities)

    #println("Baseline circuits: $baseline_raw_circuit, \n\n\n $baseline_exec_circuit")
    
    ########## Initialise population and run Genetic Algorithm #################
    
    fitness_evolution = Float64[]
    #winner_winner_chicken_dinner = 0 # Place holder for best circuit
    
    # TODO: simplify this init at zero then loop structure
    gen = 0
    @assert !(genetic_params.standard_encoding && genetic_params.warm_start)
    population = initialise_population(code_params, network_specs, genetic_params, gate_set; standard_encoding = genetic_params.standard_encoding, warm_start = genetic_params.warm_start, warm_start_gates = warm_start_gates_mcts)# length(standard_logical_zero_encoding_circuit(genetic_params.qec_code)[2]))
    
    println("############ Initial Population ##########: Generation size: $(length(population))")
    println("")

    # # evaluate population
    # fidelities, circuit_sizes = evaluate_population(population, code_params, network_specs, opt_params)
    # #print("fidelity is $fidelities")
    # fitness_scores = fitness_function(fidelities, circuit_sizes, gen, genetic_params) # TODO: Need to find a fair weighting here
    
    # println("\n Generation $gen of size $(length(population)): Best fitness is $(maximum(fitness_scores)), where fidelity = $(fidelities[argmax(fitness_scores)]) and DQC circuit size = $(circuit_sizes[argmax(fitness_scores)]) \n")
    
    # push!(fitness_evolution, maximum(fitness_scores))
    
    GA_result_raw_circuit = nothing
    GA_result_circuit_sizes = nothing
    while gen<=genetic_params.num_generations

        fidelities, circuit_sizes = evaluate_population(population, code_params, network_specs, opt_params, gate_set)
        fitness_scores = fitness_function(fidelities, circuit_sizes, gen, genetic_params)  # TODO: Need to find a fair weighting here
        if gen % 75 == 0
            println("Generation $gen (/$(genetic_params.num_generations)) of size $(length(population)): Best fitness is $(maximum(fitness_scores)), where fidelity = $(fidelities[argmax(fitness_scores)]) and circuit size:")
            println("Single qubit gates: $(circuit_sizes[argmax(fitness_scores)][1]), Two qubit gates: $(circuit_sizes[argmax(fitness_scores)][2]), Telegates: $(circuit_sizes[argmax(fitness_scores)][3]) \n ")
        end
        push!(fitness_evolution, maximum(fitness_scores))

        if gen == genetic_params.num_generations
            best_individual = argmax(fitness_scores)
            GA_result_raw_circuit = population[best_individual]
            GA_result_circuit_sizes = circuit_sizes[best_individual]
            break
        end
        #@btime evaluate_population($population, $networking_params, $genetic_params, $mapping, $inv_perm, $register_lookup_array, $data_qubits, $comm_qubits, $num_comm_qubits_per_register, $num_qubits, $target_bit_matrix, $data_qubit_capacities, $num_registers)
        # perform selection
        best_individuals = selection(population, fitness_scores; tournament_size = genetic_params.tournament_size, selection_ratio = genetic_params.selection_ratio, num_elite = genetic_params.num_elite)
        # perform crossover (incl. mutations)
        population = crossover(best_individuals, genetic_params.mutation_rate, network_specs.num_data_qubits, genetic_params.max_len, gate_set, code_params.num_X_checks)
        # apply mutations
        #mutated = mutations(new_generation, genetic_params)
        #population = new_generation
        # evaluate population
        gen += 1
        
    end

    # Extract best-performing individual
    #GA_result_raw_circuit_size = circuit_size(gates_to_circuit(GA_result_raw_circuit.gates))
    #print_gate_matrix(winner_winner_chicken_dinner)
    
    #println(winner_winner_chicken_dinner_circuit)
    #println(GA_result_raw_circuit.gates)
    #print(baseline_exec_circuit)
    GA_result_exec_circuit, _,_,_= construct_executable_circuit(GA_result_raw_circuit.gates, gate_set, network_specs, telegate_overhead = true)
    #GA_result_exec_circuit_size = circuit_size(GA_result_exec_circuit)
    #baseline_raw_circuit_size = circuit_size(gates_to_circuit(baseline_raw_circuit.gates))
    #baseline_exec_circuit_size = circuit_size(baseline_exec_circuit)

    println("\n Optimised circuit length (DQC setting):: Single-qubit gates: $(GA_result_circuit_sizes[1]), Two-qubit gates: $(GA_result_circuit_sizes[2]), Telegates: $(GA_result_circuit_sizes[3])") #  vs. $baseline_exec_circuit_size in baseline")

    verification_logical_state = false
    if code_params isa CodeParameters
        verification_logical_state = verify_success(GA_result_exec_circuit, code_params.target_state, network_specs)
    # ^NOTE: this appends a verifyop operation, but we count before so this is irrelevant
        verification_logical_state = verification_logical_state == 1.0 ? true : false
        println("\nVerification Genetic Algorithm successful (target state fidelity; only expressive (binary) in noiseless setting): $verification_logical_state")
    elseif code_params isa CodeParametersLog 
        verification_logical_state1 = verify_success(GA_result_exec_circuit, code_params.initial_states[1], code_params.target_states[1], network_specs)
        println("00 state correct?: $verification_logical_state1")
        verification_logical_state2 = verify_success(GA_result_exec_circuit, code_params.initial_states[2], code_params.target_states[2], network_specs)
        println("01 state correct?: $verification_logical_state2")
        verification_logical_state3 = verify_success(GA_result_exec_circuit, code_params.initial_states[3], code_params.target_states[3], network_specs)
        println("10 state correct?: $verification_logical_state3")
        verification_logical_state4 = verify_success(GA_result_exec_circuit, code_params.initial_states[4], code_params.target_states[4], network_specs)
        println("11 state correct?: $verification_logical_state4")
        if verification_logical_state1 == 1.0 && verification_logical_state2 == 1.0 && verification_logical_state3 == 1.0  && verification_logical_state4 == 1.0 
            verification_logical_state = true
        end
    end

    ###########################################
    ############# DATA STORAGE ################
    ###########################################

    base_ga_dir = joinpath(@__DIR__, "results", string(code_dirname(code_params.qec_code)), "GA")
    GA_dir = next_run_dir(base_ga_dir)

    println("Saving results to $(GA_dir)")
    save_circuit_diagram(GA_result_raw_circuit.gates, gate_set, GA_dir, "GA_raw_circuit__size_$(sum(GA_result_circuit_sizes)).png")
    #save_circuit_diagram(GA_result_exec_circuit, GA_dir, "GA_exec_circuit__size_$(GA_result_circuit_sizes).png")
    
    serialize(joinpath(GA_dir, "raw_circuit.jls"), GA_result_raw_circuit.gates)
    
    open(joinpath(GA_dir, "network_specs.txt"), "w") do io
        println(io, "Network Specifications")
        for fname in fieldnames(Types.NetworkSpecifications)
            println(io, fname, " = ", repr(getfield(network_specs, fname)))
        end
    end
    if code_params isa CodeParameters
        open(joinpath(GA_dir, "code_params.txt"), "w") do io
            println(io, "Code parameters")
            for fname in fieldnames(Types.CodeParameters)
                println(io, fname, " = ", repr(getfield(code_params, fname)))
            end
        end

    elseif code_params isa CodeParametersLog
        open(joinpath(GA_dir, "code_params.txt"), "w") do io
            println(io, "Code parameters")
            for fname in fieldnames(Types.CodeParametersLog)
                println(io, fname, " = ", repr(getfield(code_params, fname)))
            end
        end
    end

    open(joinpath(GA_dir, "genetic_params.txt"), "w") do io
        println(io, "Genetic parameters")
        for fname in fieldnames(Types.GeneticParameters)
            println(io, fname, " = ", repr(getfield(genetic_params, fname)))
        end
    end

    open(joinpath(GA_dir, "GA_raw_circuit.txt"), "w") do io
        println(io, "# Raw gate sequence of size $(sum(GA_result_circuit_sizes))")
        for (i, g) in enumerate(GA_result_raw_circuit.gates)
            println(io, i, "\t", repr(g))
        end
    end

    open(joinpath(GA_dir, "GA_exec_circuit.txt"), "w") do io
        println(io, "# Executable (DQC) circuit operations of size:: $GA_result_circuit_sizes")
        for (i, op) in enumerate(GA_result_exec_circuit)
            println(io, i, "\t", repr(op))
        end
    end

    open(joinpath(GA_dir, "summary.txt"), "w") do io
        println(io, "# Encoding successful: $verification_logical_state")
        println(io, "# Raw gate sequence of size $(sum(GA_result_circuit_sizes))")
        println(io, "# Executable (DQC) circuit operations of size $(GA_result_circuit_sizes) (excl. SWAPS)")
    end

    # Plot the evolution of fitness values
    plot_fitness_evol(GA_dir, fitness_evolution, genetic_params, verification_logical_state)
    


    #=
     # Two methods of verifying the creation of the encoded state (Method 1 is preferable since simpler)

    # 1. VerifyOp 
    push!(circuit, VerifyOp(target_state, data_qubits)) 


    # 2. Apply inverse unitary and measure physical all-zero state
    #    Assumes the existence of an extra qubit for verification
    
    encoding_circuit = naive_encoding_circuit_mapping(code, num_comm_qubits_per_register, register_lookup_array)

    for gate in reverse(encoding_circuit)
        push!(circuit, gate)
    end
   
   
    # push!(circuit, VerifyOp(S"ZIIIIII IZIIIII IIZIIII IIIZIII IIIIZII IIIIIZI IIIIIIZ", data_qubits)) 
   
    # or 

    circuit, pauli_string = measure_zero(circuit, data_qubits, num_qubits) # one verification qubit is appended

    #  The measurment of zero yields true in case the eigenvalue of the state is -1. That is, when the zero state actually lives in the register, 
    #  the measurement outcome will be +1, so the boolean will be False. In that case (when we measured zero), we apply the flip 
    push!(circuit,ConditionalGate(sId1(num_qubits),sX(num_qubits), pauli_string.bit))
    println("circuit after conditional: $circuit")
    # If we have measured zero and thus applied the flip, the qubit at index num_qubits will be in the |1> state, which is stabilise by -Z
    push!(circuit, VerifyOp(S"-Z", [num_qubits])) # 
    =#
        
    ### Circuit Plitting 
    # @with classicalbitslayout => :expanded begin
    #     savecircuit(circuit, "src/plots/circuit_sim/circuit_noise.png")
    # end

    ### Benchmarking

    #@btime tensor_to_circuit($code, $params.depolarising_noise, $params.gate_noise, $circuit_tensor, $mapping, $inv_perm, $register_lookup_array, $data_qubits, $num_comm_qubits_per_register, $num_qubits, $target_state, $data_qubit_capacities)
    #@btime execute_circuit($circuit, $num_qubits, $num_registers, num_traj = $num_traj)

end


function plot_fitness_evol(dir, fitness_evolution, genetic_params, success)
    title_str = "Fitness Evolution : $(genetic_params.num_individuals) individuals over $(genetic_params.num_generations) generations"     
    subtitle_str = "fitness_evolution__success_$(success)"
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="Generation", ylabel="Fitness", title=title_str, subtitle = subtitle_str)
    lines!(ax, 1:length(fitness_evolution), fitness_evolution)

    outpath = joinpath(dir, "Success_$(success).png")
    save(outpath, fig)
end

# function print_gate_matrix(circ::Circuit)
#     show(stdout, "text/plain", circ.gates)
#     println()
# end

end
