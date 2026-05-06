module Genetic

using ..Types
using ..Helper
#using ..EncodingGott

using Random
using QuantumClifford
using QuantumClifford: MixedDestabilizer, @S_str, true_success_stat, false_success_stat, continue_stat, failure_stat, PauliMeasurement, VerifyOp, AbstractOperation
using BenchmarkTools
using Serialization, CSV, DataFrames


export genetic_search


function fitness_function(fidelities, circuit_sizes, gen, g)
    penalties = map(cs -> sum(g.fitness_weights[2:4] .* cs) , circuit_sizes) #w[1]*cs[1] + w[2]*cs[2] + w[3]*cs[3]
    return g.fitness_weights[1] .* fidelities .- exp10(-gen/g.num_generations)*penalties  # fitness can decrease over time since weighting is time-dependent if (gen/genetic_params.num_generations)*
end
#(1.2-gen/g.num_generations)
# the prefactor allows for some variability in the first third of the generations, if the parameters are tuned well (keep also some non-1 fidelity in the beginning)
# in the end, we are only optimising while maintaining one fidelity

function initialise_population(num_individuals, warm_start_gates)#num_hadamards, num_data_qubits, warm_start)
    
    println("\nInitialising Population... \n")
    population = Vector{Vector{AbstractOperation}}(undef, num_individuals) # Vector{Circuit}(undef, num_individuals)
    
    #println("Naive encoding circuit: $( standard_logical_zero_encoding_circuit(qec_code))) of size $(length(standard_logical_zero_encoding_circuit(qec_code)[2]))")
    #@assert code_params.qec_code !== nothing 
    #standard_circuit = standard_logical_zero_encoding_circuit(qec_code)
    #standard_circuit_length = length(standard_circuit[2]) 
    #permutation = transpositions_to_perm(reverse(standard_circuit[3]), num_data_qubits)
    #if standard_encoding || warm_start
    #    baseline_gates = baseline_encoding(code_params, network_specs, gate_set, data_storage = false)
    #end
    #print(warm_start_gates)
    
    for i in eachindex(population)

        #gates = Vector{AbstractOperation}()#Circuit(num_data_qubits, depth)  
        #if warm_start 

        gates = copy(warm_start_gates)
        
        #gates = copy(baseline_encoding(code_params, network_specs, data_storage = false))
        # Make sure that the Hadamards are in the first layer to mitigate conflicts with the genetic search logic (which is a CNOT search)
        # this will not change the circuit since moving H to the beginning doesnt inflict commutation by the construction of baseline and MCTS (not in general of course)
        
        # we assume all SingleQubitGates in the warm start are Hadamard gates, since we are experimenting with CSS codes only
        hadamard_indices = [g.q for g in gates if g isa AbstractSingleQubitOperator]
        filter!(g -> (g isa AbstractTwoQubitOperator), gates) # only keep the TwoQubitGates,and prepend the Hadamrd gates after
        #println(gates)
        #println(hadamard_indices)
        for index in hadamard_indices
            pushfirst!(gates, sHadamard(index))
        end

        # else 
        #     @assert 0 <= num_hadamards <= num_data_qubits "# Hadamards > # Data Qubits"
        #     # choice of hadamards is random, which is acceptable for the raw GA run (hypothesis is that there is not much difference)
        #     hadamard_indices = randperm(num_data_qubits)[1:num_hadamards]
        #     for index in hadamard_indices
        #         push!(gates, sHadamard(index))
        #     end
        #     two_qubit_indices = randperm(num_data_qubits-1)
        #     for idx in two_qubit_indices
        #         push!(gates, sCNOT(idx, idx+1))
        #     end
        # end
        population[i] = gates#CircuitIndividual(gates)  # Default constructor of CircuitIndividual

    end
    
    println("Individual: $(population[1])")
    return population
end


function evaluate_population(population, code_params, network_specs, genetic_params)

    fidelities = Vector{Float64}(undef, length(population))
    circuit_sizes = Vector{Vector{Int64} }(undef, length(population))
    
    #if code_params isa CodeParameters # logical zero state genetic search
    for (idx, circuit) in enumerate(population)
        #circuit, gate_counts = gates_to_circuit(circ_individual.gates, network_specs)
        #circuit = circ_individual
        gcounts = gate_counts(circuit, network_specs)
        @assert gcounts[1] == code_params.num_X_checks
        
        # gate_counts = [0,0,0]
        # #quantum_clifford_circuit, num_single_qubit_gates, num_two_qubit_gates, num_telegates = construct_executable_circuit(circ_individual.gates, gate_set, network_specs)
        # @assert network_specs.num_shots == 1
        # circuit = Vector{QuantumClifford.AbstractOperation}()  
        # for gate in circ_individual.gates
        #     T = typeof(gate)
        #     #println(gate)
        #     if T in gate_set.single_qubit_gates# isa Union{PauliXGate, PauliYGate, PauliZGate, HadamardGate, SGate} 
        #         #reward -= mdp.mcts_params.fitness_weights[2]
        #         gate_counts += [1,0,0]
        #         #new_quantum_state = execute_circuit([gate_to_apply(T, mdp.network_specs.inv_perm[qubit]) ], initial_quantum_state, num_traj = 1)
        #         push!(circuit, gate_to_apply(T, gate.index) ) 

        #     elseif T in gate_set.two_qubit_gates 
        #         control = gate.control #network_specs.inv_perm[gate.control]
        #         target = gate.target #network_specs.inv_perm[gate.target]
        #         control_register = network_specs.register_lookup_array[network_specs.inv_perm[control]] 
        #         target_register = network_specs.register_lookup_array[network_specs.inv_perm[target]]
                
        #         if control_register == target_register # the lookup array does not account for the communication qubits
        #             #reward -= mdp.mcts_params.fitness_weights[3]
        #             gate_counts += [0,1,0]
        #         else
        #             #reward -= mdp.mcts_params.fitness_weights[4]
        #             gate_counts += [0,0,1]
                
        #         end
        #         push!(circuit, gate_to_apply(T, control, target ))
                
                
        #         #new_quantum_state = execute_circuit([gate_to_apply(T, mdp.network_specs.inv_perm[control], mdp.network_specs.inv_perm[target]) ], initial_quantum_state, num_traj = 1)
        #     end
        # end


        new_quantum_state = execute_circuit(circuit, network_specs.num_data_qubits)#, 0; num_traj=1)
        #new_quantum_state = only(new_quantum_state)
        new_quantum_state_tab = tab(canonicalize_rref!( stabilizerview(new_quantum_state) )[1])
    #print(new_quantum_state)
    #stab_view = traceout!(copy(stab_view), mdp.network_specs.comm_qubits) # TODO: This can be refactored to ptrace upon stable QS release
    # NOTE: if we swtich to ptrace, then also tableau_distance in the helper.jl needs to be adapted!
    #stab_canon = canonicalize_rref!( new_quantum_state )
    #tableau = tab(stab_canon[1])
    #println("Tableau: $tableau")
        new_quantum_state_bit_matrix = tableau_to_bitmatrix(new_quantum_state_tab) # extract the stabiliser tableau from MixedDestabilizer object
        tab_distance = tableau_distance(new_quantum_state_bit_matrix, code_params.target_bit_matrix; metric = genetic_params.tableau_metric)#, mdp.network_specs.data_qubits, mdp.network_specs.comm_qubits, mdp.opt_params.tableau_metric)
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
        circuit_sizes[idx] =  gcounts#(num_single_qubit_gates, num_two_qubit_gates, num_telegates) #circuit_size(quantum_clifford_circuit) #  length(quantum_clifford_circuit)
        
        
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

function _ensure_min_size!(ind::Vector{AbstractOperation}, min_len::Int, num_data_qubits::Int, num_hadamards)
    @assert length(ind) >= num_hadamards "Hadamard layer inflicted, something went wrong in the crossover"
    while length(ind) < min_len
        push!(ind, _random_two_qubit_gate(num_data_qubits))
    end
    return ind
end

function _cap_individual_size(ind::Vector{AbstractOperation}, max_len::Int)
    if length(ind) > max_len
        ind = ind[1:max_len]
    end
    return ind
end

function crossover(num_individuals, selected_individuals, mutation_rate, num_data_qubits, max_len, num_hadamards)
    new_generation = Vector{Vector{AbstractOperation}}()#undef, length(selected_individuals)) #copy(best_individuals) # TODO: verify there is no shadowing
    append!(new_generation, copy(selected_individuals)) #keep the selected individuals in the population
    #then make the selected individuals parents of the second half 
    #println("Current length: $(length(new_generation)) (should be 500)")
    parents = selected_individuals
    # could shuffle parents here

    # ADD CONDITION FOR ODD NUMBER: Keep size of generatios constnat 
    i = 1
    #println(length(parents))
    while i < length(parents)
        # first parts of respective partens are intentionally maintained in order to preserve H-CNOT structure
        p1 = parents[i]
        p2 = parents[i+1]

        p1_size = length(p1)
        p2_size = length(p2)
        #print("num hadamard is $num_hadamards")
        #println("p1 is $p1 of size $p1_size,\n p2 is $p2 of size $p2_size")
        # parent generation will always have at least num_hadamards+2 gates (see ensure_min_size), hence the below random choice is valid
        cp_1 = rand(num_hadamards+1:p1_size-1)  # crossover point (in vector)
        cp_2 = rand(num_hadamards+1:p2_size-1) 

        ##
        #child = CircuitIndividual(length(parents[i].gates))
        #child.gates = copy(parents[i].gates)
        ##

        child1 = Vector{AbstractOperation}(undef, cp_1 + p2_size-cp_2) #CircuitIndividual(cp_1 + p2_size-cp_2) # [sId1(1) for _ in 1:num_gates]
        child2 = Vector{AbstractOperation}(undef, cp_2 + p1_size-cp_1) #CircuitIndividual(cp_2 + p1_size-cp_1)

        child1[1:num_hadamards] = p1[1:num_hadamards]
        child1[num_hadamards+1:cp_1] = p1[num_hadamards+1:cp_1]
        child1[cp_1+1:end] = p2[cp_2+1:end]

        child2[1:num_hadamards] = p2[1:num_hadamards]
        child2[num_hadamards+1:cp_2] = p2[num_hadamards+1:cp_2]
        child2[cp_2+1:end] = p1[cp_1+1:end]

        child1 = mutation(child1, mutation_rate, num_data_qubits)
        child2 = mutation(child2, mutation_rate, num_data_qubits)
        
        child1 = _clean_circuit(child1)
        child2 = _clean_circuit(child2)

        _ensure_min_size!(child1, num_hadamards+2, num_data_qubits, num_hadamards)
        _ensure_min_size!(child2, num_hadamards+2, num_data_qubits, num_hadamards)

        _cap_individual_size(child1, max_len)
        _cap_individual_size(child2, max_len)

        push!(new_generation, child1, child2)

        i += 2
    end

    if length(selected_individuals)%2 != 0
        child = copy(parents[end])
        child = mutation(child, mutation_rate, num_data_qubits)
        child = _clean_circuit(child)
        _ensure_min_size!(child, num_hadamards+2, num_data_qubits, num_hadamards)
        _cap_individual_size(child, max_len)
        push!(new_generation, child)
    end
    #println("Current length: $(length(new_generation)) (should be 1000)")

    @assert length(new_generation) == num_individuals

    return new_generation
end

function _clean_circuit(circuit)
    # Removes gate duplicates
    clean_circuit = Vector{AbstractOperation}()#empty(circuit.gates) # will contain the clean gate sequence
    for gate in circuit
        if !isempty(clean_circuit) && clean_circuit[end] == gate
            pop!(clean_circuit)
        else
            push!(clean_circuit, gate)
        end
    end
    #circuit.gates = clean_circuit
    #return circuit
    return clean_circuit
end

# function _random_single_qubit_gate(num_data_qubits, gate_set)
#     # choose from 1‑qubit gates you already define
#     #gates = (HadamardGate, SGate)# IdentityGate, PauliXGate, PauliYGate, PauliZGate)
#     #return gates[rand(1:length(gates))](index)
#     gate_type = gate_set.single_qubit_gates[rand(1:length(gate_set.single_qubit_gates))]
#     index = rand(1:num_data_qubits)
#     return gate_type(index)
# end

function _random_two_qubit_gate(num_data_qubits)
    #gate_type = gate_set.two_qubit_gates[rand(1:length(gate_set.two_qubit_gates))]
    control = rand(1:num_data_qubits)
    target = rand(1:num_data_qubits)
    while target == control
        target = rand(1:num_data_qubits)
    end
    return sCNOT(control, target)
    #return gate_type(control, target)
end

# function _random_two_qubit_gate(gate_set, control, target)
#     gate_type = gate_set.two_qubit_gates[rand(1:length(gate_set.two_qubit_gates))]    
#     return gate_type(control, target)
# end

# function _random_gate(num_data_qubits, gate_set::GateSet; p_two_qubit=0.7)
#     if rand() < p_two_qubit
#         return _random_two_qubit_gate(num_data_qubits, gate_set)
#     else
#         return _random_single_qubit_gate(rand(1:num_data_qubits), gate_set)
#     end
# end

function mutation(individual, mutation_rate, num_data_qubits)
    
    if rand() < mutation_rate
        # Pick one matrix elemenet randomly
        circuit_length = length(individual)
        index = rand(1:circuit_length)
        
        #T = typeof(individual.gates[index])
        #if T in gate_set.two_qubit_gates# mutate existing CNOT gates
        if individual[index] isa AbstractTwoQubitOperator
            control, target = affectedqubits(individual[index])#.control
            #target = individual.gates[index].target
            r = rand()
            if r > 0.8 # SWAP control and target
                individual[index] = sCNOT(target, control)
            elseif r > 0.6
                individual[index] = _random_two_qubit_gate(num_data_qubits)
                #individual.gates[index] = _random_single_qubit_gate(rand(1:num_data_qubits), gate_set)
            elseif r > 0.4
                push!(individual, _random_two_qubit_gate(num_data_qubits))
            else
                deleteat!(individual, index)
            end


        # Don't permute Hadamards in order to maintain some circuit structure (this can be added if results are not as desired)

        # the probability to sample a single qubit gate is inherently lower already, 
        #elseif T in gate_set.single_qubit_gates
        #elseif individual.gates[index] isa SingleQubitGate
            # r = rand() 
            # hadamard_indices = [gate.index for gate in individual.gates[1:num_hadamards]]
            # if r>1 # mutate single qubit gates into single-qubit gates on random qubit
            #     new_idx = rand(1:num_data_qubits)
            #     if new_idx ∉ hadamard_indices
            #         individual.gates[index] = HadamardGate(new_idx)#_random_single_qubit_gate(new_idx, gate_set)
            #     end
            # # elseif r>0.5 && index != circuit_length
            # #     individual.gates[index], individual.gates[index+1] = individual.gates[index+1], individual.gates[index]
            # # elseif r>0.25   # ... or multi-qubit gates
            # #     control_index = rand(1:num_data_qubits)
            # #     target_index = rand(1:num_data_qubits)
            # #     tries = 0
            # #     #while ( (ind.gates[target_index,c] isa CNOT_Gate) || (target_index == r) ) && tries < 10
            # #     while target_index == control_index
            # #         target_index = rand(1:num_data_qubits)
            # #         tries +=1
            # #     end
            # #     if tries >= 10
            # #         return individual
            # #     end    
            # #     individual.gates[index] = _random_two_qubit_gate(gate_set, control_index, target_index)
            # # else
            # #     deleteat!(individual.gates, index)
            # end
        end
    end

    return individual
end


function genetic_search(code_params, network_specs, genetic_params, warm_start_gates, folder)#; warm_start = false, warm_start_gates = [], label = nothing)

    
    ########## Initialise population and run Genetic Algorithm #################
    
    @info "Genetic Search started..."
    # TODO: simplify this init at zero then loop structure
    gen = 0
    #@assert !(genetic_params.standard_encoding && genetic_params.warm_start)
    #if warm_start
    #    @assert !isempty(warm_start_gates) "Please provide warm start gates to use the warm start initialisation to the genetic search"
    #end 

    population = initialise_population(genetic_params.num_individuals, warm_start_gates) #code_params.num_X_checks, network_specs.num_data_qubits, warm_start, warm_start_gates)# length(standard_logical_zero_encoding_circuit(genetic_params.qec_code)[2]))
    
    println("############ Initial Population ##########: Generation size: $(length(population))")
    println("")

    # fidelities, circuit_sizes = evaluate_population(population, code_params, network_specs, opt_params)
    # #print("fidelity is $fidelities")
    # fitness_scores = fitness_function(fidelities, circuit_sizes, gen, genetic_params) # TODO: Need to find a fair weighting here
    
    # println("\n Generation $gen of size $(length(population)): Best fitness is $(maximum(fitness_scores)), where fidelity = $(fidelities[argmax(fitness_scores)]) and DQC circuit size = $(circuit_sizes[argmax(fitness_scores)]) \n")
    
    # push!(fitness_evolution, maximum(fitness_scores))
    
    best_circ_ind = Vector{AbstractOperation}() #CircuitIndividual([])
    best_gcounts = (typemax(Int), typemax(Int), typemax(Int))

    gate_count_evolution = Vector{Vector{Int64}}()
    fidelity_evolution = Float64[]
    fitness_evolution = Float64[]

    while gen<=genetic_params.num_generations

        # ----- Population Evaluation ----------
        fidelities, circuit_sizes = evaluate_population(population, code_params, network_specs, genetic_params)
        fitness_scores = fitness_function(fidelities, circuit_sizes, gen, genetic_params)  # TODO: Need to find a fair weighting here
        
        # ----- Logging and Plotting ----------
        if gen % 75 == 0
            println("Generation $gen (/$(genetic_params.num_generations)) of size $(length(population)): Best fitness is $(maximum(fitness_scores)), where fidelity = $(fidelities[argmax(fitness_scores)]) and circuit size:")
            println("Single qubit gates: $(circuit_sizes[argmax(fitness_scores)][1]), Two qubit gates: $(circuit_sizes[argmax(fitness_scores)][2]), Telegates: $(circuit_sizes[argmax(fitness_scores)][3]) \n ")
        end
        
        # By elite retention, the best individual in the offspring generation will always be the best individual overall, up to ties. 
        # however, this does not mean that it will have fidelity 1 (thus, we only overwrite best_circ_ind and gcounts if fidelity is 1)
        # Since we initialise with a fidelity 1 circuit and have elites, we are guaranteed to have a fidelity 1 circ in the end!
        best_index = argmax(fitness_scores)
        
        best_fidelity = fidelities[best_index]

        if best_fidelity == 1.0
            best_circ_ind = population[best_index]
            best_gcounts = circuit_sizes[best_index]
        end

        # For plitting, we want to see the best individual in terms of fitness, irregardless of fidelity
        push!(gate_count_evolution, circuit_sizes[best_index])
        push!(fidelity_evolution, best_fidelity)
        push!(fitness_evolution, fitness_scores[best_index])

        # ----- Stopping Condition ----------
        # the evolution of fitness values is monotonously (potentially not strictly) increasing
        if gen == genetic_params.num_generations
            break
        end
        #@btime evaluate_population($population, $networking_params, $genetic_params, $mapping, $inv_perm, $register_lookup_array, $data_qubits, $comm_qubits, $num_comm_qubits_per_register, $num_qubits, $target_bit_matrix, $data_qubit_capacities, $num_registers)
        
        # ----- Selection and Crossover ----------

        # perform selection
        selected_individuals = selection(population, fitness_scores; tournament_size = genetic_params.tournament_size, selection_ratio = genetic_params.selection_ratio, num_elite = genetic_params.num_elite)
        # perform crossover (incl. mutations)
        population = crossover(genetic_params.num_individuals, selected_individuals, genetic_params.mutation_rate, network_specs.num_data_qubits, genetic_params.max_len, code_params.num_X_checks)
        # apply mutations
        #mutated = mutations(new_generation, genetic_params)
        #population = new_generation
        # evaluate population
        gen += 1
        
    end

    GA_circ = best_circ_ind
    @info "GA Optimised circuit length in DQC setting:: Single-qubit gates: $(best_gcounts[1]), Two-qubit gates: $(best_gcounts[2]), Telegates: $(best_gcounts[3])" 
    @info "Note: This may differ from the best-fidelity individual of the last generation, since we are looking for fidelity 1 circuits here."
    # ----- Verification ------

    verification_logical_state = verify_success(GA_circ, code_params.target_state, network_specs)
    # ^NOTE: this appends a verifyop operation, but we count before so this is irrelevant
    @info "Verification of Genetic Algorithm circuit successful: $verification_logical_state"
    

    # ----- Data Storage ----------

    #if isnothing(label)
    #    dir = joinpath(folder, "Plain_Genetic_Algorithm")
    #elseif label == "DQC_Compiled_Gottesman"
    dir = joinpath(folder, "dqc_comp_gott>GA")
    # elseif label == "MCTS"
    #     dir = joinpath(folder, "MCTS>GA")
    #end

    mkpath(dir)
    @info "Saving results of GA run to $(dir)"

    serialize(joinpath(dir, "GA_circuit.jls"), GA_circ)

    # if label == "MCTS"
    #     serialize(joinpath(dir, "MCTS_circuit.jls"), warm_start_gates)
    # end

    # open(joinpath(dir, "GA_circuit.txt"), "w") do io
    #     println(io, "# Raw gate sequence of size $(sum(best_gcounts))")
    #     for (i, g) in enumerate(GA_circ)
    #         println(io, i, "\t", repr(g))
    #     end
    # end

    save_circuit_diagram(GA_circ, dir, "GA_circuit.png")

    open(joinpath(dir, "genetic_params.txt"), "w") do io
        println(io, "Genetic parameters")
        for fname in fieldnames(Types.GeneticParameters)
            println(io, fname, " = ", repr(getfield(genetic_params, fname)))
        end
        # if label == "DQC_Compiled_Gottesman"
        #     println(io, "Warm-Start circuit can be accessed in Gottesman encoding folder" )
        # end
    end

    open(joinpath(dir, "summary.txt"), "w") do io
        println(io, "# Encoding successful: $verification_logical_state")
        println(io, "# Raw gate sequence of size $(sum(best_gcounts))")
        println(io, "# Executable DQC circuit with $(best_gcounts[1]) single qubit gates, $(best_gcounts[2]) two qubit gates and $(best_gcounts[3]) telegates")
    end
    
    #df = DataFrame(Evolutions = ["fitness_evolution", "fidelity_evolution", "gate_count_evolution"], Values = [fitness_evolution, fidelity_evolution, gate_count_evolution])

    df = DataFrame(
        fitness_evolution   = fitness_evolution,
        fidelity_evolution  = fidelity_evolution,
        single_count = [gc[1] for gc in gate_count_evolution],
        two_qubit_count = [gc[2] for gc in gate_count_evolution],
        telegate_count = [gc[3] for gc in gate_count_evolution]
    )
    CSV.write(joinpath(dir, "genetic_evolution.csv"), df)
    #optimiser_label = "Warm-Start Genetic Algorithm" # isnothing(label) ? "Genetic Algorithm" : 

    plot_evolution(dir, "Warm-Start Genetic Algorithm", fitness_evolution, fidelity_evolution, gate_count_evolution, genetic_params, verification_logical_state)

    return verification_logical_state, best_gcounts
end


end
