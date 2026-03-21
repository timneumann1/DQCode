module Genetic

using ..Types
using ..CircuitSimulator
using ..Helper
using ..LogicalEnc

using Random
using QuantumClifford
using QuantumClifford: MixedDestabilizer, sHadamard, sCNOT, sSWAP, @S_str, true_success_stat, false_success_stat, continue_stat, failure_stat, PauliMeasurement, VerifyOp
using BenchmarkTools
using CairoMakie

export genetic_search


function fitness_function(fidelities, circuit_sizes, gen, genetic_params)
    return 1e7*fidelities-circuit_sizes  # fitness can decrease over time since weighting is time-dependent if (gen/genetic_params.num_generations)*
end

function initialise_population(code_params, network_specs, genetic_params; standard_encoding = false, warm_start = false, min_len=10)
    
    println("\nInitialising Population... \n")
    population = Vector{CircuitIndividual}(undef, genetic_params.num_individuals) # Vector{Circuit}(undef, num_individuals)
    
    #println("Naive encoding circuit: $( standard_logical_zero_encoding_circuit(qec_code))) of size $(length(standard_logical_zero_encoding_circuit(qec_code)[2]))")
    @assert code_params.qec_code !== nothing 
    #standard_circuit = standard_logical_zero_encoding_circuit(qec_code)
    #standard_circuit_length = length(standard_circuit[2]) 
    #permutation = transpositions_to_perm(reverse(standard_circuit[3]), num_data_qubits)
    if standard_encoding || warm_start
        baseline_gates = baseline_encoding(code_params, network_specs, data_storage = false)
    end
    
    for i in eachindex(population)
        
        gates = Gate[]#Circuit(num_data_qubits, depth)  
        if standard_encoding 
            
            gates = baseline_gates
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
                baseline_size = length(baseline_gates)
                @assert baseline_size > 10
                start_idx = rand(1:(baseline_size ÷ 2))
                end_idx = start_idx + (baseline_size - (baseline_size ÷ 2))
                gates = baseline_gates[start_idx:end_idx]
            
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
                    push!(gates,_random_gate(network_specs.num_data_qubits))
                end
            end
        end
        population[i] = CircuitIndividual(gates)  # Default constructor of CircuitIndividual

        #println("\n$gates")
    end
    
    #println("Individual: $(population[1])")
    return population
end


function evaluate_population(population, code_params, network_specs, opt_params)

    fidelities = Vector{Float64}(undef, length(population))
    circuit_sizes = Vector{Float64}(undef, length(population))
    for (idx, circ_individual) in enumerate(population)
        #print(ind_tensor)
        quantum_clifford_circuit = construct_executable_circuit(circ_individual.gates, network_specs)
        #push!(quantum_clifford_circuit, VerifyOp(target_state, data_qubits)) 
       # println(quantum_clifford_circuit)
        #println("Size of circuit:$(length(quantum_clifford_circuit))")

        # for perturbative expansion, only the leading order is kept, so probabilities can be smaller than 1, 
        # also, PauliMeasurement don't work with pert. expansion currently
        mc_result = execute_circuit(quantum_clifford_circuit, network_specs.num_qubits, network_specs.num_registers; num_traj=network_specs.num_shots)#, keepstates = true) # if specifying num_traj, we use MC sampling, otherwise perturbation.
        #print("MC Result:$mc_result")
        # is of type Vector{ QuantumClifford.MixedDestabilizer{ QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}} } }
                
        tableau_distances = Float64[]
        
        for stab in collect(mc_result) # each component here is a MixedDestabilizer

            stab_view = stabilizerview(stab)
            #println("$comm_qubits")
            #println("Stab view: $stab_view")
            #println(typeof(stab_view))
            stab_view = traceout!(copy(stab_view), network_specs.comm_qubits) # TODO: This can be refactored to ptrace upon stable QS release
            # NOTE: if we swtich to ptrace, then also tableau_distance in the helper.jl needs to be adapted!
            #println("Stab view traceout: $stab_view")
            stab_canon = canonicalize_rref!( stab_view )
            tableau = tab(stab_canon[1])
            #println("tableau:$tableau")
            # convert to stabiliser
            #println("current Tableau after canon:$tableau")
            current_bit_matrix = tableau_to_bitmatrix(tableau) # extract the stabiliser tableau from MixedDestabilizer object
            push!(tableau_distances, tableau_distance(current_bit_matrix, code_params.target_bit_matrix, network_specs.data_qubits, network_specs.comm_qubits, opt_params.tableau_metric))
            
            #println(stab_bit_matrix)
            #Determine tableau distance with target_state
        end
        
        #println("\nFinal Steane-7 dict: $(mc_result) \n")
        # if (mc_result[true_success_stat]  + mc_result[false_success_stat]) != genetic_params.num_shots
        #     throw(ErrorException("Some runs were invalid"))
        # end

        #fidelity = (round(mc_result[true_success_stat] / (mc_result[true_success_stat]+mc_result[false_success_stat]),digits=10))
        fidelities[idx] = 1 - sum(tableau_distances)/length(tableau_distances) # 1 is perfect alignment, average over different executions
        #print(quantum_clifford_circuit)
        circuit_sizes[idx] = circuit_size(quantum_clifford_circuit) #  length(quantum_clifford_circuit)
        #println("Hamming distances for individual $idx is in [$(minimum(hamming_distances)),$(maximum(hamming_distances))] ")
        #println("Fitness score for individual $idx is in [$(1-maximum(hamming_distances)),$(1-minimum(hamming_distances))] -> avg. fitness is $(fitness_scores[idx]).  ")

    end
    
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

function _ensure_min_size!(ind::CircuitIndividual, min_len::Int, num_data_qubits::Int)
    while length(ind.gates) < min_len
        push!(ind.gates, _random_gate(num_data_qubits))
    end
    return ind
end

function _cap_individual_size(ind::CircuitIndividual, max_len::Int)
    if length(ind.gates) > max_len
        ind.gates = ind.gates[1:max_len]
    end
    return ind
end

function crossover(best_individuals, mutation_rate, num_data_qubits, max_len)
    new_generation = copy(best_individuals) # TODO: verify there is no shadowing

    parents = copy(best_individuals)
    # could shuffle parents here

    # ADD CONDITION FOR ODD NUMBER: Keep size of generatios constnat 
    i = 1
    while i < length(parents)
        # first parts of respective partens are intentionally maintained in order to preserve H-CNOT structure
        p1 = parents[i]
        p2 = parents[i+1]

        p1_size = length(p1.gates)
        p2_size = length(p2.gates)

        cp_1 = rand(1:p1_size-1)  # crossover point (in vector)
        cp_2 = rand(1:p2_size-1) 

        ##
        #child = CircuitIndividual(length(parents[i].gates))
        #child.gates = copy(parents[i].gates)
        ##

        child1 = CircuitIndividual(cp_1 + p2_size-cp_2)
        child2 = CircuitIndividual(cp_2 + p1_size-cp_1)

        child1.gates[1:cp_1] = p1.gates[1:cp_1]
        child1.gates[cp_1+1:end] = p2.gates[cp_2+1:end]

        child2.gates[1:cp_2] = p2.gates[1:cp_2]
        child2.gates[cp_2+1:end] = p1.gates[cp_1+1:end]

        child1 = mutation(child1, mutation_rate, num_data_qubits)
        child2 = mutation(child2, mutation_rate, num_data_qubits)
        
        _ensure_min_size!(child1, 3, num_data_qubits)
        _ensure_min_size!(child2, 3, num_data_qubits)

        _cap_individual_size(child1, max_len)
        _cap_individual_size(child2, max_len)

        push!(new_generation, child1, child2)

        i += 2
    end

    if length(best_individuals)%2 != 0
        child = CircuitIndividual(copy(parents[end].gates))
        child = mutation(child, mutation_rate, num_data_qubits)
        _ensure_min_size!(child, 3, num_data_qubits)
        _cap_individual_size(child, max_len)
        push!(new_generation, child)
    end

    return new_generation
end

function _random_single_qubit_gate(index)
    # choose from 1‑qubit gates you already define
    gates = (HadamardGate, SGate)# IdentityGate, PauliXGate, PauliYGate, PauliZGate)
    return gates[rand(1:length(gates))](index)
end

function _random_two_qubit_gate(num_data_qubits)
    control = rand(1:num_data_qubits)
    target = rand(1:num_data_qubits)
    while target == control
        target = rand(1:num_data_qubits)
    end
    return CNOT_Gate(control, target)
end

function _random_gate(num_data_qubits; p_two_qubit=0.7)
    if rand() < p_two_qubit
        return _random_two_qubit_gate(num_data_qubits)
    else
        return _random_single_qubit_gate(rand(1:num_data_qubits))
    end
end

function mutation(individual, mutation_rate, num_data_qubits)
    
    if rand() < mutation_rate
        # Pick one matrix elemenet randomly
        circuit_length = length(individual.gates)
        index = rand(1:circuit_length)

        if individual.gates[index] isa CNOT_Gate # mutate existing CNOT gates
            control = individual.gates[index].control
            target = individual.gates[index].target
            r = rand()
            if r > 0.5 # SWAP control and target
                individual.gates[index] = CNOT_Gate(target, control)
            elseif r>0.1
                individual.gates[index] = _random_single_qubit_gate(rand(1:num_data_qubits))
            else
                deleteat!(individual.gates, index)
            end
        else 
            r = rand() 
            if r>0.9 # mutate single qubit gates into single-qubit gates on random qubit
                individual.gates[index] = _random_single_qubit_gate(rand(1:num_data_qubits))
            elseif r>0.5 && index != circuit_length
                individual.gates[index], individual.gates[index+1] = individual.gates[index+1], individual.gates[index]
            elseif r>0.25   # ... or multi-qubit gates
                control_index = rand(1:num_data_qubits)
                target_index = rand(1:num_data_qubits)
                tries = 0
                #while ( (ind.gates[target_index,c] isa CNOT_Gate) || (target_index == r) ) && tries < 10
                while target_index == control_index
                    target_index = rand(1:num_data_qubits)
                    tries +=1
                end
                if tries >= 10
                    return individual
                end    
                individual.gates[index] = CNOT_Gate(control_index, target_index)
            else
                deleteat!(individual.gates, index)
            end
        end
    end

    return individual
end


function genetic_search(code_params, network_specs, opt_params, genetic_params)

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
    population = initialise_population(code_params, network_specs, genetic_params, standard_encoding = false, warm_start=genetic_params.warm_start, min_len = 3)# length(standard_logical_zero_encoding_circuit(genetic_params.qec_code)[2]))
    
    println("############ Generation #$gen ##########: Generation size: $(length(population))")
    println("")

    # evaluate population
    fidelities, circuit_sizes = evaluate_population(population, code_params, network_specs, opt_params)
    #print("fidelity is $fidelities")
    fitness_scores = fitness_function(fidelities, circuit_sizes, gen, genetic_params) # TODO: Need to find a fair weighting here
    
    println("\n Generation $gen of size $(length(population)): Best fitness is $(maximum(fitness_scores)), where fidelity = $(fidelities[argmax(fitness_scores)]) and DQC circuit size = $(circuit_sizes[argmax(fitness_scores)]) \n")
    
    push!(fitness_evolution, maximum(fitness_scores))
    
    while gen<genetic_params.num_generations

        gen += 1
        #@btime evaluate_population($population, $networking_params, $genetic_params, $mapping, $inv_perm, $register_lookup_array, $data_qubits, $comm_qubits, $num_comm_qubits_per_register, $num_qubits, $target_bit_matrix, $data_qubit_capacities, $num_registers)
        # perform selection
        best_individuals = selection(population, fitness_scores, tournament_size = genetic_params.tournament_size, selection_ratio = genetic_params.selection_ratio, num_elite = genetic_params.num_elite)
        # perform crossover (incl. mutations)
        population = crossover(best_individuals, genetic_params.mutation_rate, network_specs.num_data_qubits, genetic_params.max_len)
        # apply mutations
        #mutated = mutations(new_generation, genetic_params)
        #population = new_generation
        # evaluate population
        fidelities, circuit_sizes = evaluate_population(population, code_params, network_specs, opt_params)
        fitness_scores = fitness_function(fidelities, circuit_sizes, gen, genetic_params)  # TODO: Need to find a fair weighting here
        if gen % 25== 0
        println("Generation $gen (/$(genetic_params.num_generations)) of size $(length(population)): Best fitness is $(maximum(fitness_scores)), where fidelity = $(fidelities[argmax(fitness_scores)]) and DQC circuit size = $(circuit_sizes[argmax(fitness_scores)]) \n")
        end
        push!(fitness_evolution, maximum(fitness_scores))
        
    end

    # Extract best-performing individual
    GA_result_raw_circuit = population[argmax(fitness_scores)]
    GA_result_raw_circuit_size = circuit_size(gates_to_circuit(GA_result_raw_circuit.gates))
    #print_gate_matrix(winner_winner_chicken_dinner)
    
    #println(winner_winner_chicken_dinner_circuit)
    #println(GA_result_raw_circuit.gates)
    #print(baseline_exec_circuit)
    GA_result_exec_circuit = construct_executable_circuit(GA_result_raw_circuit.gates, network_specs)
    GA_result_exec_circuit_size = circuit_size(GA_result_exec_circuit)
    #baseline_raw_circuit_size = circuit_size(gates_to_circuit(baseline_raw_circuit.gates))
    #baseline_exec_circuit_size = circuit_size(baseline_exec_circuit)

    println("\n Optimised circuit length (DQC setting): $(GA_result_exec_circuit_size)") #  vs. $baseline_exec_circuit_size in baseline")
    
    verification_logical_state = verify_success(GA_result_exec_circuit, code_params.target_state, network_specs)
    # ^NOTE: this appends a verifyop operation, but we count before so this is irrelevant
    println("\nVerification Genetic Algorithm successful (target state fidelity; only expressive (binary) in noiseless setting): $verification_logical_state")
    verification_logical_state = verification_logical_state == 1.0 ? true : false

    ###########################################
    ############# DATA STORAGE ################
    ###########################################

    base_ga_dir = joinpath(@__DIR__, "results", string(code_dirname(code_params.qec_code)), "GA")
    GA_dir = next_run_dir(base_ga_dir)

    println("Saving results to $(GA_dir)")
    save_circuit_diagram(GA_result_raw_circuit.gates, GA_dir, "GA_raw_circuit__size_$(GA_result_raw_circuit_size).png")
    save_circuit_diagram(GA_result_exec_circuit, GA_dir, "GA_exec_circuit__size_$(GA_result_exec_circuit_size).png")

    open(joinpath(GA_dir, "network_specs.txt"), "w") do io
        println(io, "Network Specifications")
        for fname in fieldnames(Types.NetworkSpecifications)
            println(io, fname, " = ", repr(getfield(network_specs, fname)))
        end
    end

    open(joinpath(GA_dir, "code_params.txt"), "w") do io
        println(io, "Code parameters")
        for fname in fieldnames(Types.CodeParameters)
            println(io, fname, " = ", repr(getfield(code_params, fname)))
        end
    end

    open(joinpath(GA_dir, "genetic_params.txt"), "w") do io
        println(io, "Genetic parameters")
        for fname in fieldnames(Types.GeneticParameters)
            println(io, fname, " = ", repr(getfield(genetic_params, fname)))
        end
    end

    open(joinpath(GA_dir, "GA_raw_circuit.txt"), "w") do io
        println(io, "# Raw gate sequence of size $GA_result_raw_circuit_size")
        for (i, g) in enumerate(GA_result_raw_circuit.gates)
            println(io, i, "\t", repr(g))
        end
    end

    open(joinpath(GA_dir, "GA_exec_circuit.txt"), "w") do io
        println(io, "# Executable (DQC) circuit operations of size $GA_result_exec_circuit_size")
        for (i, op) in enumerate(GA_result_exec_circuit)
            println(io, i, "\t", repr(op))
        end
    end

    open(joinpath(GA_dir, "summary.txt"), "w") do io
        println(io, "# Encoding successful: $verification_logical_state")
        println(io, "# Raw gate sequence of size $(GA_result_raw_circuit_size)")
        println(io, "# Executable (DQC) circuit operations of size $(GA_result_exec_circuit_size) (excl. SWAPS)")
    end

    # Plot the evolution of fitness values
    plot_fitness_evol(GA_dir, fitness_evolution, genetic_params, verification_logical_state)
    
    # @with classicalbitslayout => :expanded begin
    #     try
    #     savecircuit(
    #         gates_to_circuit(GA_result_raw_circuit.gates),
    #         joinpath(results_dir, "GA_raw_circuit__size_$(GA_result_raw_circuit_size).png");
    #         scale = 1
            
    #     )
    #     catch err
    #         @warn "savecircuit failed (circuit likely too large)" err
    #     end
    # end

    # @with classicalbitslayout => :expanded begin
    #     try
    #     savecircuit(
    #         GA_result_exec_circuit,
    #         joinpath(results_dir, "GA_exec_circuit__size$(GA_result_exec_circuit_size)__success_$(verification_logical_state).png");
    #         scale = 1
            
    #     )
    #     catch err
    #         @warn "savecircuit failed (circuit likely too large)" err
    #     end
    # end

    # @with classicalbitslayout => :expanded begin
    #     try
    #     savecircuit(
    #         gates_to_circuit(baseline_raw_circuit.gates),
    #         joinpath(results_dir, "baseline_raw_circuit__size_$(baseline_raw_circuit_size).png");
    #         scale = 1
            
    #     )
    #     catch err
    #         @warn "savecircuit failed (circuit likely too large)" err
    #     end
    # end

    # @with classicalbitslayout => :expanded begin
    #     try
    #     savecircuit(
    #         baseline_exec_circuit,
    #         joinpath(results_dir, "baseline_exec_circuit__size_$(baseline_exec_circuit_size)__success_$(verification_logical_state).png");
    #         scale = 1
            
    #     )
    #     catch err
    #         @warn "savecircuit failed (circuit likely too large)" err
    #     end
    # end

    #^GOOD CODE
    


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
