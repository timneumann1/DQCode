module Genetic

using ..Types
using ..CircuitSimulator
using ..Helper
using ..LogicalEnc

using Random
using Quantikz: savecircuit, @with, classicalbitslayout
using QECCore
using QECCore: Steane7, QuantumTannerGraphProduct, CyclicQuantumTannerGraphProduct
using QuantumClifford.ECC: DistanceMIPAlgorithm
using HiGHS
using QuantumClifford
using QuantumClifford: MixedDestabilizer, sHadamard, sCNOT, sSWAP, @S_str, true_success_stat, false_success_stat, continue_stat, failure_stat, PauliMeasurement, VerifyOp
using BenchmarkTools
using CairoMakie

export run_genetic_search

function define_parameters()

    #TODO: Rename or introduce a SimulationParameters type for this sim as well (in addition to DTS)
    networking_params = NetworkingParameters(
        [18], #register sizes of Type-I architecture (here: only memory qubits per core), CircuitSim automatically adds comm. qubits (ancillas are only added in DTS)
        0.0, # depolarising_prob 
        0.0, # gate_noise_prob 
        0.0, # Telegate noise (depolarising channel)
    )

    genetic_params = GeneticParameters(
        500, # individuals
        1000, # generations
        1, # shots
        1,  # mutation rate
        5, # tournament size
        0.5, # selection_ratio
        2, # num_elite
        true, # warm_start
        BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]) # qec code
        )
    return networking_params, genetic_params
end


function initialise_population(num_individuals, num_data_qubits; warm_start = false, qec_code = nothing)
    population = Vector{CircuitIndividual}(undef, num_individuals) # Vector{Circuit}(undef, num_individuals)
    
    #println("Naive encoding circuit: $( standard_logical_zero_encoding_circuit(qec_code))) of size $(length(standard_logical_zero_encoding_circuit(qec_code)[2]))")

    for i in eachindex(population)
        
        gates = Gate[]#Circuit(num_data_qubits, depth)  
        if warm_start 
            @assert qec_code !== nothing  
            # TODO: convert to my own type but add indices there!
            
            permutation = transpositions_to_perm(reverse(standard_logical_zero_encoding_circuit(qec_code)[3]), num_data_qubits)
            #circ.gates = cyclic_tanner_encoding_1(circ)#steane_encoding_circuit(circ) # warm start
            for op in standard_logical_zero_encoding_circuit(qec_code)[2]
                if op isa QuantumClifford.sHadamard
                    push!(gates, HadamardGate(permutation[op.q]))
                elseif op isa QuantumClifford.sX
                    push!(gates, PauliXGate(permutation[op.q]))
                elseif op isa QuantumClifford.sZ
                    push!(gates, PauliZGate(permutation[op.q]))
                elseif op isa QuantumClifford.sZCX
                    control, target = Tuple(affectedqubits(op))
                    push!(gates, CNOT_Gate(permutation[control], permutation[target]))
                # elseif op isa QuantumClifford.sZCY
                #     control, target = Tuple(affectedqubits(op))
                #     push!(gates, CNOT_Gate(permutation[control], permutation[target]))
                # elseif op isa QuantumClifford.sZCZ
                #     control, target = Tuple(affectedqubits(op))
                #     push!(gates, CNOT_Gate(permutation[control], permutation[target]))
                elseif op isa QuantumClifford.sSWAP
                    continue
                else
                    error("Unsupported warm-start gate type: $(typeof(op))")
                end  
               
            end
        end
        #println("\n\n\n\n\n\n")
        #print(gates)
        population[i] = CircuitIndividual(gates)
    end
    
    println("Indiviudal: $(population[1])")
    return population
end

function tableau_to_bitmatrix(tableau::QuantumClifford.Tableau{<:AbstractVector{UInt8}, <:AbstractMatrix{<:Unsigned}})
    rows, cols = size(tableau)
    bits = Matrix{Int}(undef, rows, cols+1)#  falses(rows, cols)
    @inbounds for r in 1:rows
        for c in 1:cols
            x, z = tableau[r, c]      # every entry of the Tableau contains a tuple (x,z): (0,0) is I, (1,0) is X, (0,1) is Z, (1,1) is Y
            bits[r, c] = x + 2*z # we map I, X, Z, Y to 0, 1, 2 and 3 to later determine the Hamming distance (in how many entries they disagree)
        end
        # phases
        if (tableau.phases[r] ∉ (0, 2))
            throw("Phase of the tableau is imaginary. Please investigate this case.")
        end
        bits[r, cols + 1] = 1/2*tableau.phases[r]  # phase +1 is represented as 0, phase -1 is represented as 2
        # -> positive phase is represented as 0, negative phase as 1
        #println(bits[r,:])
    end
    return bits
    
end

function hamming_distance(matrix::Matrix{Int}, target_matrix::Matrix{Int}, data_qubits, comm_qubits)
    #check that both marices have same dimensions
   
    # we want to compare the data qubits, so we only keep the columns (qubits) corresponding to data qubits AND the phase column
    cols_keep = vcat(data_qubits, size(matrix, 2))
    matrix = matrix[:, cols_keep] 
    # we also want to eliminate the first #comm_qubits rows, which contain the stabilisers of the communication qubits 
    # the assumption is that the tableau is always in a product state of dataqubits and comm qubits, which is valid 
    # since the comm qubits get measured and then reset to zero
    matrix = matrix[(length(comm_qubits)+1):end,:]
    
    # the target matrix already has the lexicographical ordering of data qubits, so no need to filter or change here
    @assert size(matrix) == size(target_matrix)
    # hamming count = 0
    # for rows
    #     for columns
    #         if numbers at repective positions different
    #             increae hamming count by one
    # return hammingcount/(rows*cols)
    #println("Final matrix: $matrix")
    #println("Final target: $target_matrix")

    return count(!iszero, matrix .!= target_matrix) / length(matrix)

end

function circuit_size(circuit)
    return count(op ->
        !(op isa QuantumClifford.NoiseOp) &&
        !(op isa QuantumClifford.sSWAP) &&
        !(op isa QuantumClifford.sId1),
        circuit
    )
end

function evaluate_population(population, networking_params, genetic_params, mapping, inv_perm, register_lookup_array, data_qubits, comm_qubits, num_comm_qubits_per_register, num_qubits, target_bit_matrix, data_qubit_capacities, num_registers)
    fidelities = Vector{Float64}(undef, length(population))
    circuit_sizes = Vector{Float64}(undef, length(population))
    for (idx, circ_individual) in enumerate(population)
        #print(ind_tensor)
        quantum_clifford_circuit = construct_executable_circuit(networking_params.depolarising_noise, networking_params.gate_noise, networking_params.telegate_noise, circ_individual.gates, mapping, inv_perm, register_lookup_array, data_qubits, num_comm_qubits_per_register, num_qubits, data_qubit_capacities)
        #push!(quantum_clifford_circuit, VerifyOp(target_state, data_qubits)) 
       # println(quantum_clifford_circuit)
        #println("Size of circuit:$(length(quantum_clifford_circuit))")

        # for perturbative expansion, only the leading order is kept, so probabilities can be smaller than 1, 
        # also, PauliMeasurement don't work with pert. expansion currently
        mc_result = execute_circuit(quantum_clifford_circuit, num_qubits, num_registers; num_traj=genetic_params.num_shots)#, keepstates = true) # if specifying num_traj, we use MC sampling, otherwise perturbation.
        #print("MC Result:$mc_result")
        # is of type Vector{ QuantumClifford.MixedDestabilizer{ QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}} } }
                
        hamming_distances = Float64[]
        
        for stab in collect(mc_result) # each component here is a MixedDestabilizer

            stab_view = stabilizerview(stab)
            #println("$comm_qubits")
            #println("Stab view: $stab_view")
            #println(typeof(stab_view))
            stab_view = traceout!(copy(stab_view), comm_qubits) # TODO: This can be refactored to ptrace upon stable QS release
            #println("Stab view traceout: $stab_view")
            stab_canon = canonicalize_rref!( stab_view )
            tableau = tab(stab_canon[1])
            #println("tableau:$tableau")
            # convert to stabiliser
            #println("current Tableau after canon:$tableau")
            current_bit_matrix = tableau_to_bitmatrix(tableau) # extract the stabiliser tableau from MixedDestabilizer object
            push!(hamming_distances, hamming_distance(current_bit_matrix, target_bit_matrix, data_qubits, comm_qubits))
            
            #println(stab_bit_matrix)
            #Determine tableau distance with target_state
        end
        
        #println("\nFinal Steane-7 dict: $(mc_result) \n")
        # if (mc_result[true_success_stat]  + mc_result[false_success_stat]) != genetic_params.num_shots
        #     throw(ErrorException("Some runs were invalid"))
        # end

        #fidelity = (round(mc_result[true_success_stat] / (mc_result[true_success_stat]+mc_result[false_success_stat]),digits=10))
        fidelities[idx] = 1 - sum(hamming_distances)/length(hamming_distances) # 1 is perfect alignment
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

function crossover(best_individuals, genetic_params, num_data_qubits)
    new_generation = deepcopy(best_individuals)

    parents = deepcopy(best_individuals)
    # could shuffle parents here

    # ADD CONDITION FOR ODD NUMBER: Keep size of generatios constnat 
    i = 1
    while i < length(parents)
        p1 = parents[i]
        p2 = parents[i+1]

        p1_size = length(p1.gates)
        p2_size = length(p2.gates)

        cp_1 = rand(1:p1_size-1)  # crossover point (in vector)
        cp_2 = rand(1:p2_size-1)  

        child1 = CircuitIndividual(cp_1 + p2_size-cp_2)
        child2 = CircuitIndividual(cp_2 + p1_size-cp_1)

        child1.gates[1:cp_1] = deepcopy(p1.gates[1:cp_1])
        child1.gates[cp_1+1:end] = deepcopy(p2.gates[cp_2+1:end])

        child2.gates[1:cp_2] = deepcopy(p2.gates[1:cp_2])
        child2.gates[cp_2+1:end] = deepcopy(p1.gates[cp_1+1:end])

        
        push!(new_generation, mutation(child1, genetic_params, num_data_qubits), mutation(child2, genetic_params, num_data_qubits))

        i += 2
    end

    if length(best_individuals)%2 != 0
        p1 = parents[1]
        p2 = parents[length(parents)]

        p1_size = length(p1.gates)
        p2_size = length(p2.gates)

        cp_1 = rand(1:p1_size-1)  
        cp_2 = rand(1:p2_size-1)  

        child1 = CircuitIndividual(cp_1 + p2_size-cp_2)#Circuit(nrows, ncols)
        child1.gates[1:cp_1] = deepcopy(p1.gates[1:cp_1])
        child1.gates[cp_2+1:end] = deepcopy(p2.gates[cp_2+1:end])

        push!(new_generation, mutation(child1, genetic_params, num_data_qubits))
    end

    return new_generation
end

function _random_single_qubit_gate(index)
    # choose from 1‑qubit gates you already define
    gates = (HadamardGate, IdentityGate, PauliXGate, PauliYGate, PauliZGate)
    return gates[rand(1:length(gates))](index)
end

function mutation(individual, genetic_params, num_data_qubits)
    rate = genetic_params.mutation_rate

    if rand() < rate
        # Pick one matrix elemenet randomly
        circuit_length = length(individual.gates)
        index = rand(1:circuit_length)

        if individual.gates[index] isa CNOT_Gate # mutate existing CNOT gates
            control = individual.gates[index].control
            target = individual.gates[index].target
            if rand() > 0.5 # SWAP control and target
                individual.gates[index] = CNOT_Gate(target, control)
            else
                individual.gates[index] = _random_single_qubit_gate(rand(1:num_data_qubits))
            end
        else  
            if rand()>0.5 # mutate single qubit gates into single-qubit gates on random qubit
                individual.gates[index] = _random_single_qubit_gate(rand(1:num_data_qubits))
            else  # ... or multi-qubit gates
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
            end
        end
    end

    return individual
end

function fitness_function(fidelities, circuit_sizes, gen, genetic_params)
    return 1000*fidelities - circuit_sizes #(gen/genetic_params.num_generations)*  # fitness can decrease over time since weighting is time-dependent
end

function run_genetic_search()

    ############## Define environment for GA ##########################

    # TODO: Define noise for the telegates (MOST IMPORTANT)
    # TODO: Outsource this to a function
    networking_params, genetic_params = define_parameters()                             # retrieve parameters
    # TODO: Mapping stage -> use dictionary to map indices to one another
    # As extracted from Hypergraph Partitoning
    #permutation = [1,7,4,2,3,5,6]
    qec_code_required_qubits = code_n(genetic_params.qec_code)
    permutation = collect(1:qec_code_required_qubits)#[1,2,3,4,5,6,7]
    inv_perm = invperm(permutation)
    mapping = perm_to_transpositions(deepcopy(permutation)) # careful: without deepcopy, this does in-place substitution of permutation    
    # NOTE: When generating the infromation for hypergraph part., we need to consult the naive encoding function in the logical encoding script to obtain the logical oeprators.
    # For the inversion of the circuit, we have a custoim function in circsim.jl since this requries applicaiton of correct indices, accounting for communication qubits.
    data_qubit_capacities = networking_params.register_sizes
    @assert sum(data_qubit_capacities) === qec_code_required_qubits
    num_registers  = length(data_qubit_capacities)
    register_lookup_array, data_qubits, num_data_qubits = create_lookup_array_cliff(data_qubit_capacities)      # create lookup array

    num_comm_qubits_per_register = num_registers-1
    num_qubits = num_data_qubits + num_comm_qubits_per_register*(num_registers) # one verification qubit
    all_qubits = collect(1:num_qubits)
    comm_qubits = setdiff(all_qubits, data_qubits)
    println("number of qubits is $num_qubits, number of comm. qubits per register is $num_comm_qubits_per_register")
    println("Lookup Array: $register_lookup_array")
    println("Data qubits: $data_qubits")

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
    code = MixedDestabilizer(genetic_params.qec_code)#S"XIXIXIX IXXIIXX IIIXXXX ZIZZIZI ZZIIZZI ZZIZIIZ IZIZIZI"
    code_stabilizer = stabilizerview(code)
    logical_z = logicalzview(code)
    println("Logical operators are $(logical_z)")
    target_state = vcat(code_stabilizer, logical_z)
    println("\nTarget state:$target_state and qubit size: $(qec_code_required_qubits) as well as logical qubit size: $(code_k(genetic_params.qec_code))")
    code_distance = distance(genetic_params.qec_code, DistanceMIPAlgorithm(solver=HiGHS))
    println("Code distance is $code_distance.\n\n")
    
    # instead, can also do MixedDestabiliser(Steane7()) and then extract the stabiliser tableau

    # We want to compare tableaus, so we canonicalize
    target_canon = canonicalize_rref!(target_state)
        #println("target canon: $target_canon")
        #println(typeof(target_canon))
    target_tableau = tab(target_canon[1])
    target_bit_matrix = tableau_to_bitmatrix(target_tableau)
        #println()
        #println("Target tableau: $target_tableau")
        #println("Bit Matrix: $target_bit_matrix")
    
    ########## Initialise population and run Genetic Algorithm #################
    
    fitness_evolution = Float64[]
    #winner_winner_chicken_dinner = 0 # Place holder for best circuit
    
    gen = 0
    population = initialise_population(genetic_params.num_individuals, num_data_qubits, warm_start=genetic_params.warm_start, qec_code = genetic_params.qec_code)
    
    # TODO: Add verification that warm start indeed yields a fidelity of 1!
    println("############ Generation #$gen ##########: Generation size: $(length(population))")
    println("")

    # evaluate population
    fidelities, circuit_sizes = evaluate_population(population, networking_params, genetic_params, mapping, inv_perm, register_lookup_array, data_qubits, comm_qubits, num_comm_qubits_per_register, num_qubits, target_bit_matrix, data_qubit_capacities, num_registers)
    #print("fidelity is $fidelities")
    fitness_scores = fitness_function(fidelities, circuit_sizes, gen, genetic_params) # TODO: Need to find a fair weighting here
    
    println("\n Generation $gen of size $(length(population)): Best fitness is $(maximum(fitness_scores)), where fidelity = $(fidelities[argmax(fitness_scores)]) and circuit size = $(circuit_sizes[argmax(fitness_scores)]) \n")
    
    push!(fitness_evolution, maximum(fitness_scores))
    
    while gen<genetic_params.num_generations

        gen += 1
        #@btime evaluate_population($population, $networking_params, $genetic_params, $mapping, $inv_perm, $register_lookup_array, $data_qubits, $comm_qubits, $num_comm_qubits_per_register, $num_qubits, $target_bit_matrix, $data_qubit_capacities, $num_registers)
        # perform selection
        best_individuals = selection(population, fitness_scores, tournament_size = genetic_params.tournament_size, selection_ratio = genetic_params.selection_ratio, num_elite = genetic_params.num_elite)
        # perform crossover (incl. mutations)
        new_generation = crossover(best_individuals, genetic_params, num_data_qubits)
        # apply mutations
        #mutated = mutations(new_generation, genetic_params)
        population = new_generation
        # evaluate population
        fidelities, circuit_sizes = evaluate_population(population, networking_params, genetic_params, mapping, inv_perm, register_lookup_array, data_qubits, comm_qubits, num_comm_qubits_per_register, num_qubits, target_bit_matrix, data_qubit_capacities, num_registers)
        fitness_scores = fitness_function(fidelities, circuit_sizes, gen, genetic_params)  # TODO: Need to find a fair weighting here
        if gen % 25== 0
        println("\n Generation $gen (/$(genetic_params.num_generations)) of size $(length(population)): Best fitness is $(maximum(fitness_scores)), where fidelity = $(fidelities[argmax(fitness_scores)]) and circuit size = $(circuit_sizes[argmax(fitness_scores)]) \n")
        end
        push!(fitness_evolution, maximum(fitness_scores))
        
    end

    # Extract best-performing individual
    winner_winner_chicken_dinner = population[argmax(fitness_scores)]
    winner_winner_chicken_dinner_size = circuit_sizes[argmax(fitness_scores)]
    #print_gate_matrix(winner_winner_chicken_dinner)

    winner_winner_chicken_dinner_circuit = construct_executable_circuit(networking_params.depolarising_noise, networking_params.gate_noise, networking_params.telegate_noise, winner_winner_chicken_dinner.gates, mapping, inv_perm, register_lookup_array, data_qubits, num_comm_qubits_per_register, num_qubits, data_qubit_capacities)
    #println(winner_winner_chicken_dinner_circuit)

    verification_logical_state = verify_success(winner_winner_chicken_dinner_circuit, target_state, num_qubits, data_qubits, num_registers)
    println("\nVerification successful (target state fidelity; only expressive (binary) in noiseless setting): $verification_logical_state")
    verification_logical_state = verification_logical_state == 1.0 ? true : false

    # Plot the evolution of fitness values
    #println("\n\nEvolution of fitness values: $fitness_evolution")
    plot_fitness_evol(fitness_evolution, networking_params, genetic_params, verification_logical_state)

    # @with classicalbitslayout => :expanded begin
    #     try
    #     savecircuit(
    #         winner_winner_chicken_dinner_circuit,
    #         "src/plots/GA/circuit_noise_GA_winner_size_$(winner_winner_chicken_dinner_size)_$(networking_params.telegate_noise)_ind$(genetic_params.num_individuals)_gen$(genetic_params.num_generations)_depth$(genetic_params.depth)_success_$(verification_logical_state)_warmstart_$(genetic_params.warm_start).png";
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

function verify_success(circuit, target_state, num_qubits, data_qubits, num_registers)
    push!(circuit, VerifyOp(target_state, data_qubits))
    
    initial_state = Register(one(MixedDestabilizer,num_qubits),num_registers*(num_registers-1))
    #print(mctrajectories(initial_state, circuit, trajectories=10000))
    mc_result = mctrajectories(initial_state, circuit, trajectories=10000)
    if (mc_result[true_success_stat]  + mc_result[false_success_stat]) != 10000
            throw(ErrorException("Some runs were invalid"))
    end
    fidelity = (round(mc_result[true_success_stat] / (mc_result[true_success_stat]+mc_result[false_success_stat]),digits=10))
    return fidelity
end


function plot_fitness_evol(fitness_evolution, networking_params, genetic_params, success)
    title_str = "Fitness Evolution : $(genetic_params.num_individuals) individuals over $(genetic_params.num_generations) generations"     
    subtitle_str = "fitness_evolution_telenoise_$(networking_params.telegate_noise)_success_$(success)_warmstart_$(genetic_params.warm_start)"
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="Generation", ylabel="Fitness", title=title_str, subtitle = subtitle_str)
    lines!(ax, 1:length(fitness_evolution), fitness_evolution)
    save("src/plots/GA/fitness_evolution_telenoise_$(networking_params.telegate_noise)_ind$(genetic_params.num_individuals)_gen$(genetic_params.num_generations)_success_$(success)_warmstart_$(genetic_params.warm_start).png", fig)
end

function print_gate_matrix(circ::Circuit)
    show(stdout, "text/plain", circ.gates)
    println()
end

end
