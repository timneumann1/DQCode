module Genetic

using ..Types
using ..CircuitSimulator
using ..Helper
using ..LogicalEnc

using Random
using Quantikz: savecircuit, @with, classicalbitslayout
using QECCore: Steane7
using QuantumClifford: MixedDestabilizer, sHadamard, sCNOT, sSWAP, @S_str, true_success_stat, false_success_stat, continue_stat, failure_stat, PauliMeasurement, VerifyOp
using BenchmarkTools

export run_genetic_search


function define_parameters()

    #TODO: Rename or introduce a SimulationParameters type for this sim as well (in addition to DTS)
    networking_params = NetworkingParameters(
        [3,4], #register sizes
        0.05, # depolarising_prob 
        0.01, # gate_noise_prob 
    )

    genetic_params = GeneticParameters(
        100, # individuals
        100, # generations
        25000, # shots
        0.1,  # mutation rate
        5, # tournament size
        0.5, # selection_ratio
        )
    return networking_params, genetic_params
end


function steane_encoding_circuit(circuit)

    # could start with circuit from logical_encoding.jl here
    # (need to be recast from Vector{AbstractOperation} to tensor)
    
    circuit.gates[1,1] = HadamardGate()
    circuit.gates[2,1] = HadamardGate()
    circuit.gates[3,1] = HadamardGate()

    circuit.gates[7,2] = circuit.gates[4,2] = CNOT_Gate(7,4)

    circuit.gates[1,3] = circuit.gates[4,3] = CNOT_Gate(1,4)
    circuit.gates[7,3] = circuit.gates[5,3] = CNOT_Gate(7,5)

    circuit.gates[1,4] = circuit.gates[5,4] = CNOT_Gate(1,5)

    circuit.gates[1,5] = circuit.gates[6,5] = CNOT_Gate(1,6)

    circuit.gates[2,6] = circuit.gates[4,6] = CNOT_Gate(2,4)

    circuit.gates[2,7] = circuit.gates[6,7] = CNOT_Gate(2,6)

    circuit.gates[2,8] = circuit.gates[7,8] = CNOT_Gate(2,7)

    circuit.gates[3,9] = circuit.gates[5,9] = CNOT_Gate(3,5)

    circuit.gates[3,10] = circuit.gates[6,10] = CNOT_Gate(3,6)

    circuit.gates[3,11] = circuit.gates[7,11] = CNOT_Gate(3,7)

    circuit.gates[3,12] = circuit.gates[4,12] = SWAP_Gate(3,4)
    circuit.gates[6,12] = circuit.gates[7,12] = SWAP_Gate(6,7)

    #savecircuit(circuit, "src/plots/circuit_sim/circuit.png") # plotting is performed by enabling the reset function
    return circuit.gates
end


function initialise_population(num_individuals, num_data_qubits)
    population = Vector{Circuit}(undef, num_individuals)
    for i in eachindex(population)
        circ = Circuit(num_data_qubits, 12)         
        circ.gates = steane_encoding_circuit(circ)
        population[i] = circ
    end
    return population
end

function evaluate_population(population, networking_params, genetic_params, code, mapping, inv_perm, register_lookup_array, data_qubits, num_comm_qubits_per_register, num_qubits, target_state, data_qubit_capacities, num_registers)
    fitness_scores = Vector{Float64}(undef, length(population))
    for (idx, ind_tensor) in enumerate(population)
        quantum_clifford_circuit = tensor_to_circuit(code, networking_params.depolarising_noise, networking_params.gate_noise, ind_tensor.gates, mapping, inv_perm, register_lookup_array, data_qubits, num_comm_qubits_per_register, num_qubits, target_state, data_qubit_capacities)
        push!(quantum_clifford_circuit, VerifyOp(target_state, data_qubits)) 

        mc_result = execute_circuit(quantum_clifford_circuit, num_qubits, num_registers; num_traj=genetic_params.num_shots) # if specifying num_traj, we use MC sampling, otherwise perturbation.
        # for perturbative expansion, only the leading order is kept, so probabilities can be smaller than 1, 
        # also, PauliMeasurement don't work with pert. expansion currently

        #println("\nFinal Steane-7 dict: $(mc_result) \n")
        if (mc_result[true_success_stat]  + mc_result[false_success_stat]) != genetic_params.num_shots
            throw(ErrorException("Some runs were invalid"))
        end

        fidelity = (round(mc_result[true_success_stat] / (mc_result[true_success_stat]+mc_result[false_success_stat]),digits=3))
        fitness_scores[idx] = fidelity
    end
    return fitness_scores
end

function selection(generation, fitness_scores; tournament_size::Int=5, selection_ratio::Float64=1.0)
   length_generation = length(generation) 
    @assert length_generation == length(fitness_scores)
    num_selected = Int(floor(length_generation * selection_ratio))
    best_individuals = Vector{eltype(generation)}()
    remaining = collect(eachindex(generation))

    for _ in 1:num_selected
        tsize = min(tournament_size, length(remaining))
        tournament = remaining[randperm(length(remaining))[1:tsize]]
        # pick best fitness (max)
        best_idx = tournament[argmax(fitness_scores[tournament])]
        push!(best_individuals, generation[best_idx])
        deleteat!(remaining, findfirst(==(best_idx), remaining))
    end
    return best_individuals
end

function crossover(best_individuals, genetic_params)
    new_generation = deepcopy(best_individuals)

    parents = deepcopy(best_individuals)
    # could shuffle parents here

    i = 1
    while i < length(parents)
        p1 = parents[i]
        p2 = parents[i+1]

        nrows, ncols = size(p1.gates)
        cp = rand(1:ncols-1)  # crossover point (between columns)

        child1 = Circuit(nrows, ncols)
        child2 = Circuit(nrows, ncols)

        child1.gates = hcat(p1.gates[:, 1:cp], p2.gates[:, cp+1:end])
        child2.gates = hcat(p2.gates[:, 1:cp], p1.gates[:, cp+1:end])

        push!(new_generation, child1, child2)

        i += 2
    end

    return new_generation
end

function _random_single_qubit_gate()
    # choose from 1‑qubit gates you already define
    gates = (HadamardGate, IdentityGate, PauliXGate, PauliYGate, PauliZGate)
    return gates[rand(1:length(gates))]()
end

function mutations(new_generation, genetic_params)
    mutated = deepcopy(new_generation)
    rate = genetic_params.mutation_rate

    for ind in mutated
        if rand() < rate
            nrows, ncols = size(ind.gates)

            # pick a location that is NOT a 2‑qubit gate
            r, c = rand(1:nrows), rand(1:ncols)
            tries = 0
            while ind.gates[r, c] isa Union{CNOT_Gate, SWAP_Gate} && tries < 50
                r, c = rand(1:nrows), rand(1:ncols)
                tries += 1
            end

            # replace with a random 1‑qubit gate
            ind.gates[r, c] = _random_single_qubit_gate()
        end
    end

    return mutated
end

function run_genetic_search()

    ############## Define environment for GA ##########################

    # TODO: Define noise for the telegates (MOST IMPORTANT)
    # TODO: Outsource this to a function
    networking_params, genetic_params = define_parameters()                             # retrieve parameters
    # TODO: Mapping stage -> use dictionary to map indices to one another
    # As extracted from Hypergraph Partitoning
    permutation = [1,7,4,2,3,5,6]
    inv_perm = invperm(permutation)
    mapping = perm_to_transpositions(deepcopy(permutation)) # careful: without deepcopy, this does in-place substitution of permutation    
    # NOTE: When generating the infromation for hypergraph part., we need to consult the naive encoding function in the logical encoding script to obtain the logical oeprators.
    # For the inversion of the circuit, we have a custoim function in circsim.jl since this requries applicaiton of correct indices, accounting for communication qubits.
    data_qubit_capacities = networking_params.register_sizes
    num_registers  = length(data_qubit_capacities)
    register_lookup_array, data_qubits, num_data_qubits = create_lookup_array_cliff(data_qubit_capacities)      # create lookup array
    num_comm_qubits_per_register = num_registers-1
    num_qubits = num_data_qubits + num_comm_qubits_per_register*(num_registers) # one verification qubit
    println("number of qubits is $num_qubits, number of comm. qubits per register is $num_comm_qubits_per_register")
    println("Lookup Array: $register_lookup_array")
    println("Data qubits: $data_qubits")

    target_state = S"XIXIXIX IXXIIXX IIIXXXX ZIZZIZI ZZIIZZI ZZIZIIZ IZIZIZI"
    code = Steane7()
    
    ########## Initialise population and run Genetic Algorithm #################
    population = initialise_population(genetic_params.num_individuals, num_data_qubits)
    
    gen = 0
    while gen<genetic_params.num_generations
        # evaluate population
        fitness_scores = evaluate_population(population, networking_params, genetic_params, code, mapping, inv_perm, register_lookup_array, data_qubits, num_comm_qubits_per_register, num_qubits, target_state, data_qubit_capacities, num_registers)
        # perform selection

        best_individuals = selection(population, fitness_scores, tournament_size = genetic_params.tournament_size, selection_ratio = genetic_params.selection_ratio)
        # perform crossover
        new_generation = crossover(best_individuals, genetic_params)
        # apply mutations
        mutated = mutations(new_generation, genetic_params)
        println("\n Generation $gen: Best fidelity is $(maximum(fitness_scores))\n")
        population = mutated

        gen += 1
    end
    
    # Extract best-performing individual

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

end
