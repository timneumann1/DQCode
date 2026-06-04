# genetic.jl

"""
Genetic search for efficient encoding circuits of CSS QEC codes.
"""
module Genetic

export genetic_search

using ..Types
using ..Helper

using Random
using QuantumClifford
using QuantumClifford: MixedDestabilizer, @S_str, true_success_stat, false_success_stat, continue_stat, failure_stat, PauliMeasurement, VerifyOp, AbstractOperation
using Serialization, CSV, DataFrames
using ProgressMeter


"""    
    genetic_search(code_params::CodeParameters, network_specs::NetworkSpecifications, genetic_params::GeneticParameters, warm_start_gates::Vector{AbstractOperation})

Execute a genetic algorithm to search for an optimal encoding circuit for a given quantum error correction code.

### Input
- `code_params` -- properties of the target QEC CSS code
- `network_specs` -- underlying hardware network specifications and connectivity
- `genetic_params` -- optimisation hyperparameters for the genetic algorithm 
- `warm_start_gates` -- initial list of gates to seed the starting population

### Output
Returns a tuple containing:
- `GA_circ` -- the best circuit found by the algorithm
- `verification_logical_state` -- string indicating success or failure of logical state verification
- `best_gcounts` -- tuple of gate counts (single-, two-qubit and telegates) of the optimal circuit
- `fitness_evolution` -- tracking of maximum fitness score per generation
- `fidelity_evolution` -- tracking of fidelity score associated with highest-fitness individual per generation
- `gate_count_evolution` -- tracking of minimal gate counts associated with the best fitness per generation

### Notes
The encoding circuit of a logical zero state of a CSS code can always be described as a H-CNOT template.
Therefore, we extract the Hadamard layer from the seed/warm-start population of DQC-optimised encoding circuits,
and then perform an evolutionary search over the CNOT space.
Elite retention guarantees that at least one fidelity-1.0 circuit from the initial population is maintained 
throughout the entire evolution. For each generation, we track the best-performing individual with 
respect to the overall fitness (which is a function of fidelity and gate count); however, the overall best-performing
individual is only overwritten if has fidelity 1.0, guaranteeing we have a valid encoding circuit in the end.
"""
function genetic_search(code_params::CodeParameters, network_specs::NetworkSpecifications, 
                            genetic_params::GeneticParameters, warm_start_gates::Vector{AbstractOperation})    
    gen = 0
    population = initialise_population(genetic_params.num_individuals, warm_start_gates)   
    @info "Genetic search initialised with $(length(population)) individuals."
    best_circ_ind = Vector{AbstractOperation}() 
    best_gcounts = (typemax(Int), typemax(Int), typemax(Int))
    gate_count_evolution = Vector{Vector{Int64}}()
    fidelity_evolution = Float64[]
    fitness_evolution = Float64[]
    p = Progress(genetic_params.num_generations + 1; desc = "Genetic search", showspeed = true)
    #--------------- Evolution ----------------
    for gen in 0:genetic_params.num_generations
        fidelities, circuit_sizes = evaluate_population(population, code_params, network_specs, genetic_params)
        fitness_scores = fitness_function(fidelities, circuit_sizes, gen, genetic_params)  
        best_index = argmax(fitness_scores)
        best_fidelity = fidelities[best_index]
        if best_fidelity == 1.0
            best_circ_ind = population[best_index]
            best_gcounts = circuit_sizes[best_index]
        end
        next!(p; showvalues = [
        (:generation, gen),
        (:best_fitness, maximum(fitness_scores)),
        (:best_fidelity, fidelities[best_index]),
        (:best_gate_counts, circuit_sizes[best_index])
        ])
        push!(gate_count_evolution, circuit_sizes[best_index])
        push!(fidelity_evolution, best_fidelity)
        push!(fitness_evolution, fitness_scores[best_index])
        selected_individuals = selection(population, fitness_scores; tournament_size = genetic_params.tournament_size,
                                            num_elite = genetic_params.num_elite)
        population = crossover(genetic_params.num_individuals, selected_individuals, genetic_params.mutation_rate, 
                                network_specs.num_data_qubits, genetic_params.max_len, code_params.num_X_checks) 
    end
    finish!(p)
    GA_circ = best_circ_ind    
    verification_logical_state = verify_success(GA_circ, code_params.target_state, network_specs)
    @assert verification_logical_state "Verification of Genetic Algorithm circuit unsuccessful"
    return GA_circ, verification_logical_state, best_gcounts, fitness_evolution, fidelity_evolution, gate_count_evolution
end


"""
    initialise_population(num_individuals, warm_start_gates)::Vector{Vector{AbstractOperation}}

Create the starting population for the genetic algorithm with the warm-start circuit.

### Input
- `num_individuals` -- the number of identical individuals to generate for the initial population
- `warm_start_gates` -- sequence of quantum gates serving as the seed circuit for the search

### Output
A vector of circuit individuals representing the full initial population, where each 
individual is derived from the evaluated `warm_start_gates`.

### Notes
To maintain compatibility with the genetic algorithm's task of performing a CNOT circuit optimisation.
we enforce an H-CNOT structural template. The ability to separate the single-qubit from the two-qubit 
gates in the warm-start circuit relies on the specific construction of +1-phased CSS codes, as implemented
in `encoding_gott.jl``.
"""
function initialise_population(num_individuals::Int, warm_start_gates::Vector{AbstractOperation})::Vector{Vector{AbstractOperation}}
    @info "Initialising Population..."
    population = Vector{Vector{AbstractOperation}}(undef, num_individuals) 
    for i in eachindex(population)
        gates = copy(warm_start_gates)
        hadamard_indices = Vector{Int}()
        for g in gates
            if g isa AbstractSingleQubitOperator
                if !(g isa sHadamard)
                    @error "evolutionary search is currently only implemented for all-positive phases in the canonical tableau"
                else
                    push!(hadamard_indices, g.q)
                end
            end
        end
        filter!(g -> (g isa AbstractTwoQubitOperator), gates) # only keep the TwoQubitGates,and prepend the Hadamrd gates after
        for index in hadamard_indices
            pushfirst!(gates, sHadamard(index))
        end
        population[i] = gates
    end
    return population
end


"""
    evaluate_population(population::Vector{Vector{AbstractOperation}}, code_params::CodeParameters, 
                            network_specs::NetworkSpecifications, genetic_params::GeneticParameters)::Tuple{Vector{Float64}, Vector{Vector{Int64}}}

Evaluate the fitness of each circuit individual in the given population.

### Input
- `population` -- vector of circuits forming the current generation
- `code_params` -- parameters of the target QEC CSS code
- `network_specs` -- underlying hardware network connectivity and specifications
- `genetic_params` -- hyperparameter settings for the genetic search (including tableau distance metric)

### Output
The fidelity (defined as `1 - tableau_distance`) and a vector capturing the gate counts for each circuit.

### Notes
The fidelity of a circuit is calculated as `1 - tab_distance`, such that a tableau distance of zero indicates
a fidelity of 1.0. The tableau distance evaluates the bitmatrix agreement between desired logical zero state and 
the state obtained from applying the circuit onto the all-zero state.
"""
function evaluate_population(population::Vector{Vector{AbstractOperation}}, code_params::CodeParameters, 
                                network_specs::NetworkSpecifications, genetic_params::GeneticParameters)::Tuple{Vector{Float64}, Vector{Vector{Int64}}}
    fidelities = Vector{Float64}(undef, length(population))
    circuit_sizes = Vector{Vector{Int64}}(undef, length(population))
    for (idx, circuit) in enumerate(population)
        gcounts = gate_counts(circuit, network_specs)
        new_quantum_state = execute_circuit(circuit, network_specs.num_data_qubits)
        new_quantum_state_tab = tab(canonicalize_rref!( stabilizerview(new_quantum_state) )[1])
        new_quantum_state_bit_matrix = tableau_to_bitmatrix(new_quantum_state_tab)
        tab_distance = tableau_distance(new_quantum_state_bit_matrix, code_params.target_bit_matrix; metric = genetic_params.tableau_metric)
        fidelities[idx] = 1 - tab_distance 
        circuit_sizes[idx] =  gcounts
    end
    return fidelities, circuit_sizes
end


"""
    execute_circuit(circuit::Vector{AbstractOperation}, num_qubits::Int)::MixedDestabilizer

Execute the given sequence of gates on the initial all-zero quantum state.

### Input
- `circuit` -- the specified sequence of quantum gates to apply
- `num_qubits` -- total number of qubits of the quantum hardware

### Output
The final state of the quantum system, captured as `MixedDestabilizer` object.

### Notes
For circuit execution, we leverage the Monte Carlo trajectory simulation function `mctrajectory!`
provided by `QuantumClifford`.
"""
function execute_circuit(circuit::Vector{AbstractOperation}, num_qubits::Int)::MixedDestabilizer
    initial_state = Register(one(MixedDestabilizer,num_qubits),0)
    state, stat = mctrajectory!(initial_state, circuit)
    return state.stab
end


"""
    fitness_function(fidelities::Vector{Float64}, circuit_sizes::Vector{Vector{Int64}}, gen::Int, g::GeneticParameters)::Vector{Float64}

Calculate the fitness scores for the entire population based on fidelity and gate count penalties.

### Input
- `fidelities` -- vector of circuit fidelities (between 0.0 and 1.0)
- `circuit_sizes` -- vector containing the single-gate, two-gate, and telegate counts for each circuit
- `gen` -- the current generation index in the evolutionary process
- `g` -- genetic hyperparameters, including the vector of fitness weights

### Output
A vector containing the computed fitness score for each individual in the generation.

### Notes
The fitness score is a weighted combination of a circuit's fidelity and size (gate penalties). 
To promote structural exploration in the early stages of evolution, the penalty term is scaled 
by a coefficient that decays exponentially as a function of the current `gen`, thus decreasing 
the effective weight of the penalty throughout the course of optimisation.
This allows low-fidelity, small circuits to survive initially, while enforcing stricter fidelity
constraints towards the end of the evolution.
"""
function fitness_function(fidelities::Vector{Float64}, circuit_sizes::Vector{Vector{Int64} }, gen::Int, g::GeneticParameters)::Vector{Float64}
    penalties = map(cs -> sum(g.fitness_weights[2:4] .* cs) , circuit_sizes) 
    return g.fitness_weights[1] .* fidelities .- exp10(-gen/g.num_generations)*penalties 
end


"""
    selection(generation::Vector{Vector{AbstractOperation}}, fitness_scores::Vector{Float64}; 
                    tournament_size::Int=5, num_elite::Int = 1)::Vector{Vector{AbstractOperation}}

Select the highest-performing individuals from the current generation to form the parents of the next generation.

### Input
- `generation` -- vector of circuits representing the current population pool
- `fitness_scores` -- vector of corresponding numerical fitness scores for the population
- `tournament_size` -- (optional, default: `5`) number of randomly chosen individuals to compete in each tournament round
- `num_elite` -- (optional, default: `1`) number of top-performing individuals guaranteed to survive to the next generation

### Output
A vector containing the selected `best_individuals` (circuits) that will act as parents for the subsequent crossover stage.
We select half of the current population to serve as the parent population for the next generation.
"""
function selection(generation::Vector{Vector{AbstractOperation}}, fitness_scores::Vector{Float64}; 
                    tournament_size::Int=5, num_elite::Int = 1)::Vector{Vector{AbstractOperation}}
    length_generation = length(generation) 
    @assert length_generation == length(fitness_scores)
    num_selected = Int(floor(1/2*length_generation))
    elite_idx = sortperm(fitness_scores, rev=true)[1:num_elite]    
    best_individuals = generation[elite_idx]
    remaining = setdiff(collect(eachindex(generation)), elite_idx)
    for _ in 1:(num_selected-num_elite)
        tsize = min(tournament_size, length(remaining))
        tournament = remaining[randperm(length(remaining))[1:tsize]]
        best_idx = tournament[argmax(fitness_scores[tournament])]
        push!(best_individuals, generation[best_idx])
        deleteat!(remaining, findfirst(==(best_idx), remaining))
    end
    @assert length(best_individuals) == num_selected
    return best_individuals
end


"""
    crossover(num_individuals::Int, selected_individuals::Vector{Vector{AbstractOperation}}, mutation_rate::Float64, 
                num_data_qubits::Int, max_len::Int, num_hadamards::Int)::Vector{Vector{AbstractOperation}}

Generate a new population pool by performing crossover operations on the selected parent circuits.

### Input
- `num_individuals` -- required total size of the new generation
- `selected_individuals` -- vector of circuits chosen as parents during the selection phase
- `mutation_rate` -- probability of an individual mutating post-crossover
- `num_data_qubits` -- number of data qubits on the hardware, defining bounds for random gates
- `max_len` -- absolute maximum gate length allowed for generated circuits
- `num_hadamards` -- exact length of the preserved initial Hadamard layer 

### Output
A vector of circuit individuals containing the next generation of size `num_individuals`, 
including both the selected parents and their newly generated and mutated offspring.

### Notes
Crossover strictly guards the initial Hadamard layer defined by `num_hadamards`. Single point 
crossover occurs strictly within the trailing CNOT sequence, selecting a random (and potentially different)
crossover point per parent, and recombining the thus created vector parts to form the offsprings. 
We define helper functions to enforce that the offspring contain no back-to-back duplicate gates and 
do not cross the `max_len` threshold.
"""
function crossover(num_individuals::Int, selected_individuals::Vector{Vector{AbstractOperation}}, mutation_rate::Float64, 
                    num_data_qubits::Int, max_len::Int, num_hadamards::Int)::Vector{Vector{AbstractOperation}}
    new_generation = Vector{Vector{AbstractOperation}}()
    append!(new_generation, copy(selected_individuals))
    parents = copy(selected_individuals)
    i = 1
    while i < length(parents)
        p1 = parents[i]
        p2 = parents[i+1]
        p1_size = length(p1)
        p2_size = length(p2)
        cp_1 = rand(num_hadamards+1:p1_size-1)  # crossover point 
        cp_2 = rand(num_hadamards+1:p2_size-1)  # each parent has >= num_hadamards+2 gates (see `ensure_min_size`)

        child1 = Vector{AbstractOperation}(undef, cp_1 + p2_size-cp_2) 
        child2 = Vector{AbstractOperation}(undef, cp_2 + p1_size-cp_1) 
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
        child1 = _ensure_min_size!(child1, num_hadamards+2, num_data_qubits, num_hadamards)
        child2 = _ensure_min_size!(child2, num_hadamards+2, num_data_qubits, num_hadamards)
        child1 = _cap_individual_size(child1, max_len)
        child2 = _cap_individual_size(child2, max_len)

        push!(new_generation, child1, child2)
        i += 2
    end
    if num_individuals%2 != 0 # `2*length(parents) = num_individuals-1`
        child = copy(parents[end])
        child = mutation(child, mutation_rate, num_data_qubits)
        child = _clean_circuit(child)
        child = _ensure_min_size!(child, num_hadamards+2, num_data_qubits, num_hadamards)
        child = _cap_individual_size(child, max_len)
        push!(new_generation, child)
    end
    @assert length(new_generation) == num_individuals "generation size is not constant"
    return new_generation
end


"""
    _ensure_min_size!(ind::Vector{AbstractOperation}, min_len::Int, num_data_qubits::Int, num_hadamards::Int)::Vector{AbstractOperation}

Pad a circuit with random CNOT gates until it reaches a specified minimum length.

### Input
- `ind` -- the circuit to pad, modified in place
- `min_len` -- minimum gate length required
- `num_data_qubits` -- number of data qubits on the hardware
- `num_hadamards` -- exact length of the initial Hadamard layer, used to assert structural integrity

### Output
The modified vector of size greater than `min_len`.
"""
function _ensure_min_size!(ind::Vector{AbstractOperation}, min_len::Int, num_data_qubits::Int, num_hadamards::Int)::Vector{AbstractOperation}
    @assert length(ind) >= num_hadamards "Hadamard layer inflicted, something went wrong in the crossover"
    while length(ind) < min_len
        push!(ind, _random_two_qubit_gate(num_data_qubits))
    end
    return ind
end


"""
    _cap_individual_size(ind::Vector{AbstractOperation}, max_len::Int)::Vector{AbstractOperation}

Truncate a circuit if it exceeds a specified maximum gate length.

### Input
- `ind` -- the target circuit sequence to process
- `max_len` -- maximum gate length allowed

### Output
A vector of abstract operations strictly bounded by `max_len`.
"""
function _cap_individual_size(ind::Vector{AbstractOperation}, max_len::Int)::Vector{AbstractOperation}
    if length(ind) > max_len
        ind = ind[1:max_len]
    end
    return ind
end


"""
    _clean_circuit(circuit::Vector{AbstractOperation})::Vector{AbstractOperation}

Remove redundant operations from a given circuit sequence. Iteratively scan and remove
adjacent pairs of repeated gates to minimize the gate count of the individual.

### Input
- `circuit` -- the sequence of quantum gates to clean

### Output
An equivalent/cleaned circuit with no consecutive gate redundancy.

### Notes
Since `H` and `CNOT` are not only unitary but also involutory, they are self-inverse; thus, a sequence
H-H-CNOT-CNOT reduces to the identity.
"""
function _clean_circuit(circuit::Vector{AbstractOperation})::Vector{AbstractOperation}
    # Removes gate duplicates
    clean_circuit = Vector{AbstractOperation}() # will contain the clean gate sequence
    for gate in circuit
        if !isempty(clean_circuit) && clean_circuit[end] == gate
            pop!(clean_circuit)
        else
            push!(clean_circuit, gate)
        end
    end
   
    return clean_circuit
end


"""
    mutation(individual::Vector{AbstractOperation}, mutation_rate::Float64, num_data_qubits::Int)::Vector{AbstractOperation}

Apply random mutation operators to an individual circuit, given a specific probabilistic threshold.

### Input
- `individual` -- the circuit subjected to potential mutation steps
- `mutation_rate` -- probability (0.0 to 1.0) defining how likely the individual is to mutate
- `num_data_qubits` -- defining the range `1:num_data_qubits` for affected qubits of randomly sampled gates

### Output
The potentially modified `individual` circuit sequence.

### Notes
If the probabilistic threshold `mutation_rate` is surpassed by the random number generator, a single gate
is uniformly selected from the circuit. If this selected gate is a two-qubit logical operation, one of four disjoint actions 
occurs with varying probabilities.
"""
function mutation(individual::Vector{AbstractOperation}, mutation_rate::Float64, num_data_qubits::Int)::Vector{AbstractOperation}    
    if rand() < mutation_rate
        circuit_length = length(individual)
        index = rand(1:circuit_length)
        if individual[index] isa AbstractTwoQubitOperator
            control, target = affectedqubits(individual[index])
            r = rand()
            if r > 0.8 # SWAP control and target
                individual[index] = sCNOT(target, control)
            elseif r > 0.6 # replace gate with a random 2-qubit gate
                individual[index] = _random_two_qubit_gate(num_data_qubits)
            elseif r > 0.4 # append a random two qubit gate
                push!(individual, _random_two_qubit_gate(num_data_qubits))
            else # delete the gate
                deleteat!(individual, index)
            end
        end
    end
    return individual
end

"""
    _random_two_qubit_gate(num_data_qubits::Int)::AbstractTwoQubitOperator

Generate a random two-qubit gate with control and target selected from the available hardware qubits.

### Input
- `num_data_qubits` -- the total number of physical data qubits in the system block

### Output
A `sCNOT` operation with uniquely selected control and target indices.

### Notes
This function is currently restricted to CNOT operations to maintain integrity of the H-CNOT
encoding circuit template.
"""
function _random_two_qubit_gate(num_data_qubits::Int)::AbstractTwoQubitOperator
    control = rand(1:num_data_qubits)
    target = rand(1:num_data_qubits)
    while target == control
        target = rand(1:num_data_qubits)
    end
    return sCNOT(control, target)
end


end