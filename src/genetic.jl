module Genetic

using ..Types
using ..Helper

using Random
using QuantumClifford
using QuantumClifford: MixedDestabilizer, @S_str, true_success_stat, false_success_stat, continue_stat, failure_stat, PauliMeasurement, VerifyOp, AbstractOperation
using Serialization, CSV, DataFrames
using ProgressMeter

export genetic_search


function fitness_function(fidelities, circuit_sizes, gen, g)
    penalties = map(cs -> sum(g.fitness_weights[2:4] .* cs) , circuit_sizes) 
    # the exponential prefactor allows for some variability in the first third of the generations, if the parameters are tuned well (keep also some non-1 fidelity in the beginning)
    # in the end, we are only optimising while maintaining one fidelity
    return g.fitness_weights[1] .* fidelities .- exp10(-gen/g.num_generations)*penalties  # fitness can decrease over time since weighting is time-dependent
end

function initialise_population(num_individuals, warm_start_gates)#num_hadamards, num_data_qubits, warm_start)
    
    println("\nInitialising Population... \n")
    population = Vector{Vector{AbstractOperation}}(undef, num_individuals) # Vector{Circuit}(undef, num_individuals)
    
    for i in eachindex(population)
        gates = copy(warm_start_gates)
        # Make sure that the Hadamards are in the first layer to mitigate conflicts with the genetic search logic (which is a CNOT search)
        # this will not change the circuit since moving H to the beginning doesnt inflict commutation by the construction of baseline and MCTS (not in general of course)
        # we assume all SingleQubitGates in the warm start are Hadamard gates, since we are experimenting with CSS codes only
        hadamard_indices = [g.q for g in gates if g isa AbstractSingleQubitOperator]
        filter!(g -> (g isa AbstractTwoQubitOperator), gates) # only keep the TwoQubitGates,and prepend the Hadamrd gates after
        for index in hadamard_indices
            pushfirst!(gates, sHadamard(index))
        end
        population[i] = gates
    end
    return population
end


function evaluate_population(population, code_params, network_specs, genetic_params)

    fidelities = Vector{Float64}(undef, length(population))
    circuit_sizes = Vector{Vector{Int64} }(undef, length(population))
    
    for (idx, circuit) in enumerate(population)
       
        gcounts = gate_counts(circuit, network_specs)
        @assert gcounts[1] == code_params.num_X_checks

        new_quantum_state = execute_circuit(circuit, network_specs.num_data_qubits)
        new_quantum_state_tab = tab(canonicalize_rref!( stabilizerview(new_quantum_state) )[1])
 
        new_quantum_state_bit_matrix = tableau_to_bitmatrix(new_quantum_state_tab) # extract the stabiliser tableau from MixedDestabilizer object
        tab_distance = tableau_distance(new_quantum_state_bit_matrix, code_params.target_bit_matrix; metric = genetic_params.tableau_metric)
        
        fidelities[idx] = 1 - tab_distance # 1 is perfect alignment, here we are in the noiseless setting (one shot)
        circuit_sizes[idx] =  gcounts
    end

    return fidelities, circuit_sizes
end

function selection(generation, fitness_scores; tournament_size::Int=5, selection_ratio::Float64=1.0, num_elite = 1)
    length_generation = length(generation) 
    @assert length_generation == length(fitness_scores)
    num_selected = Int(floor(length_generation * selection_ratio))

    elite_idx = sortperm(fitness_scores, rev=true)[1:num_elite]    
    best_individuals = generation[elite_idx]
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
    new_generation = Vector{Vector{AbstractOperation}}()
    append!(new_generation, copy(selected_individuals)) #keep the selected individuals in the population
    #then make the selected individuals parents of the second half 
    parents = selected_individuals

    i = 1
    while i < length(parents)
        # first parts of respective partens are intentionally maintained in order to preserve H-CNOT structure
        p1 = parents[i]
        p2 = parents[i+1]

        p1_size = length(p1)
        p2_size = length(p2)

        # parent generation will always have at least num_hadamards+2 gates (see ensure_min_size), hence the below random choice is valid
        cp_1 = rand(num_hadamards+1:p1_size-1)  # crossover point (in vector)
        cp_2 = rand(num_hadamards+1:p2_size-1) 

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

    @assert length(new_generation) == num_individuals

    return new_generation
end

function _clean_circuit(circuit)
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


function _random_two_qubit_gate(num_data_qubits)
    control = rand(1:num_data_qubits)
    target = rand(1:num_data_qubits)
    while target == control
        target = rand(1:num_data_qubits)
    end
    return sCNOT(control, target)
end

function mutation(individual, mutation_rate, num_data_qubits)
    
    if rand() < mutation_rate
        # Pick one gate
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

        # We currently don't permute Hadamards in order to maintain some circuit structure

        end
    end

    return individual
end


function genetic_search(code_params, network_specs, genetic_params, warm_start_gates)#; warm_start = false, warm_start_gates = [], label = nothing)

    
    #------------ Initialise population and run Genetic Algorithm ------------
    
    gen = 0
    population = initialise_population(genetic_params.num_individuals, warm_start_gates)   
    @info "Genetic search initialised with $(length(population)) individuals"
    
    best_circ_ind = Vector{AbstractOperation}() #CircuitIndividual([])
    best_gcounts = (typemax(Int), typemax(Int), typemax(Int))

    gate_count_evolution = Vector{Vector{Int64}}()
    fidelity_evolution = Float64[]
    fitness_evolution = Float64[]

    #----------------------------- Evolution --------------------------------
    p = Progress(genetic_params.num_generations + 1; desc = "Genetic search", showspeed = true)
    
    # the evolution of fitness values is monotonously (potentially not strictly) increasing
    for gen in 0:genetic_params.num_generations
        # ----- Population Evaluation ----------
        fidelities, circuit_sizes = evaluate_population(population, code_params, network_specs, genetic_params)
        fitness_scores = fitness_function(fidelities, circuit_sizes, gen, genetic_params)  # TODO: Need to find a fair weighting here
        
        # By elite retention, the best individual in the offspring generation will always be the best individual overall, up to ties. 
        # however, this does not mean that it will have fidelity 1 (thus, we only overwrite best_circ_ind and gcounts if fidelity is 1)
        # Since we initialise with a fidelity 1 circuit and have elites, we are guaranteed to have a fidelity 1 circ in the end!
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

        # For plotting, we want to see the best individual in terms of fitness, irregardless of fidelity
        push!(gate_count_evolution, circuit_sizes[best_index])
        push!(fidelity_evolution, best_fidelity)
        push!(fitness_evolution, fitness_scores[best_index])
        
        # ----- Selection and Crossover ----------
        selected_individuals = selection(population, fitness_scores; tournament_size = genetic_params.tournament_size, selection_ratio = genetic_params.selection_ratio, num_elite = genetic_params.num_elite)
        population = crossover(genetic_params.num_individuals, selected_individuals, genetic_params.mutation_rate, network_specs.num_data_qubits, genetic_params.max_len, code_params.num_X_checks) 
    end
    GA_circ = best_circ_ind
    
    # ----- Verification ------
    verification_logical_state = verify_success(GA_circ, code_params.target_state, network_specs)
    @info "Verification of Genetic Algorithm circuit successful: $verification_logical_state"
    

    return GA_circ, verification_logical_state, best_gcounts, fitness_evolution, fidelity_evolution, gate_count_evolution
end


end
