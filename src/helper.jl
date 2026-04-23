module Helper

using ..Types

using QuantumClifford
using QuantumClifford: AbstractOperation, true_success_stat
using KaHyPar
using SparseArrays
using Quantikz: savecircuit, @with, classicalbitslayout
using CairoMakie

export create_lookup_array, comm_qubits_array, perm_to_transpositions, transpositions_to_perm
export data_qubit_partitioning, circuit_size, tableau_distance, tableau_to_bitmatrix
export execute_circuit, gate_to_apply, gates_to_circuit, verify_success, save_circuit_diagram
export next_run_dir, code_dirname
export compare_states
export plot_evolution

# lookup arrays needs to be created only once before executing the genetic search
#TODO: replace with cleaner version
# function create_lookup_array(register_sizes)
#     register_lookup_array = Vector{Int}(undef, sum(register_sizes))
#     register_start_indices = Vector{Int}(undef, length(register_sizes))
#     register = 1
#     register_start_index = 1
    
#     for i in eachindex(register_sizes)
#         j = register_sizes[i]
#         register_lookup_array[register_start_index:register_start_index+j-1] .= register
#         register_start_indices[register] = register_start_index
#         register_start_index+=j
#         register +=1
#     end
#     #print(register_lookup_array, register_start_indices)
#     return register_lookup_array, register_start_indices
# end


function mctrajectories_states(initialstate, circuit; trajectories::Int=500) # returns Vector{ QuantumClifford.MixedDestabilizer{ QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}} } }
    #counts = Dict{Tuple{typeof(initialstate), QuantumClifford.CircuitStatus}, Int}()
    stabilisers = Vector{ QuantumClifford.MixedDestabilizer{ QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}} } }()
    for i in 1:trajectories
        st, stat = QuantumClifford.mctrajectory!(copy(initialstate), circuit)
        #println("Type: $(typeof(st)), $(typeof(st.stab))")
        #println("$(fieldnames(typeof(st)))")
        if (stat==continue_stat)
            push!(stabilisers, st.stab)
            #print(typeof(st.stab))
        else
            throw("There were faulty circuit executions")
        end
        # println("Trajectory $i: status $stat")
        # println("Correct status?: $(stat==continue_stat)")
        # #counts[key] = get(counts, key, 0) + 1
        # #push!(stabilisers, st)
        # println("Current array of tableaus: $stabilisers")
    end
    return stabilisers
end


function execute_circuit(circuit, num_qubits::Int, num_registers::Int; num_traj::Int=500)#, keepstates::Bool=false)#, mode = "mc")
    initial_state = Register(one(MixedDestabilizer,num_qubits),num_registers*(num_registers-1))# S" IIIIIIZ IIIIIZI IIIIZII IIIZIII IIZIIII IZIIIII ZIIIIII"  # zero state  # we need num_communication_qubits slots in the classical register
    #print(fieldnames(typeof(Register(one(MixedDestabilizer, 1), 1))))
    #println(typeof(mctrajectories(initial_state, circuit, trajectories=num_traj)))
    #println("No. of trajectories:$(num_traj)")
    #print(circuit)
    return mctrajectories_states(initial_state, circuit, trajectories=num_traj)
end

function execute_circuit(circuit, initial_state::MixedDestabilizer{QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}}}; num_traj::Int=1)#, keepstates::Bool=false)#, mode = "mc")
    # Circuit execution for MCTS
    initial_state = Register(initial_state, 0)# S" IIIIIIZ IIIIIZI IIIIZII IIIZIII IIZIIII IZIIIII ZIIIIII"  # zero state  # we don't need num_communication_qubits slots in this simulation
    return mctrajectories_states(initial_state, circuit, trajectories=num_traj)
end


function execute_circuit(circuit, num_qubits::Int)#, mode = "pert")
    initial_state = one(MixedDestabilizer,num_qubits)# S" IIIIIIZ IIIIIZI IIIIZII IIIZIII IIZIIII IZIIIII ZIIIIII"  # zero state
    return petrajectories(initial_state, circuit)
end



function create_lookup_array(num_data_qubits_per_register)
    
    # Lookup array is only targeting data qubits (no mapping due to partitioning),
    # but also keeps track of the indices for the communication assignment mapping later
    # TODO: Reduce redundancy with genetic.jl file here
    num_data_qubits = sum(num_data_qubits_per_register)
    num_registers = length(num_data_qubits_per_register)
    num_comm_qubits_per_register = num_registers-1

    # Takes in an array containing the number of data qubits per module
    register_lookup_array = Vector{Int}(undef, num_data_qubits)
    data_qubits = Vector{Int}(undef, num_data_qubits)

    register = 1
    register_start_index = 1
    data_qubit_start_index = 1
    
    for j in num_data_qubits_per_register
        register_lookup_array[register_start_index:register_start_index+j-1] .= register 
        data_qubits[register_start_index: (register_start_index+j-1)] = data_qubit_start_index: (data_qubit_start_index+j-1)  

        register_start_index+=j
        register +=1
        data_qubit_start_index += (j+num_comm_qubits_per_register)
    end

    return register_lookup_array, data_qubits, num_data_qubits, num_comm_qubits_per_register
end

#function create_lookup_array(params)
# register_lookup_array = Int[]
# for (register, size) in enumerate(params.register_sizes)
#     append!(register_lookup_array, fill(register, size))
# end
#end

function comm_qubits_array(params)
    comm_qubits_array = Vector{Int}(undef,length(params.register_sizes))
    index = 1
    for i in eachindex(params.register_sizes)
        j = params.register_sizes[i]
        comm_qubits_array[i] = index
        index += j
    end
    return comm_qubits_array
end

function perm_to_transpositions(perm)
    n = length(perm)
    transpositions = Tuple{Int, Int}[]
    for i in n:-1:1
        if perm[i]!=i
            j = findfirst(==(i), perm)
            @assert !isnothing(j)
            push!(transpositions, (i, j))
            perm[j] = perm[i]
        end
    end
    return transpositions
end

function transpositions_to_perm(transpositions::Vector{Tuple{Int,Int}}, n::Int)
    perm = collect(1:n)
    for (i, j) in transpositions
        perm[i], perm[j] = perm[j], perm[i]
    end
    return perm
end

function circuit_size(circuit)
    return count(op ->
        !(op isa QuantumClifford.NoiseOp) &&
        !(op isa QuantumClifford.sSWAP) &&
        !(op isa QuantumClifford.sId1),
        circuit
    )
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
        # if (tableau.phases[r] ∉ (0,2))
        #     throw("Phase of the tableau is imaginary. Please investigate this case.")
        # end
        # bits[r, cols + 1] = 1/2*tableau.phases[r]  # phase +1 is represented as 0, phase -1 is represented as 2
        # -> positive phase is represented as 0, negative phase as 1
        #println(bits[r,:])
        bits[r, cols + 1] = tableau.phases[r] 
    end
    return bits
    
end

function tableau_distance(matrix::Matrix{Int}, target_matrix::Matrix{Int}; metric = "jaccard")
    #check that both marices have same dimensions
   
    # tableau distance without communication qubit overhead

    #cols_keep = vcat(data_qubits, size(matrix, 2))
    #matrix = matrix[:, cols_keep] 
    # we also want to eliminate the first #comm_qubits rows, which contain the stabilisers of the communication qubits after rref!().
    # the assumption is that the tableau is always in a product state of dataqubits and comm qubits, which is valid 
    # since the comm qubits get measured and then reset to zero
    #matrix = matrix[(length(comm_qubits)+1):end,:]
    
    # the target matrix already has the lexicographical ordering of data qubits, so no need to filter or change here
    @assert size(matrix) == size(target_matrix) "Mismatch in dimensions"
    # hamming count = 0
    # for rows
    #     for columns
    #         if numbers at repective positions different
    #             increae hamming count by one
    # return hammingcount/(rows*cols)
    #println("Final matrix: $matrix")
    #println("Final target: $target_matrix")
    difference_mask = matrix .!= target_matrix

    if metric == "hamming"
        return count(difference_mask) / length(matrix)
    elseif metric == "jaccard"
        support_mask = (matrix .!= 0) .| (target_matrix .!= 0)
        return count(difference_mask .& support_mask) / count(support_mask) # edge case of denom == 0 is trivial (n identity operators stabilise the state) and cannot occur for any valid stabilizer code
    end
end


function tableau_distance(matrix::Matrix{Int}, target_matrix::Matrix{Int}, data_qubits, comm_qubits, metric = "jaccard")
    #check that both marices have same dimensions
   
    # we want to compare the data qubits, so we only keep the columns (qubits) corresponding to data qubits AND the phase column
    cols_keep = vcat(data_qubits, size(matrix, 2))
    matrix = matrix[:, cols_keep] 
    # we also want to eliminate the first #comm_qubits rows, which contain the stabilisers of the communication qubits after rref!().
    # the assumption is that the tableau is always in a product state of dataqubits and comm qubits, which is valid 
    # since the comm qubits get measured and then reset to zero
    matrix = matrix[(length(comm_qubits)+1):end,:]
    
    # the target matrix already has the lexicographical ordering of data qubits, so no need to filter or change here
    @assert size(matrix) == size(target_matrix) "Check whether tracing out communication qubits (before passing to the tableau_distance() function) added identities or deleted rows"
    # hamming count = 0
    # for rows
    #     for columns
    #         if numbers at repective positions different
    #             increae hamming count by one
    # return hammingcount/(rows*cols)
    #println("Final matrix: $matrix")
    #println("Final target: $target_matrix")
    difference_mask = matrix .!= target_matrix

    if metric == "hamming"
        return count(difference_mask) / length(matrix)
    elseif metric == "jaccard"
        support_mask = (matrix .!= 0) .| (target_matrix .!= 0)
        return count(difference_mask .& support_mask) / count(support_mask) # edge case of denom == 0 is trivial (n identity operators stabilise the state) and cannot occur for any valid stabilizer code
    end
end

function _repair_partition_capacities(assignments::Vector{Int}, capacities::Vector{Int})
    k = length(capacities)
    
    # Count actual block sizes
    block_sizes = [count(==(b), assignments) for b in 1:k]
    println("Actual block sizes: $block_sizes")
    println("Required capacities: $capacities")
    
    # Find permutation of block labels that matches capacities
    # i.e. find which block should be relabelled as which register
    label_map = zeros(Int, k)
    remaining_blocks = collect(1:k)
    
    for (reg, cap) in enumerate(capacities)
        # Find a block whose size matches this register's capacity
        idx = findfirst(b -> block_sizes[b] == cap, remaining_blocks)
        if isnothing(idx)
            throw(ArgumentError(
                "No block with exactly $cap qubits found for register $reg. Please change the register sizes or relax the constraint on equal-weight partitions."))
            # fallback: assign whatever is left
            idx = 1
        end
        label_map[reg] = remaining_blocks[idx]
        deleteat!(remaining_blocks, idx)
    end
    
    println("Label map (register => old block): $label_map")

    # Apply relabeling
    repaired = similar(assignments)
    for (reg, old_block) in enumerate(label_map)
        repaired[assignments .== old_block] .= reg
    end
    
    return repaired
end


function data_qubit_partitioning(capacities, stabilizers)
    
    # traverse through rows in stabilisers and build up the connectivity graph

    k = length(capacities)
    nqubits = size(stabilizers, 2)
    @assert sum(capacities) == nqubits "Register capacities must sum to code length."
    
    println("Stabilizers: $stabilizers")
    # Build incidence matrix: rows=qubits (vertices), cols=stabilizers (hyperedges)
    I = Int[]
    J = Int[]
    e = 0
    for r in 1:size(stabilizers, 1)
        support = Int[]
        for q in 1:nqubits
            x, z = stabilizers[r, q]
            if x || z
                push!(support, q)
            end
        end
        # ignore trivial/singleton rows for communication objective
        if length(support) >= 2
            e += 1
            for q in support
                push!(I, q)
                push!(J, e)
            end
        end
        println("Support of stabilizer $r is $support\n")
    end

    if e == 0
        @warn "No multi-qubit stabilizer supports found; returning identity permutation."
        return collect(1:nqubits)
    end

    # apply KaHyPa to partition the graph
    println("Hyperedges:\nI:$I \nJ:$J")
    A = sparse(I, J, ones(Int, length(I)), nqubits, e)
    h = KaHyPar.HyperGraph(A)
    #println("Hypergraph:: $h")

    # exact balance for equal capacities; otherwise allow small imbalance then repair
    #equal_caps = all(c -> c == capacities[1], capacities)
    imbalance = 0.0 # equal_caps ? 0.0 : 0.2
    cfg = joinpath(@__DIR__, "km1_rKaHyPar_sea20.ini")
    partition = KaHyPar.partition(h, k; configuration=cfg)#  :edge_cut, imbalance=imbalance, seed = 42)
    #improved_partition = KaHyPar.improve_partition(h, k, partition; num_iterations=10, imbalance=imbalance, seed = 42)
    # normalize block ids to 1..k 
    min_id = minimum(partition)
    assignments = min_id == 0 ? Int.(partition .+ 1) : Int.(partition)
    #println("Partitioning is $assignments")
    assignments = _repair_partition_capacities(assignments, capacities)
    #println("Partitioning is $assignments")

    permutation = Int[]
    for b in 1:k
        append!(permutation, sort(findall(==(b), assignments)))
    end

    @assert length(permutation) == nqubits
    @assert sort(permutation) == collect(1:nqubits)

    println("KaHyPar assignments: $assignments")
    println("Data-qubit permutation: $permutation")
    return permutation
    # extract permutation from partitioning, e.g. [1,7,4,2,3,5,6] would mean that the 7th data qubit is in the first core if we have a [3,4] core assignment
end

gate_to_apply(::Type{HadamardGate}, i::Int) = sHadamard(i)  # CliffordRepr  #TODO: verify that this is fixed in the next release
gate_to_apply(::Type{PauliXGate}, i::Int) = sX(i)
gate_to_apply(::Type{PauliYGate}, i::Int) = sY(i)
gate_to_apply(::Type{PauliZGate}, i::Int) = sZ(i)
gate_to_apply(::Type{SGate}, i::Int) = sPhase(i)
gate_to_apply(::Type{InvSGate}, i::Int) = sInvPhase(i)
gate_to_apply(::Type{SqrtXGate}, i::Int) = sSQRTX(i)
gate_to_apply(::Type{InvSqrtXGate}, i::Int) = sInvSQRTX(i)
gate_to_apply(::Type{CX_Gate}, i::Int, j::Int) = sCNOT(i,j)
gate_to_apply(::Type{CZ_Gate}, i::Int, j::Int) = sCPHASE(i,j)

function gates_to_circuit(gates, n)
    # Converts gates to a circuit (same indexing), but counts gate overhead (in contrast to the below function which only constructs the circuit)
    gate_counts = [0,0,0]
    #@assert n.num_shots == 1
    circuit = Vector{QuantumClifford.AbstractOperation}()  
    for gate in gates
        T = typeof(gate)
        
        if T <: SingleQubitGate # isa Union{PauliXGate, PauliYGate, PauliZGate, HadamardGate, SGate} 
            #reward -= mdp.mcts_params.fitness_weights[2]
            gate_counts += [1,0,0]
            #new_quantum_state = execute_circuit([gate_to_apply(T, mdp.network_specs.inv_perm[qubit]) ], initial_quantum_state, num_traj = 1)
            push!(circuit, gate_to_apply(T, gate.index) ) 

        elseif T <: TwoQubitGate
            control = gate.control #network_specs.inv_perm[gate.control]
            target = gate.target #network_specs.inv_perm[gate.target]
            control_register = n.register_lookup_array[n.inv_perm[control]] 
            target_register = n.register_lookup_array[n.inv_perm[target]]
            
            if control_register == target_register # the lookup array does not account for the communication qubits
                #reward -= mdp.mcts_params.fitness_weights[3]
                gate_counts += [0,1,0]
            else
                #reward -= mdp.mcts_params.fitness_weights[4]
                gate_counts += [0,0,1]
            
            end
            push!(circuit, gate_to_apply(T, control, target ))
        
        else
            error("Unsupported gate type in gates_to_circuit: $(typeof(gate))")
            
            #new_quantum_state = execute_circuit([gate_to_apply(T, mdp.network_specs.inv_perm[control], mdp.network_specs.inv_perm[target]) ], initial_quantum_state, num_traj = 1)
        end
    end
    return circuit, gate_counts
end

function gates_to_circuit(gates)
   
    circuit = Vector{QuantumClifford.AbstractOperation}()  
    for gate in gates
        T = typeof(gate)
        if T <: SingleQubitGate # isa Union{PauliXGate, PauliYGate, PauliZGate, HadamardGate, SGate} 
            push!(circuit, gate_to_apply(T, gate.index) ) 
        elseif T <: TwoQubitGate
            control = gate.control #network_specs.inv_perm[gate.control]
            target = gate.target #network_specs.inv_perm[gate.target]
            push!(circuit, gate_to_apply(T, control, target ))
        else
            throw(ArgumentError("Unsupported gate type in gates_to_circuit: $(typeof(gate))"))
        end
    end
    return circuit
end

# function gates_to_circuit(gates, gate_set)
#     "Function to convert array of Main.DQCircuitSearch.Types.Gate gates to AbstractOperations object (e.g., for plotting)"

#     circuit = Vector{QuantumClifford.AbstractOperation}()

#     for gate in gates
#         T = typeof(gate)
#         if T in gate_set.single_qubit_gates
#             push!(circuit, gate_to_apply(T, gate.index))
#         elseif T in gate_set.two_qubit_gates
#             push!(circuit, gate_to_apply(T, gate.control, gate.target))
#         else
#             throw(ArgumentError("Unsupported gate type in gates_to_circuit: $(typeof(gate))"))
#         end
#     end

#     return circuit
# end


function save_circuit_diagram(gates::Vector{Gate}, directory, label)
    @with classicalbitslayout => :expanded begin
        try
        savecircuit(
            gates_to_circuit(gates),
            joinpath(directory, label);
            scale = 1
            
        )
        catch err
            @warn "savecircuit failed (circuit likely too large)" err
        end
    end
end

function save_circuit_diagram(circuit::Vector{AbstractOperation}, directory, label)
    @with classicalbitslayout => :expanded begin
        try
        savecircuit(
            circuit,
            joinpath(directory, label);
            scale = 1
            
        )
        catch err
            @warn "savecircuit failed (circuit likely too large)" err
        end
    end
end

function compare_states(test_state, target_state, n)
    circ = [VerifyOp(target_state, n.data_qubits)]
#    initial_state = Register(initial_state,n.num_registers*(n.num_registers-1))
    state, stat = mctrajectory!(test_state, circ)#, trajectories=1)
    if stat == true_success_stat
        return true
    elseif stat == false_success_stat
        return false
    else
        throw(ErrorException("Run was invalid due to status: $stat"))
    end
    # identity = ((round(mc_result[true_success_stat] / (mc_result[true_success_stat]+mc_result[false_success_stat]),digits=10)) == 1.0) ? true : false
    # return identity
end

function verify_success(circuit, initial_state, target_state, n)
    verification_circuit = copy(circuit)
    push!(verification_circuit, VerifyOp(target_state, n.data_qubits))
    
    initial_state = Register(initial_state,n.num_registers*(n.num_registers-1))
    state, stat = mctrajectory!(initial_state, verification_circuit)#, trajectories=n.num_shots)
    # if (mc_result[true_success_stat]  + mc_result[false_success_stat]) != 1# n.num_shots
    #         throw(ErrorException("Run was invalid"))
    # end
    # fidelity = mc_result[true_success_stat]# (round(mc_result[true_success_stat] / (mc_result[true_success_stat]+mc_result[false_success_stat]),digits=10))
    # return fidelity
    if stat == true_success_stat
        return true
    elseif stat == false_success_stat
        return false
    else
        throw(ErrorException("Run was invalid due to status: $stat"))
    end
end

function verify_success(circuit, target_state, n; comm_setting=false)
    verification_circuit = copy(circuit)
    if comm_setting
        push!(verification_circuit, VerifyOp(target_state, n.data_qubits))
    else
        push!(verification_circuit, VerifyOp(target_state, collect(1:n.num_data_qubits)))
    end
    initial_state = Register(one(MixedDestabilizer, n.num_qubits),n.num_registers*(n.num_registers-1))
    #print(mctrajectories(initial_state, circuit, trajectories=10000))
    state, stat = mctrajectory!(initial_state, verification_circuit)#, trajectories=n.num_shots)
    if stat == true_success_stat
        return true
    elseif stat == false_success_stat
        return false
    else
        throw(ErrorException("Run was invalid due to status: $stat"))
    end
    # print(mc_result)
    # print(mc_result[2])
    # if (mc_result[2][true_success_stat]  + mc_result[2][false_success_stat]) != 1# n.num_shots
    #         throw(ErrorException("Run was invalid"))
    # end
    # fidelity = mc_result[true_success_stat]#(round(mc_result[true_success_stat] / (mc_result[true_success_stat]+mc_result[false_success_stat]),digits=10))
    # return fidelity
end

function code_dirname(code)
    # s = lowercase(string(typeof(code)))                     # e.g. trivariatebicyclecode.trivariatebicycleviacirculantmat
    # s = replace(s, r"^main\.dqcircuitsearch\." => "")      # drop Main.DQCircuitSearch.
    # s = split(s, '.')[1]                                    # keep only "trivariatebicyclecode"
    # return s
    name = string(nameof(typeof(code)))   # "Steane7" or "TrivariateBicycleViaCirculantMat"
    name = replace(name, r"Via.*$" => "") # "TrivariateBicycle"
    name = replace(name, r"\d+$" => "")   # "Steane"
    return name
end

function next_run_dir(base_dir::AbstractString)
    mkpath(base_dir)
    runs = Int[]
    for name in readdir(base_dir)
        p = joinpath(base_dir, name)
        if isdir(p)
            v = tryparse(Int, name)
            if v !== nothing
                push!(runs, v)
            end
        end
    end
    next_id = isempty(runs) ? 1 : maximum(runs) + 1
    run_dir = joinpath(base_dir, string(next_id))
    mkpath(run_dir)
    return run_dir
end

function plot_evolution(dir, optimiser_label::String, fitness_scores, fidelities, gate_counts, genetic_params::GeneticParameters, success)
    title_str = "Evolution of optimiser metrics for $optimiser_label"     
    subtitle_str = "$(genetic_params.num_individuals) individuals over $(genetic_params.num_generations) generations -- Optimisation Success: $(success)"
    
    fig = Figure(size = (800, 900))

    ax_fit   = Axis(fig[1, 1], ylabel="Fitness", title=title_str, subtitle = subtitle_str)
    ax_gates = Axis(fig[2, 1], ylabel="Gate Counts")
    ax_fid   = Axis(fig[3, 1], xlabel="Generation", ylabel="Fidelity")

    generations = 1:length(fitness_scores)
    lines!(ax_fit, generations, fitness_scores, color=:chartreuse3, linewidth=2)

    single_q_counts = [g[1] for g in gate_counts]
    two_q_counts    = [g[2] for g in gate_counts]
    telegate_counts = [g[3] for g in gate_counts]

    lines!(ax_gates, generations, single_q_counts, label="Single-qubit", color=:goldenrod, linewidth=2)
    lines!(ax_gates, generations, two_q_counts,    label="Two-qubit",    color=:steelblue,  linewidth=2)
    lines!(ax_gates, generations, telegate_counts, label="Telegates",    color=:pucrimsonrple, linewidth=2)
    axislegend(ax_gates, position=:rt) 

    lines!(ax_fid, generations, fidelities, color=:green, linewidth=2)
    ylims!(ax_fid, -0.05, 1.05) 

    linkxaxes!(ax_fit, ax_gates, ax_fid)
    hidexdecorations!(ax_fit, grid=false)
    hidexdecorations!(ax_gates, grid=false)

    rowgap!(fig.layout, 10)

    outpath = joinpath(dir, "Optimisation_Evolution.png")
    save(outpath, fig)
end

function plot_evolution(dir, optimiser_label::String, fidelities, gate_counts, mcts_params::MCTSParameters, success)
    
    title_str = "Evolution of optimiser metrics for $optimiser_label"     
    subtitle_str = "$(mcts_params.n_iterations) Iterations over depth $(mcts_params.depth)-- Optimisation Success: $(success)"
    
    fig = Figure(size = (800, 900))

    ax_gates = Axis(fig[1, 1], ylabel="Gate Counts",  title=title_str, subtitle = subtitle_str)
    ax_fid   = Axis(fig[2, 1], xlabel="Step", ylabel="Fidelity")

    generations = 1:length(fidelities)
    lines!(ax_fid, generations, fidelities, color=:chartreuse3, linewidth=2)
    ylims!(ax_fid, -0.05, 1.05) 

    tick_labels = string.(generations)
    ax_gates.xticks = (generations, tick_labels)
    ax_fid.xticks = (generations, tick_labels)

    ax_gates.xticklabelrotation = pi/4
    ax_fid.xticklabelrotation = pi/4

    single_q_counts = [g[1] for g in gate_counts]
    two_q_counts    = [g[2] for g in gate_counts]
    telegate_counts = [g[3] for g in gate_counts]

    lines!(ax_gates, generations, single_q_counts, label="Single-qubit", color=:goldenrod, linewidth=2)
    lines!(ax_gates, generations, two_q_counts,    label="Two-qubit",    color=:steelblue,  linewidth=2)
    lines!(ax_gates, generations, telegate_counts, label="Telegates",    color=:crimson, linewidth=2)
    axislegend(ax_gates, position=:rt) 

    linkxaxes!(ax_gates, ax_fid)
    hidexdecorations!(ax_gates, grid=false)

    rowgap!(fig.layout, 10)

    outpath = joinpath(dir, "Optimisation_Evolution.png")
    save(outpath, fig)
end




end