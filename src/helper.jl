module Helper

using ..Types

using QuantumClifford
using QuantumClifford: AbstractOperation, true_success_stat, mctrajectory!
using KaHyPar
using SparseArrays
using Quantikz: savecircuit, @with, classicalbitslayout
using CairoMakie

export create_lookup_array, comm_qubits_array, perm_to_transpositions, transpositions_to_perm
export data_qubit_partitioning, circuit_size, tableau_distance, tableau_to_bitmatrix
export execute_circuit, gate_to_apply, gate_counts, verify_success, save_circuit_diagram
export next_run_dir, code_dirname, save_txt
export compare_states
export plot_evolution
export qc_circuit_to_qasm


function mctrajectories_states(initial_state, circuit; trajectories::Int=500) # returns Vector{ QuantumClifford.MixedDestabilizer{ QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}} } }
    #counts = Dict{Tuple{typeof(initialstate), QuantumClifford.CircuitStatus}, Int}()
    stabilisers = Vector{ QuantumClifford.MixedDestabilizer{ QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}} } }()
    for i in 1:trajectories
        st, stat = mctrajectory!(copy(initial_state), circuit)
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

# # NOT USED SO FAR
# function execute_circuit(circuit, num_qubits::Int, num_registers::Int; num_traj::Int=500)#, keepstates::Bool=false)#, mode = "mc")
#     initial_state = Register(one(MixedDestabilizer,num_qubits),num_registers*(num_registers-1))# S" IIIIIIZ IIIIIZI IIIIZII IIIZIII IIZIIII IZIIIII ZIIIIII"  # zero state  # we need num_communication_qubits slots in the classical register
#     #print(fieldnames(typeof(Register(one(MixedDestabilizer, 1), 1))))
#     #println(typeof(mctrajectories(initial_state, circuit, trajectories=num_traj)))
#     #println("No. of trajectories:$(num_traj)")
#     #print(circuit)
#     return mctrajectories_states(initial_state, circuit, trajectories=num_traj)
# end

# USED IN MCTS SEARCH
function execute_circuit(circuit, initial_state::MixedDestabilizer{QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}}})#; num_traj::Int=1)#, keepstates::Bool=false)#, mode = "mc")
    # Circuit execution for MCTS
    initial_state = Register(initial_state, 0)# S" IIIIIIZ IIIIIZI IIIIZII IIIZIII IIZIIII IZIIIII ZIIIIII"  # zero state  # we don't need num_communication_qubits slots in this simulation
    state, stat = mctrajectory!(copy(initial_state), circuit)
    return state.stab
    #mctrajectories_states(initial_state, circuit, trajectories=num_traj)
end

# # NOT USED SO FAR
# function execute_circuit(circuit, num_qubits::Int, num_registers::Int; num_traj::Int=500)#, keepstates::Bool=false)#, mode = "mc")
#     initial_state = Register(one(MixedDestabilizer,num_qubits),num_registers*(num_registers-1))# S" IIIIIIZ IIIIIZI IIIIZII IIIZIII IIZIIII IZIIIII ZIIIIII"  # zero state  # we need num_communication_qubits slots in the classical register
#     #print(fieldnames(typeof(Register(one(MixedDestabilizer, 1), 1))))
#     #println(typeof(mctrajectories(initial_state, circuit, trajectories=num_traj)))
#     #println("No. of trajectories:$(num_traj)")
#     #print(circuit)
#     return mctrajectories_states(initial_state, circuit, trajectories=num_traj)
# end

#USED IN GA EVALUATION
function execute_circuit(circuit, num_qubits::Int)#, keepstates::Bool=false)#, mode = "mc")
    initial_state = Register(one(MixedDestabilizer,num_qubits),0)# S" IIIIIIZ IIIIIZI IIIIZII IIIZIII IIZIIII IZIIIII ZIIIIII"  # zero state  # we need num_communication_qubits slots in the classical register
    #print(fieldnames(typeof(Register(one(MixedDestabilizer, 1), 1))))
    #println(typeof(mctrajectories(initial_state, circuit, trajectories=num_traj)))
    #println("No. of trajectories:$(num_traj)")
    #print(circuit)
    state, stat = mctrajectory!(initial_state, circuit)
    return state.stab
    #return mctrajectories_states(initial_state, circuit, trajectories=num_traj)
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
    #data_qubit_dqc_indices = Vector{Int}(undef, num_data_qubits)

    register = 1
    register_start_index = 1
    #data_qubit_start_index = 1
    
    for j in num_data_qubits_per_register
        register_lookup_array[register_start_index:register_start_index+j-1] .= register 
        #data_qubit_dqc_indices[register_start_index: (register_start_index+j-1)] = data_qubit_start_index: (data_qubit_start_index+j-1)  

        register_start_index+=j
        register +=1
        #data_qubit_start_index += (j+num_comm_qubits_per_register)
    end

#    return register_lookup_array, data_qubit_dqc_indices, num_data_qubits, num_comm_qubits_per_register
    return register_lookup_array, num_data_qubits, num_comm_qubits_per_register

end



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
    for i in 1:n
        if perm[i]!=i
            j = findlast(==(i), perm)
            @assert !isnothing(j)
            push!(transpositions, (i, j))
            perm[j] = perm[i]
        end
    end
    return transpositions
end

function transpositions_to_perm(transpositions::Vector{Tuple{Int,Int}}, n::Int)
    perm = collect(1:n)
    for (i, j) in reverse(transpositions)
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

function qc_circuit_to_qasm(circ::Vector{AbstractOperation}) 

    # We convert a QuantumClifford circuit to QASM code

    """
    QASM code is of the form

        OPENQASM 2.0;
        include "qelib1.inc";
        qreg q[2];
        h q[0];
        cx q[0],q[1];
    
    Thus, we need to parse through the circuit identify the number of qubits needed,
        and then append the respective operations. For CSS codes, we only need to add Hadamard
        and CNOT gates.
    """
    
    max_q = 1
    for g in circ
        if g isa AbstractSingleQubitOperator
            if g.q > max_q 
                max_q = g.q
            end
        elseif g isa AbstractTwoQubitOperator
            if max(g.q1, g.q2) > max_q 
                max_q = max(g.q1, g.q2)
            end
        end
    end

    qasm_code = String[ "OPENQASM 2.0;", "include \"qelib1.inc\";" , "qreg q[$max_q];"    ]

    for g in circ
        if g isa AbstractSingleQubitOperator
            push!(qasm_code, "h q[$(g.q-1)];")
        elseif g isa AbstractTwoQubitOperator
            push!(qasm_code, "cx q[$(g.q1-1)],q[$(g.q2-1)];")
        end
    end

    return join(qasm_code, "\n")
end

function qasm_to_qc_circuit(qasm::String) 

    # We convert a QuantumClifford circuit to QASM code

    """
    QASM code is of the form

        OPENQASM 2.0;
        include "qelib1.inc";
        qreg q[2];
        h q[0];
        cx q[0],q[1];
    
    Thus, we need to parse through the qasm string line by line to identify which gate to add to the
    vector of QuantumClifford.AbstractOperation.
    
    For the CSS code encoding circuit, we only need to add Hadamard and CNOT gates.
    Additionally, there will be ancilla qubits (and CX to and from them), as well as measurements.
    """
    
    max_q = 1
    for g in circ
        if g isa AbstractSingleQubitOperator
            if g.q > max_q 
                max_q = g.q
            end
        elseif g isa AbstractTwoQubitOperator
            if max(g.q1, g.q2) > max_q 
                max_q = max(g.q1, g.q2)
            end
        end
    end

    qasm_code = String[ "OPENQASM 2.0;", "include \"qelib1.inc\";" , "qreg q[$max_q];"    ]

    for g in circ
        if g isa AbstractSingleQubitOperator
            push!(qasm_code, "h q[$(g.q-1)];")
        elseif g isa AbstractTwoQubitOperator
            push!(qasm_code, "cx q[$(g.q1-1)],q[$(g.q2-1)];")
        end
    end

    return return join(qasm_code, "\n")
end

function tableau_to_bitmatrix(tableau::QuantumClifford.Tableau{<:AbstractVector{UInt8}, <:AbstractMatrix{<:Unsigned}})
    rows, cols = size(tableau)
    bits = Matrix{Int}(undef, rows, cols+1)
    @inbounds for r in 1:rows
        for c in 1:cols
            x, z = tableau[r, c]      # every entry of the Tableau contains a tuple (x,z): (0,0) is I, (1,0) is X, (0,1) is Z, (1,1) is Y
            bits[r, c] = x + 2*z      # we map I, X, Z, Y to 0, 1, 2 and 3 to later determine the Hamming distance (in how many entries they disagree)
        end
        
        # if (tableau.phases[r] ∉ (0,2))
        #     throw("Phase of the tableau is imaginary. Please investigate this case.")
        # end
        
        bits[r, cols + 1] = tableau.phases[r] 
    end
    return bits
    
end



function tableau_distance(matrix::Matrix{Int}, target_matrix::Matrix{Int}; metric = "jaccard")
    # tableau distance without communication qubit overhead

    #check that both marices have same dimensions
    @assert size(matrix) == size(target_matrix) "Mismatch in dimensions"

    difference_mask = matrix .!= target_matrix

    if metric == "hamming"
        return count(difference_mask) / length(matrix)
    elseif metric == "jaccard"
        support_mask = (matrix .!= 0) .| (target_matrix .!= 0)
        return count(difference_mask .& support_mask) / count(support_mask) # edge case of denom == 0 is trivial (n identity operators stabilise the state) and cannot occur for any valid stabilizer code
    end
end


# function tableau_distance(matrix::Matrix{Int}, target_matrix::Matrix{Int}, data_qubits, comm_qubits, metric = "jaccard")
#     #check that both marices have same dimensions
   
#     # we want to compare the data qubits, so we only keep the columns (qubits) corresponding to data qubits AND the phase column
#     cols_keep = vcat(data_qubits, size(matrix, 2))
#     matrix = matrix[:, cols_keep] 
#     # we also want to eliminate the first #comm_qubits rows, which contain the stabilisers of the communication qubits after rref!().
#     # the assumption is that the tableau is always in a product state of dataqubits and comm qubits, which is valid 
#     # since the comm qubits get measured and then reset to zero
#     matrix = matrix[(length(comm_qubits)+1):end,:]
    
#     # the target matrix already has the lexicographical ordering of data qubits, so no need to filter or change here
#     @assert size(matrix) == size(target_matrix) "Check whether tracing out communication qubits (before passing to the tableau_distance() function) added identities or deleted rows"
    
#     difference_mask = matrix .!= target_matrix

#     if metric == "hamming"
#         return count(difference_mask) / length(matrix)
#     elseif metric == "jaccard"
#         support_mask = (matrix .!= 0) .| (target_matrix .!= 0)
#         return count(difference_mask .& support_mask) / count(support_mask) # edge case of denom == 0 is trivial (n identity operators stabilise the state) and cannot occur for any valid stabilizer code
#     end
# end


function data_qubit_partitioning(capacities, stabilizers)
    
    # traverse through rows in stabilisers and build up the connectivity graph

    k = length(capacities)
    nqubits = size(stabilizers, 2)
    @assert sum(capacities) == nqubits "Register capacities must sum to code length."
    
    #println("Stabilizers: $stabilizers")
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
        #println("Support of stabilizer $r is $support\n")
    end

    if e == 0
        @warn "No multi-qubit stabilizer supports found; returning identity permutation."
        return collect(1:nqubits)
    end

    # apply KaHyPa to partition the graph
    #println("Hyperedges:\nI:$I \nJ:$J")
    A = sparse(I, J, ones(Int, length(I)), nqubits, e)
    h = KaHyPar.HyperGraph(A)
    #println("Hypergraph:: $h")

    # exact balance for equal capacities; otherwise allow small imbalance then repair
    #equal_caps = all(c -> c == capacities[1], capacities)
    #imbalance = 0.0 # equal_caps ? 0.0 : 0.2

    #partition = KaHyPar.partition(h,k,configuration = :connectivity)
    cfg = joinpath(@__DIR__, "km1_kKaHyPar_sea20.ini")
    partition = KaHyPar.partition(h, k; configuration=cfg)#, seed = 42)#  :edge_cut, imbalance=imbalance, )
    
    #improved_partition = KaHyPar.improve_partition(h, k, partition; num_iterations=10, imbalance=imbalance, seed = 42)
    # normalize block ids to 1..k 

    #min_id = minimum(partition)
    #println("MINID PARTITION: $min_id")
    #assignments = min_id == 0 ? Int.(partition .+ 1) : Int.(partition)
    assignments = partition.+1 # KaHyPar is a Python-based optimiser based on 0 indexing
    #println("Partitioning is $assignments")
    #assignments = _repair_partition_capacities(assignments, capacities)
    #println("Partitioning is $assignments")
    block_sizes = [count(==(b), assignments) for b in 1:k]
    #println("Actual block sizes: $block_sizes")
    #println("Required capacities: $capacities")
    if block_sizes != capacities
        error("KaHyPar failed to satisfy the capacity constraint with correct block sizes, please change the capacities or alter the mapping manually")
    end
    mapping = Int[]
    for b in 1:k
        append!(mapping, sort(findall(==(b), assignments)))
    end

    @assert length(mapping) == nqubits
    @assert sort(mapping) == collect(1:nqubits)

    # mapping here indicates for each position, which qubit will sit there
    #println("KaHyPar assignments: $assignments")
    #println("Data-qubit mapping: $mapping")
    return mapping
    # extract permutation from partitioning, e.g. [1,7,4,2,3,5,6] would mean that the 7th data qubit is in the first core if we have a [3,4] core assignment
end



function gate_counts(circuit, n)
    # Converts gates to a circuit (same indexing), but counts gate overhead (in contrast to the below function which only constructs the circuit)
    
    mapping = copy(n.inv_map)

    gate_counts = [0,0,0]

    for op in circuit
        T = typeof(op)
        if T <: AbstractSingleQubitOperator
            #push!(encoding_circuit_orig, T(perm[op.q]))
            gate_counts += [1,0,0]
        elseif T<: AbstractTwoQubitOperator
            #println("Gate type is $T")
            control_register = n.register_lookup_array[mapping[op.q1]] 
            target_register = n.register_lookup_array[mapping[op.q2]]
            
            #push!(encoding_circuit_orig, T(p_control, p_target))
            if control_register == target_register
                # the above is equivalent to network_specs.register_lookup_array[findfirst(==(p_target), network_specs.permutation)]
                gate_counts += [0,1,0]
            else
                if T == sSWAP
                    gate_counts += [0,0,3] # a SWAP can be decomposed as three CNOT gates, each of which needs to be performed inter-core
                else 
                    gate_counts += [0,0,1]
                end
            end
            # For a SWAP gate, we assume that the SWAP operation is actually performed (qubit teleportation), after which we update the mapping
        else
            error("Unsupported gate type in gates_to_circuit: $T")
        end
    end
    
    return gate_counts
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


# CURRENTLY NOT USED
# function verify_success(circuit, initial_state, target_state, n)
#     verification_circuit = copy(circuit)
#     push!(verification_circuit, VerifyOp(target_state, n.data_qubits))
    
#     initial_state = Register(initial_state,n.num_registers*(n.num_registers-1))
#     state, stat = mctrajectory!(initial_state, verification_circuit)#, trajectories=n.num_shots)
#     # if (mc_result[true_success_stat]  + mc_result[false_success_stat]) != 1# n.num_shots
#     #         throw(ErrorException("Run was invalid"))
#     # end
#     # fidelity = mc_result[true_success_stat]# (round(mc_result[true_success_stat] / (mc_result[true_success_stat]+mc_result[false_success_stat]),digits=10))
#     # return fidelity
#     if stat == true_success_stat
#         return true
#     elseif stat == false_success_stat
#         return false
#     else
#         throw(ErrorException("Run was invalid due to status: $stat"))
#     end
# end

# AS USED IN GENETIC, MCTS and GOTTESMAN as well as mqt and qiskit baseline encoding
function verify_success(circuit, target_state, n)#; comm_setting=false)
    verification_circuit = copy(circuit)
    push!(verification_circuit, VerifyOp(target_state, n.data_qubits))
   
    initial_state = Register(one(MixedDestabilizer, n.num_data_and_comm_qubits),n.num_registers*(n.num_registers-1))
    state, stat = mctrajectory!(initial_state, verification_circuit)#, trajectories=n.num_shots)
    if stat == true_success_stat
        return true
    elseif stat == false_success_stat
        return false
    else
        throw(ErrorException("Run was invalid due to status: $stat"))
    end
end

function code_dirname(code)
    name = string(nameof(typeof(code)))   # "Steane7" or "TrivariateBicycleViaCirculantMat"
    name = replace(name, r"Via.*$" => "") # "TrivariateBicycle"
    name = replace(name, r"\d+$" => "")   # "Steane"
    return name
end

#NOT USED
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


function save_txt(folder, title, obj)
    open(joinpath(folder, title), "w") do io
        for fn in fieldnames(typeof(obj))
            println(io, fn, " = ", repr(getfield(obj, fn)))
        end
    end
end


end