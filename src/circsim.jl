# circsim.jl

module CircuitSimulator

using ..Types

using QuantumClifford
#using QuantumSavory: H, CNOT, X, Y, Z, stateof
using QuantumClifford: MixedDestabilizer, sHadamard, sCNOT, @S_str#, sX, SY, SZ
using QECCore

export create_lookup_array_cliff, execute_circuit, add_verification, add_telegate, add_noise, tensor_to_circuit, perm_to_transpositions

function create_lookup_array_cliff(num_data_qubits_per_register)
    
    # Lookup array is only targeting regular deta qubits (no mapping due to partitining), even though keeping track of the indices for the comm mapping later
    num_data_qubits = sum(num_data_qubits_per_register)
    #virtual_data_qubits = collect(1:num_data_qubits)
    #physical_data_qubits = [1,7,4,2,3,5,6]
   
    num_registers = length(num_data_qubits_per_register)
    num_comm_qubits_per_register = num_registers-1

    # Takes in an array containing the number of data qubits per module
    register_lookup_array = Vector{Int}(undef, num_data_qubits)#+num_comm_qubits_per_register*num_registers)
    #register_start_indices = Vector{Int}(undef, num_registers)
    data_qubits = Vector{Int}(undef, num_data_qubits)

    register = 1
    register_start_index = 1
    data_qubit_start_index = 1
    
    for j in num_data_qubits_per_register
        register_lookup_array[register_start_index:register_start_index+j-1] .= register #(num_comm_qubits_per_register)-1
        #register_lookup_array[physical_data_qubits[register_start_index:register_start_index+j-1]] .= register #(num_comm_qubits_per_register)-1

        #register_start_indices[register] = register_start_index
        data_qubits[register_start_index: (register_start_index+j-1)] = data_qubit_start_index: (data_qubit_start_index+j-1)  
        
        register_start_index+=j#+num_comm_qubits_per_register
        register +=1
        data_qubit_start_index += (j+num_comm_qubits_per_register)
    end
    #permute the register_lookup_array to match the physical data qubits
    println("Register lookup array:$register_lookup_array")
    #register_lookup_array = register_lookup_array[physical_data_qubits]
    #println(register_lookup_array)
    return register_lookup_array, data_qubits, num_data_qubits, num_comm_qubits_per_register
end


gate_to_apply(::Type{HadamardGate}, i::Int) = sHadamard(i)  # CliffordRepr  #TODO: verify that this is fixed in the next release
gate_to_apply(::Type{PauliXGate}, i::Int) = sX(i)
gate_to_apply(::Type{PauliYGate}, i::Int) = sY(i)
gate_to_apply(::Type{PauliZGate}, i::Int) = sZ(i)

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

function tensor_to_circuit(tensor, mapping, inv_perm, register_lookup_array, data_qubits, num_comm_qubits_per_register, target_state)
    
    #print("Comparing mapping and mapping 2: $mapping vs $mapping2")
    #idx(index::Int) =  index+num_comm_qubits_per_register * (register_lookup_array[index]-1)
    comm_idx(index::Int) = index + num_comm_qubits_per_register * (register_lookup_array[index]-1)  # applies correct mapping based on communication qubit structure
    #comm_perm_idx(index::Int) = permutation[index]+num_comm_qubits_per_register * (register_lookup_array[permutation[index]]-1) # applies correct mapping based on the inverse permutation of the qubit partitioning
    comm_inv_perm_idx(index::Int) = inv_perm[index]+num_comm_qubits_per_register * (register_lookup_array[inv_perm[index]]-1) # applies correct mapping based on the inverse permutation of the qubit partitioning

    circuit = [] # dimensions are (#rows in tensor+ num_comm_qubits x #cols in tensor)
    println("DATA")
    println(data_qubits)
    #println(idx.(collect(1:7)))
    #print("Test inverse: $(idx(2)): expect 4+1 = 5 ")
    #print("comm: $num_comm_qubits_per_register, $(register_lookup_array[7])")
    
    # Apply correct mapping
    println("Mapping")
    # Apply the inverse permutation of the mapping by applying transpoistions of inverse perm in left action <-> transpoistions of perm in right action via reverse(mapping) [the mapping contains transposition derived from the permutation, implementing it in left action]
    for (i,j) in reverse(mapping)
        #println(i,j)
        push!(circuit, sSWAP(comm_idx(i),comm_idx(j)))  # We could also use comm_perm_idx or comm_inv_perm_idx, since the relabeling based on the permutation conjugtes and thus fixes the permutation induces by the transposition SWAPS
        #println(idx(i), idx(j))
        println("Done")
    end
    #for (i,j) in mapping
    
    #print(tensor)

    # Add noise somewhere
    #print(typeof(tensor))

    for col in axes(tensor, 2) # each column corresponds to one layer

        # For operations, we need to use the comm_perm_idx function to correctly permute
        #@simlog sim "Entered the iteration $col in the outer loop"
        # single_qubit_gates = Dict{Type, Vector{Int}}() # stores single qubit gates and corresponding qubit indices per layer
        # CNOT_gates = Vector{Tuple{Int, Int}}()
        # CNOT_telegates = Vector{Tuple{Int, Int}}()
        #flags = Set{Int}() # to flag the control/target qubits that can be ignored

        for qubit in axes(tensor, 1) # each row corresponds to one qubit
            #register = register_lookup_array[qubit]
            #offset = num_comm_qubits_per_register*(register-1)
            gate = tensor[qubit,col]
            
            if gate isa Union{PauliXGate, PauliYGate, PauliZGate, HadamardGate} #  Gate && !(gate isa CNOT) 
                #gate_type = typeof(gate)
                #println(gate, qubit+offset)
                #println(idx(qubit))
                #println("HEEERE:$gate")
                #print(qubit)
                #println(idx(qubit))
                push!(circuit, gate_to_apply(typeof(gate),comm_inv_perm_idx(qubit)) ) 
            elseif gate isa CNOT_Gate  
                #qubit in flags && continue
                control = gate.control
                target = gate.target
                #target_register = register_lookup_array[target]
                #target_offset = num_comm_qubits_per_register*(target_register-1)
                if qubit == target
                    continue # only process the CNOTs via the control (all CNOTs and comm qubits are mutually exclusive)
                end

                if register_lookup_array[inv_perm[control]] == register_lookup_array[inv_perm[target]]
                    #println(gate, control+offset, target+offset)
                    push!(circuit, sCNOT(comm_inv_perm_idx(control), comm_inv_perm_idx(target) ))
                else 
                    # Perform telegate
                    #println("REMOTE gate between $control and $target: becomes $(idx(control)) and $(idx(target))")
                #println(gate, control+offset, target+target_offset)
                    push!(circuit, sCNOT(comm_inv_perm_idx(control),comm_inv_perm_idx(target)) )
                    push!(circuit, sCNOT(comm_inv_perm_idx(control),comm_inv_perm_idx(target)) )
                    push!(circuit, sCNOT(comm_inv_perm_idx(control),comm_inv_perm_idx(target)) )
                end
                #push!(flags, control)
                #push!(flags, target)
            elseif gate isa SWAP_Gate
                if qubit == gate.qubit_2
                    continue
                else
                    push!(circuit, sSWAP(comm_inv_perm_idx(gate.qubit_1), comm_inv_perm_idx(gate.qubit_2)))
                end
            end
        end
    end



    #push!(circuit, sSWAP(3,4))
    #push!(circuit, sSWAP(6,7))
    # Revert swapping for measurement of target state

    for (i,j) in mapping # reverse reverse  -> we do the reverse of the orginial permutation (the reverse transposition)
        #println("SWAAAPS")
        #println(swap_idx(i),swap_idx(j))
        #println(idx(i),idx(j))
        push!(circuit, sSWAP(comm_idx(i),comm_idx(j))) # We could also use comm_perm_idx, since the relabeling based on the permutation conjugtes and thus fixes the permutation induces by the transposition SWAPS
    end

    push!(circuit, VerifyOp(target_state, data_qubits)) 
    # Here, we should extract the tableau and compute the tableau distance, can we still sample many times

    return circuit


end


function add_telegate(circuit)
    # for comm_qubits:
    # num_of_connected_register - 1_{num_register < num_of_connected_register})

    # For telegate primitive, we only need to find the respective comm qubits of the registers 
    # (in lexicographical ordering: ~ num_of_connected_register - 1_{num_register < num_of_connected_register}),
    # and apply conditional operations on those qubits + the comm_inv_perm indices we normally apply ops on
end

function add_noise(circuit)
#
end

function add_verification(circuit, target_state, data_qubits)
    push!(circuit, VerifyOp(target_state , data_qubits)) 
    return circuit
end

function execute_circuit(circuit, num_qubits, num_traj)#, mode = "mc")
    initial_state = one(MixedDestabilizer,num_qubits)# S" IIIIIIZ IIIIIZI IIIIZII IIIZIII IIZIIII IZIIIII ZIIIIII"  # zero state
    return mctrajectories(initial_state, circuit, trajectories=num_traj)
end

function execute_circuit(circuit, num_qubits)#, mode = "pert")
    initial_state = one(MixedDestabilizer,num_qubits)# S" IIIIIIZ IIIIIZI IIIIZII IIIZIII IIZIIII IZIIIII ZIIIIII"  # zero state
    return petrajectories(initial_state, circuit)
end


end