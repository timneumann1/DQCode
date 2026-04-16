# circsim.jl

module CircuitSimulator

using ..Types
using ..Helper

using QuantumClifford
using QuantumClifford: MixedDestabilizer, sHadamard, sCNOT, @S_str, Register, continue_stat#, sX, SY, SZ
using QECCore

import QuantumClifford: apply!, affectedqubits # we want to extend this with ConditionalGate
export execute_circuit, add_verification, add_telegate, add_noise, construct_executable_circuit



#We can also use NoisyGate: https://github.com/QuantumSavory/QuantumClifford.jl/blob/74ee758e87f5d7b1255d6747b346cff15ee10cea/docs/src/noisycircuits_ops.md

Pauli_gate_noise(::Type{PauliXGate}, p::Float64) = PauliNoise(p,0,0)
Pauli_gate_noise(::Type{PauliYGate}, p::Float64) = PauliNoise(0,p,0)
Pauli_gate_noise(::Type{PauliZGate}, p::Float64) = PauliNoise(0,0,p)
Pauli_gate_noise(::Type{HadamardGate}, p::Float64) = PauliNoise(p/2,0,p/2)


function apply!(state::Register, op::Types.ConditionalGate)
    if state.bits[op.controlbit]
        apply!(state, op.truegate)
    else
        apply!(state, op.falsegate)
    end
    return state
end


#function apply!(state::Register, op::projectZ)
#function projectZ!(s::AbstractStabilizer,qubit::Int;keep_result::Bool=true,phases::Bool=true)
#    project!(s, single_z(nqubits(s), qubit) ; keep_result, phases)
#end

function affectedqubits(op::Types.ConditionalGate)
    qs = Int[]
    append!(qs, collect(affectedqubits(op.truegate)))
    if op.falsegate !== nothing
       append!(qs, collect(affectedqubits(op.falsegate)))
    end
    return unique(qs)
end

# comm_idx(index::Int, n::NetworkingSpecifications) = index + n.num_comm_qubits_per_register * (n.register_lookup_array[index]-1)  # applies correct mapping based on communication qubit structure
#     #comm_perm_idx(index::Int) = permutation[index]+num_comm_qubits_per_register * (register_lookup_array[permutation[index]]-1) # applies correct mapping based on the inverse permutation of the qubit partitioning
# comm_inv_perm_idx(index::Int, n::NetworkingSpecifications) = n.inv_perm[index] + n.num_comm_qubits_per_register * (n.register_lookup_array[n.inv_perm[index]]-1) # applies correct mapping based on the inverse permutation of the qubit partitioning


function construct_executable_circuit(gates, gate_set, n; telegate_overhead = false)
    
    # n stands for networking specs

    # Depolarising channel: https://github.com/QuantumSavory/QuantumClifford.jl/blob/74ee758e87f5d7b1255d6747b346cff15ee10cea/src/noise.jl#L63-73

    #circuit = []
    circuit = Vector{QuantumClifford.AbstractOperation}()  
    num_single_qubit_gates = 0
    num_two_qubit_gates = 0
    num_telegates = 0
   
    # Apply the inverse permutation of the mapping by applying transpoistions of inverse perm in left action <-> transpoistions of perm in right action via reverse(mapping) [the mapping contains transposition derived from the permutation, implementing it in left action]
    for (i,j) in reverse(n.mapping)
        push!(circuit, sSWAP(n.comm_idx[i],n.comm_idx[j]))  # We could also use comm_perm_idx or comm_inv_perm_idx, since the relabeling based on the permutation conjugtes and thus fixes the permutation induces by the transposition SWAPS
    end
    
    # Add depolarising noise to all qubits at the beginning of the circuit
    #for data_qubit in collect(1:length(data_qubits))
    #    circuit = add_noise(circuit, depolarising_prob, comm_inv_perm_idx(data_qubit) )
    #end

    for gate in gates
        
        T = typeof(gate)
        if T in gate_set.single_qubit_gates# isa Union{PauliXGate, PauliYGate, PauliZGate, HadamardGate, SGate} 
            qubit = gate.index
            push!(circuit, gate_to_apply(T, n.comm_inv_perm_idx[qubit]) ) 
            num_single_qubit_gates += 1
            
        elseif T in gate_set.two_qubit_gates 
            control = gate.control
            target = gate.target
            control_register = n.register_lookup_array[n.inv_perm[control]] 
            target_register = n.register_lookup_array[n.inv_perm[target]] 
            

            if control_register == target_register # the lookup array does not account for the communication qubits
                push!(circuit, gate_to_apply(T, n.comm_inv_perm_idx[control], n.comm_inv_perm_idx[target] ))
                num_two_qubit_gates +=1
            else

                    # Perform telegate between control and target qubit in different registers, 
                    # i.e., push!(circuit, sCNOT(comm_inv_perm_idx(control),comm_inv_perm_idx(target)) ) remotely
                    # println("Performing a telegated between (unmapped) qubits $control and $target")
                if telegate_overhead
                    circuit = add_telegate(circuit, control, target, control_register, target_register, n)
                    num_telegates += 1
                else 
                    push!(circuit, gate_to_apply(T, n.comm_inv_perm_idx[control], n.comm_inv_perm_idx[target]))
                    num_telegates += 1
                end
            end
        end
    end

    # for col in axes(tensor, 2) # each column corresponds to one layer
    # # Add depolarising noise to all qubits
    #     # for data_qubit in collect(1:length(data_qubits))
    #     # circuit = add_noise(circuit, depolarising_prob, comm_inv_perm_idx(data_qubit) )
    #     # end
    #     for qubit in axes(tensor, 1) # each row corresponds to one qubit
            
    #         gate = tensor[qubit,col]

    #         if gate isa Union{PauliXGate, PauliYGate, PauliZGate, HadamardGate} #  Gate && !(gate isa CNOT) 
                
    #             # if qubit ==1 && gate != IdentityGate
    #             #     # Apply unfiform noise!
    #             #     # circuit = add_noise(circuit, gate_noise_prob, comm_inv_perm_idx(qubit), Main.DQCircuitSearch.Types.PauliZGate())
    #             #     # ^ Z error is harmless to the state, whereas an X error is destructive
    #             # end
             
    #             push!(circuit, gate_to_apply(typeof(gate),comm_inv_perm_idx(qubit)) ) 
            
    #         elseif gate isa CNOT_Gate  
    #             control = gate.control
    #             target = gate.target
    #             control_register = register_lookup_array[inv_perm[control]] 
    #             target_register = register_lookup_array[inv_perm[target]] 
             
    #             if qubit == target
    #                 continue # only process the CNOTs via the control (all CNOTs and comm qubits are mutually exclusive)
    #             end

    #             if control_register == target_register # the lookup array does not account for the communication qubits
    #                 push!(circuit, sCNOT(comm_inv_perm_idx(control), comm_inv_perm_idx(target) ))
    #             else

    #                 # Perform telegate between control and target qubit in different registers, 
    #                 # i.e., push!(circuit, sCNOT(comm_inv_perm_idx(control),comm_inv_perm_idx(target)) ) remotely
    #                 # println("Performing a telegated between (unmapped) qubits $control and $target")
    #                 circuit = add_telegate(circuit, control, target, control_register, target_register, num_comm_qubits_per_register, num_qubits, data_qubit_capacities, inv_perm, register_lookup_array, telegate_noise)
    #             end
            
    #         # elseif gate isa SWAP_Gate

    #         #     if qubit == gate.qubit_2
    #         #         continue
    #         #     else
    #         #         push!(circuit, sSWAP(comm_inv_perm_idx(gate.qubit_1), comm_inv_perm_idx(gate.qubit_2)))
    #         #     end
    #         end
    #     end
    # end

    # Revert swapping for measurement of target state
    for (i,j) in n.mapping # reverse reverse  -> we do the reverse of the orginial permutation (the reverse transposition)
        push!(circuit, sSWAP(n.comm_idx[i], n.comm_idx[j])) # We could also use comm_perm_idx, since the relabeling based on the permutation conjugtes and thus fixes the permutation induces by the transposition SWAPS
    end
    #circuit = add_noise(circuit, depolarising_prob, )
    # for data_qubit in collect(1:length(data_qubits))
    #     circuit = add_noise(circuit, depolarising_prob, comm_inv_perm_idx(data_qubit) )
    # end
    #add_noise(circuit, gate_noise_prob, comm_inv_perm_idx(2), Main.DQCircuitSearch.Types.PauliZGate() )
    #add_noise(circuit, gate_noise_prob, comm_inv_perm_idx(4), Main.DQCircuitSearch.Types.PauliZGate() )
    #add_noise(circuit, gate_noise_prob, comm_inv_perm_idx(6), Main.DQCircuitSearch.Types.PauliZGate() )

    return circuit, num_single_qubit_gates, num_two_qubit_gates, num_telegates
end


function add_telegate(circuit, control, target, control_register, target_register, n)

    if control_register == target_register
        throw("Ooops, this is not a telegate.")
    end
    
    #comm_inv_perm_idx(index::Int) = p.inv_perm[index] + p.num_comm_qubits_per_register * (p.register_lookup_array[p.inv_perm[index]]-1) # applies correct mapping based on the inverse permutation of the qubit partitioning
    control_comm_index = sum(n.register_sizes[1:control_register])+ (n.num_comm_qubits_per_register * (control_register-1)) + (control_register < target_register ? target_register-1 : target_register ) # start at the first communication qubit, and then see to which target register the control is connected
    target_comm_index = sum(n.register_sizes[1:target_register])+ (n.num_comm_qubits_per_register * (target_register-1)) + (target_register < control_register ? control_register-1 : control_register )
    #println("Communication qubits are at indices $control_comm_index and $target_comm_index")
    
    ### EJPP Protocol
    
    #circuit = add_noise(circuit, telegate_noise, comm_inv_perm_idx(control))
    #circuit = add_noise(circuit, telegate_noise, comm_inv_perm_idx(target))
    # Bell state entanglement
    # println("Telegate between qubits mapped_comm qubits $(comm_inv_perm_idx(control)) and $(comm_inv_perm_idx(target)).")
    
    
    push!(circuit, sHadamard(control_comm_index))
    push!(circuit, sCNOT(control_comm_index, target_comm_index))

    # println("Adding noise to comm qubit $control_comm_index and $target_comm_index")
    #circuit = add_noise(circuit, telegate_noise, control_comm_index)
    #circuit = add_noise(circuit, telegate_noise, target_comm_index)

    push!(circuit, sCNOT(n.comm_inv_perm_idx[control], control_comm_index))

    # Measurement 
    pauli_string_control = build_pauli_string_measurement(n.num_qubits, [control_comm_index])
    classical_register_index_control = control_comm_index - sum(n.register_sizes[1:control_register])
    meas_control = PauliMeasurement(pauli_string_control, classical_register_index_control)
    push!(circuit, meas_control)

    # Conditional Operations
    push!(circuit, Types.ConditionalGate(sX(target_comm_index),sId1(target_comm_index), meas_control.bit)) # perform the conditional operation
    push!(circuit, Types.ConditionalGate(sX(control_comm_index),sId1(control_comm_index), meas_control.bit)) # restore the |0> state in the control comm qubit
    
    push!(circuit, sCNOT(target_comm_index,n.comm_inv_perm_idx[target] ))
    push!(circuit, sHadamard(target_comm_index))
    
    # Measurement
    pauli_string_target = build_pauli_string_measurement(n.num_qubits, [target_comm_index])
    classical_register_index_target = target_comm_index - sum(n.register_sizes[1:target_register])
    meas_target = PauliMeasurement(pauli_string_target, classical_register_index_target)
    push!(circuit, meas_target)
    
    # Conditional Operations
    push!(circuit, Types.ConditionalGate(sZ(n.comm_inv_perm_idx[control]),sId1(n.comm_inv_perm_idx[control]), meas_target.bit))
    push!(circuit, Types.ConditionalGate(sX(target_comm_index),sId1(target_comm_index), meas_target.bit))  # restore the |0> state in the target comm qubit

    # Introduce noise to data qubits
    #circuit = add_noise(circuit, telegate_noise, comm_inv_perm_idx(control))
    #circuit = add_noise(circuit, telegate_noise, comm_inv_perm_idx(target))

    return circuit

end

function add_noise(circuit, prob::Float64) 
    """Circuit noise"""
    if prob<0 || prob > 1
        throw("Please provide a valid noise probabilty in [0,1]")
    end
    circuit_noise = NoiseOpAll(UnbiasedUncorrelatedNoise(prob));
    push!(circuit, circuit_noise)
    return circuit
end

function add_noise(circuit, prob::Float64, qubit) 
    """Circuit noise on single qubit"""
    if prob<0 || prob > 1
        throw("Please provide a valid noise probabilty in [0,1]")
    end
    circuit_noise = NoiseOp(UnbiasedUncorrelatedNoise(prob),[qubit]);
    push!(circuit, circuit_noise)
    return circuit
end

function add_noise(circuit, prob::Float64, qubit::Int, gate::Gate)
    # Special noise
    if prob<0 || prob > 1
        throw("Please provide a valid noise probabilty in [0,1]")
    end
    gate_noise_channel = Pauli_gate_noise(typeof(gate), prob)
    gate_noise = NoiseOp(gate_noise_channel, [qubit])
    push!(circuit, gate_noise)
    return circuit
end


function build_pauli_string_measurement(num_qubits::Int, qubits::Vector{Int})
    pauli = I # we can safely assume that the first qubit is a data qubit, since this is only false whenever there are zero qubits
    @inbounds for i in 2:(num_qubits) # traverses all data and comm qubits 
        pauli = (i in qubits) ? pauli⊗Z : pauli⊗I
    end
    return pauli
end



#=

function build_pauli_string_data_qubits(num_qubits::Int, data_qubits::Vector{Int})
    data_qubits_set = Set(data_qubits)
    pauli = Z # we can always assume that the first qubit is a data qubit, since this is only false whenever there are zero qubits
#Z⊗Z⊗Z⊗Z⊗I⊗Z⊗Z
    @inbounds for i in 2:(num_qubits) #traverses all data and comm qubits (not the verification one)
        pauli = (i in data_qubits_set) ? pauli⊗Z : pauli⊗I
    end
    print(pauli)
    return pauli
end

function measure_zero(circuit, data_qubits, num_qubits)
    #logical_zero_steane = P"XIXIXIX IXXIIXX IIIXXXX ZIZZIZI  ZZIIZZI ZZIZIIZ IZIZIZI"
    #print(data_qubits)
    zero_pauli_string = build_pauli_string_data_qubits(num_qubits, data_qubits)
    # pauli_string = ""
    # for i in collect(1:9) #TODO; replace with total number of qubits
    #     if i in data_qubits
    #         pauli_string += "Z"
    #     else
    #         pauli_string += "I"
    #     end
    # end
    print(zero_pauli_string)
    #println("Output of project")
    pauli_zero = PauliMeasurement(zero_pauli_string, 1)# Test: PauliMeasurement(P"-IIIZIIIIII", 1) yields true boolean when applied to a state |0>, since -Z|0> has eigenvalue -1, and  PauliMeasurement(P"IIIXIIIIII", 1) yield true or false in half of the cases
    #println(project!(data_qubits, zero_pauli_string))
    push!(circuit,pauli_zero) 

    #push!(circuit, project!(data_qubits, zero_pauli_string))

    #push!(circuit, PauliMeasurement(zero_pauli_string, 1)) 
    
    
    # #TODO: make the pauli string more general, e.g. by building it u
    # for data_qubit in data_qubits
    #     #make the pauli string more general!
    #     push!(circuit, PauliMeasurement(P"ZZZIZZZZI", 1))
    # end
    return circuit, pauli_zero
end

function naive_encoding_circuit_mapping(code, num_comm_qubits_per_register, register_lookup_array; undoperm=true)

    # TODO: Put this function somehwere else, but make sure it sees num_comm_qubits_per_register and register_lookup_array globally (same for the circuit builder function)
    comm_idx(index::Int) = index + num_comm_qubits_per_register * (register_lookup_array[index]-1)  # applies correct mapping based on communication qubit structure

    # Creatung the canonical tableau (without re-permuting the columns)
    n = code_n(code)
    k = code_k(code)
    #println("\nFor the given code, we have n=$n, k=$k. \n")
    
    code_standard_form, r, permx, permz = MixedDestabilizer(code, undoperm=false, reportperm=true); # undoperm without returns gives the orignal stabiliser tableau
    #println("Standard form of code is \n$code_standard_form")
    X = logicalxview(code_standard_form)
    Z = logicalzview(code_standard_form)
    
    # Creating te canonical tableau (with re-permuting the columns) for final state specification
    code_original_with_logicals = MixedDestabilizer(code, undoperm=true);
    # X_cleaned = logicalxview(md_cleaned)
    # println("X is $X_cleaned")
    # Z_cleaned = logicalzview(md_cleaned)
    # println("Z is $Z_cleaned")

    circ = QuantumClifford.AbstractOperation[]
    #push!(circ, Reset(initial_state, [1,2,3,4,5,6,7]))
    
    # Constructing the encoding circuit

    S = stabilizerview(code_standard_form)

    # # logical Xs
    # for i in 1:k
    #     for t in 1:n-k
    #         if X[i,t][1] == true
    #             push!(circ, sCNOT(comm_idx(n-k+i), comm_idx(t)))
    #         end
    #     end
    # end

    # projection on codespace
    for i in 1:r
        push!(circ, sHadamard(comm_idx(i)))
        if S[i,i][2] == true
            push!(circ, sPhase(comm_idx(i)))
        end
        for t in 1:n
            if i!=t
                xz = S[i,t]
                g = if xz == (true, true)  # Y
                    sZCY
                elseif xz == (true, false) # X
                    sZCX
                elseif xz == (false, true) && !(i<t<n-k+1) # Z
                    sZCZ
                end
                isnothing(g) || push!(circ, g(comm_idx(i),comm_idx(t)))
            end
        end
    end

    # correct for negative phases in the tableau
    for i in 1:n-k 
        if phases(S)[i]!=0
            if i<=r
                push!(circ, sZ(comm_idx(i)))
            else
                push!(circ, sX(comm_idx(i)))
            end
        end
    end

    # undoing the permutations to have the correct circuit for the final (re-permuted) tableau
    if undoperm
        perm = permx[permz]
        transpositions = perm_to_transpositions(perm)
        for (i,j) in transpositions
            push!(circ, sSWAP(comm_idx(i),comm_idx(j)))
        end
    end
    circ
end

=#



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


end