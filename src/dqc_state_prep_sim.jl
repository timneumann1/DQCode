module DQCLogicalStatePrepSimulator

using ..Types
using ..Helper
using ..QECTools

using QuantumClifford
using QuantumClifford: MixedDestabilizer, sHadamard, sCNOT, @S_str, Register, continue_stat, AbstractOperation, AbstractStabilizer, AbstractNoise, apply_single_x!, apply_single_y!, apply_single_z! #, sX, SY, SZ
using QECCore

import QuantumClifford: apply!, affectedqubits, applynoise! # we want to extend this with ConditionalGate
export add_telegate, add_noise, construct_executable_circuit
export dqc_state_prep


#We can also use NoisyGate: https://github.com/QuantumSavory/QuantumClifford.jl/blob/74ee758e87f5d7b1255d6747b346cff15ee10cea/docs/src/noisycircuits_ops.md

Pauli_gate_noise(::Type{PauliXGate}, p::Float64) = PauliNoise(p,0,0)
Pauli_gate_noise(::Type{PauliYGate}, p::Float64) = PauliNoise(0,p,0)
Pauli_gate_noise(::Type{PauliZGate}, p::Float64) = PauliNoise(0,0,p)
#Pauli_gate_noise(::Type{HadamardGate}, p::Float64) = PauliNoise(p/2,0,p/2)


function apply!(state::Register, op::Types.ConditionalGate)
    #println(state)
    #println("State.bits : $(state.bits)")
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

function affectedqubits(op::SingleQubitGate)
    qs = Int[]
    append!(qs, op.index)
    return qs
end

function affectedqubits(op::TwoQubitGate)
    qs = Int[]
    append!(qs, op.control)
    append!(qs, op.target)
    return qs
end


# comm_idx(index::Int, n::NetworkingSpecifications) = index + n.num_comm_qubits_per_register * (n.register_lookup_array[index]-1)  # applies correct mapping based on communication qubit structure
#     #comm_perm_idx(index::Int) = permutation[index]+num_comm_qubits_per_register * (register_lookup_array[permutation[index]]-1) # applies correct mapping based on the inverse permutation of the qubit partitioning
# comm_inv_perm_idx(index::Int, n::NetworkingSpecifications) = n.inv_perm[index] + n.num_comm_qubits_per_register * (n.register_lookup_array[n.inv_perm[index]]-1) # applies correct mapping based on the inverse permutation of the qubit partitioning

function build_layers(gates, n)

    gate_list = copy(gates)
    layers = Vector{Vector{Gate}}()
    
    while !isempty(gate_list)

        idx = 1
        layer = Vector{Gate}()
        qubit_used_in_layer = falses(n.num_data_qubits)  
        del_gates = Vector{Int}()

        while idx <= length(gate_list)
            gate = gate_list[idx]
            if !( any(qubit_used_in_layer[affectedqubits(gate)]) ) 
                push!(layer, gate)
                push!(del_gates, idx)
            end
            qubit_used_in_layer[affectedqubits(gate)] .= true
            #push!(layers, layer)
            
            idx += 1
            # println("Layer: $layer")
            #     layer = Vector{Gate}()
            #     qubit_used_in_layer = falses(n.num_data_qubits) 
                
            # else 
            #     push!(layer, gate)
            #     qubit_used_in_layer[affectedqubits(gate)] .= true
            #     idx +=1
            #     if idx > length(gates)
            #         push!(layers, layer)
            #     end
            # end
        end
        push!(layers, layer)
        deleteat!(gate_list, del_gates)
        println("Layer: $layer")
        # for (remove_count, del_idx) in enumerate(del_gates)
        #     remove_count -= 1 # the remove count should start at 0
        #     deleteat!(gate_list, del_idx-remove_count)
        # end

    end
    return layers
end

function construct_DQC_executable_circuit(gates, n, noise)
    
    # n stands for networking specs

    # Depolarising channel: https://github.com/QuantumSavory/QuantumClifford.jl/blob/74ee758e87f5d7b1255d6747b346cff15ee10cea/src/noise.jl#L63-73

    #circuit = []
    circuit = Vector{QuantumClifford.AbstractOperation}()  
   
    # Apply the inverse permutation of the mapping by applying transpoistions of inverse perm in left action <-> transpoistions of perm in right action via reverse(mapping) [the mapping contains transposition derived from the permutation, implementing it in left action]
    for (i,j) in reverse(n.mapping)
        push!(circuit, sSWAP(n.comm_idx[i],n.comm_idx[j]))  # We could also use comm_perm_idx or comm_inv_perm_idx, since the relabeling based on the permutation conjugtes and thus fixes the permutation induces by the transposition SWAPS
    end
    
    # Add depolarising noise to all qubits at the beginning of the circuit
    #for data_qubit in collect(1:length(data_qubits))
    #    circuit = add_noise(circuit, depolarising_prob, comm_inv_perm_idx(data_qubit) )
    #end
    
    ### Build ASAP layers
    
    #go through circuit and accumulate layers (we may assume that no further grouping for telegates can be done)
    println("CIRCUIT: $gates")
    # adapted from MQT circuit_utils.py file
    
    layers = build_layers(gates, n)

    #ADD init noise to all qubits init_noise
    
    add_noise(circuit, [n.comm_inv_perm_idx[data_q] for data_q in collect(1:n.num_data_qubits)], noise.init_noise)

    for layer in layers
    #for every layer, add gate noise for every normal gate and idling noise for any idle qubits (if there is telegates in the layer, increase the idle probability!)
   # only after that, insert the telegate gadgets in place (this includes the folloing noise: comm init noise, two qubit depolarising for each gate , measeuremtn noise, classical noise, one more gate noise for single quits )
    
        #find idling gates and apply idle_depolarising_noise

        affected_qubits = Set( Iterators.flatten( [affectedqubits(gate) for gate in layer] ) ) 
        idle_qubits = setdiff(1:n.num_data_qubits, affected_qubits)
        idle_qubits_DQC = [n.comm_inv_perm_idx[idle_q] for idle_q in idle_qubits]
    
        telegates_layer = false
        for gate in layer
            T = typeof(gate)
            if T <: SingleQubitGate# isa Union{PauliXGate, PauliYGate, PauliZGate, HadamardGate, SGate} 
                qubit = gate.index
                DQC_qubit = n.comm_inv_perm_idx[qubit]
                push!(circuit, gate_to_apply(T, DQC_qubit) ) 
                #num_single_qubit_gates += 1
                add_noise(circuit, [DQC_qubit], noise.single_q_gate_noise)
                
            elseif T <: TwoQubitGate
                control = gate.control
                target = gate.target
                DQC_control = n.comm_inv_perm_idx[control]
                DQC_target = n.comm_inv_perm_idx[target]
                control_register = n.register_lookup_array[n.inv_perm[control]] 
                target_register = n.register_lookup_array[n.inv_perm[target]] 
                

                if control_register == target_register # the lookup array does not account for the communication qubits
                    push!(circuit, gate_to_apply(T, DQC_control, DQC_target ))
                    add_noise(circuit, [DQC_control, DQC_target], noise.two_q_gate_noise; two_qubits = true)
                    #num_two_qubit_gates +=1
                else

                        # Perform telegate between control and target qubit in different registers, 
                        # i.e., push!(circuit, sCNOT(comm_inv_perm_idx(control),comm_inv_perm_idx(target)) ) remotely
                        # println("Performing a telegated between (unmapped) qubits $control and $target")
                    # if telegate_overhead
                    circuit = add_telegate(circuit, DQC_control, DQC_target, control_register, target_register, n, noise)
                    telegates_layer = true
                    #     num_telegates += 1
                    # else 
                    #push!(circuit, gate_to_apply(T, n.comm_inv_perm_idx[control], n.comm_inv_perm_idx[target]))
                    #num_telegates += 1
                    
                end
            else
                throw("Circuit contains gates that have not been classified as Single- or Two-Qubit gate so far.")
            end
        end
        
        telegates_layer ? add_noise(circuit, idle_qubits_DQC, noise.idle_depolarising_noise_tele) : add_noise(circuit, idle_qubits_DQC, noise.idle_depolarising_noise)

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

    return circuit
end


function add_telegate(circuit, DQC_control, DQC_target, control_register, target_register, n, noise)

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
    
    
    # This H-CNOT only mimics the way Bell state entanglement is created; in reality, this is achieved via beam splitters. 
    # In order to account for the differing physical circumstance, we apply a two-qubit depolarising channel with a specific noise probability afterwards
    push!(circuit, sHadamard(control_comm_index))
    push!(circuit, sCNOT(control_comm_index, target_comm_index))
    add_noise(circuit, [control_comm_index, target_comm_index], noise.comm_qubit_init_noise; two_qubits = true)


    # println("Adding noise to comm qubit $control_comm_index and $target_comm_index")
    #circuit = add_noise(circuit, telegate_noise, control_comm_index)
    #circuit = add_noise(circuit, telegate_noise, target_comm_index)

    push!(circuit, sCNOT(DQC_control, control_comm_index))
    add_noise(circuit, [DQC_control, control_comm_index], noise.two_q_gate_noise_diff_species; two_qubits = true)

    # Measurement 
    pauli_string_control = build_pauli_string_measurement(n.num_qubits, [control_comm_index])
    classical_register_index_control = control_comm_index - sum(n.register_sizes[1:control_register])
    meas_control = PauliMeasurement(pauli_string_control, classical_register_index_control)
    add_noise(circuit, [control_comm_index], noise.comm_qubit_measurement_noise)
    push!(circuit, meas_control)

    # Conditional Operations
    add_noise(circuit, [target_comm_index], noise.comm_idle_depolarising_noise)
    push!(circuit, Types.ConditionalGate(sX(target_comm_index),sId1(target_comm_index), meas_control.bit)) # perform the conditional operation
    push!(circuit, Types.ConditionalGate(sX(control_comm_index),sId1(control_comm_index), meas_control.bit)) # restore the |0> state in the control comm qubit
    add_noise(circuit, [target_comm_index], noise.classical_comm_noise) # no matter if X or I applied, we assume some classical communication noise
    #meas_control.bit indicates the bit in the register that will hold the measurement information
    # if state.bits[op.controlbit] is true, the measurment yielded eigenvalue -1, if it is false, it yielded +1 
    #print(meas_control.bit)
    # Ideally, we would add noise conditional on whether or not we apply a gate. However, it is acceptable to simply assume that we apply an identity gate
    # in order to compensate for the effect for decoherence here.
    # if meas_control.bit
    #     add_noise(circuit, [target_comm_index], noise.single_comm_q_gate_noise)
    # else
    #     add_noise(circuit, [target_comm_index], noise.comm_idle_depolarising_noise)
    # end
    add_noise(circuit, [target_comm_index], noise.comm_idle_depolarising_noise)

    #The restoration can probably be neglected since in reality, the Bell pair will be created anew via photonic beam splitters. for the sake of simulation however,
    # we assume lossless restoration
    
    push!(circuit, sCNOT(target_comm_index, DQC_target ))
    add_noise(circuit, [target_comm_index, DQC_target], noise.two_q_gate_noise_diff_species; two_qubits = true)

    push!(circuit, sHadamard(target_comm_index))
    add_noise(circuit, [target_comm_index], noise.single_comm_q_gate_noise)

    
    # Measurement
    pauli_string_target = build_pauli_string_measurement(n.num_qubits, [target_comm_index])
    classical_register_index_target = target_comm_index - sum(n.register_sizes[1:target_register])
    meas_target = PauliMeasurement(pauli_string_target, classical_register_index_target)
    add_noise(circuit, [target_comm_index], noise.comm_qubit_measurement_noise)
    push!(circuit, meas_target)
    
    # Conditional Operations
    add_noise(circuit, [DQC_control], noise.idle_depolarising_noise)
    push!(circuit, Types.ConditionalGate(sZ(DQC_control),sId1(DQC_control), meas_target.bit))
    push!(circuit, Types.ConditionalGate(sX(target_comm_index),sId1(target_comm_index), meas_target.bit))  # restore the |0> state in the target comm qubit
    add_noise(circuit, [DQC_control], noise.classical_comm_noise) # no matter if X or I applied, we assume some classical communication noise
    # if meas_control.bit
    #     add_noise(circuit, [DQC_control], noise.single_q_gate_noise)
    # else
    #     add_noise(circuit, [DQC_control], noise.idle_depolarising_noise)
    # end
    add_noise(circuit, [DQC_control], noise.comm_idle_depolarising_noise)

    # Introduce noise to data qubits
    #circuit = add_noise(circuit, telegate_noise, comm_inv_perm_idx(control))
    #circuit = add_noise(circuit, telegate_noise, comm_inv_perm_idx(target))

    return circuit

end


function add_noise(circuit, qubits::Vector{Int}, prob::Float64; two_qubits = false) 
    """Depolarising noise on a set of qubits"""
    if prob<0 || prob > 1
        throw("Please provide a valid noise probabilty in [0,1]")
    end
    if two_qubits
        noise = NoiseOp(TwoQubitDepolarisingNoise(prob),qubits);
    else
        noise = NoiseOp(UnbiasedUncorrelatedNoise(prob),qubits);
    end

    push!(circuit, noise)
    return circuit
end

# function add_noise(circuit, prob::Float64) 
#     """Depolarising circuit noise on all qubits"""
#     if prob<0 || prob > 1
#         throw("Please provide a valid noise probabilty in [0,1]")
#     end
#     noise = NoiseOpAll(UnbiasedUncorrelatedNoise(prob));
#     push!(circuit, noise)
#     return circuit
# end

# function add_noise(circuit, gate::SingleQubitGate, prob::Float64) 
#     """Depolarising gate noise for single-qubit gate"""
#     if prob<0 || prob > 1
#         throw("Please provide a valid noise probabilty in [0,1]")
#     end
#     gate_noise = NoiseOp(UnbiasedUncorrelatedNoise(prob), affectedqubits(gate));
#     push!(circuit, gate_noise)
#     return circuit
# end

#abstract type AbstractNoise end

#abstract type AbstractNoiseOp <: AbstractOperation end
# Recall: noiseop takes indices and applies the noise
#NoiseOp(noise, indices::AbstractVector{Int}) = NoiseOp(noise, tuple(indices...))

struct TwoQubitDepolarisingNoise{T} <: AbstractNoise
    p::T
end
TwoQubitDepolarisingNoise(p::Integer) = TwoQubitDepolarisingNoise(float(p))

function applynoise!(s::AbstractStabilizer, noise::TwoQubitDepolarisingNoise, indices::Tuple{Int, Int})
    infid = noise.p/15
    i,j = indices[1], indices[2]
    r = rand()
    if r<infid
        apply_single_x!(s,i)
    elseif r<2infid
        apply_single_z!(s,i)
    elseif r<3infid
        apply_single_y!(s,i)
    elseif r<4infid
        apply_single_x!(s,j)
    elseif r<5infid
        apply_single_z!(s,j)
    elseif r<6infid
        apply_single_y!(s,j)
    elseif r<7infid
        apply_single_x!(s,i)
        apply_single_x!(s,j)
    elseif r<8infid
        apply_single_x!(s,i)
        apply_single_z!(s,j)
    elseif r<9infid
        apply_single_x!(s,i)
        apply_single_y!(s,j)
    elseif r<10infid
        apply_single_z!(s,i)
        apply_single_x!(s,j)
    elseif r<11infid
        apply_single_z!(s,i)
        apply_single_y!(s,j)
    elseif r<12infid
        apply_single_z!(s,i)
        apply_single_z!(s,j)
    elseif r<13infid
        apply_single_y!(s,i)
        apply_single_x!(s,j)
    elseif r<14infid
        apply_single_y!(s,i)
        apply_single_y!(s,j)
    elseif r<15infid
        apply_single_y!(s,i)
        apply_single_z!(s,j)
    end
    s
end

# function add_noise(circuit, qubits::Int[], prob::Float64) 
#     """Depolarising gate noise for two-qubit gate"""
#     if prob<0 || prob > 1
#         throw("Please provide a valid noise probabilty in [0,1]")
#     end
    
#     gate_noise = NoiseOp(TwoQubitDepolarisingNoise(prob),qubits);
#     push!(circuit, gate_noise)
#     return circuit
# end


# function add_noise(circuit, prob::Float64, qubit::Int, noise_gate::Gate)
#     # Pauli noise of specific type on single qubit
#     if prob<0 || prob > 1
#         throw("Please provide a valid noise probabilty in [0,1]")
#     end
#     gate_noise_channel = Pauli_gate_noise(typeof(noise_gate), prob)
#     gate_noise = NoiseOp(gate_noise_channel, [qubit])
#     push!(circuit, gate_noise)
#     return circuit
# end


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






function dqc_state_prep(gates, code_params, network_specs, noise)
    DQC_circuit = construct_DQC_executable_circuit(gates, network_specs, noise)
    success = verify_success(DQC_circuit, code_params.target_state, network_specs; comm_setting = true)
    println("Success: $success")
    #println("DQC circuit: \n $DQC_circuit")
   
    css_lut_decoder = CSSTableDecoder(code_params.qec_code, error_weight = 1)
    println("X Lookup Table: $(css_lut_decoder.tabledecoderx.lookup_table)")
    println("Z Lookup Table: $(css_lut_decoder.tabledecoderz.lookup_table)")
    H = parity_checks(css_lut_decoder)
    # Build perfect syndrome circuit, with ancillas being appended to DQC cores (consisting of data and comm qubits)
    # and bit string measurements being written into classical register n.num_registers*(n.num_registers-1) + 1
    syndrome_circ, n_anc_syndrome, syndrome_bits = QECTools.syndrome_circuit(H, network_specs.num_qubits + 1, network_specs.num_registers*(network_specs.num_registers-1)+1, network_specs)

    logical_Z_circ, n_anc_logical_Z, logical_Z_bits = QECTools.syndrome_circuit(code_params.logical_Zs, network_specs.num_qubits + n_anc_syndrome + 1, last(syndrome_bits)+1, network_specs )

    total_number_qubits = network_specs.num_qubits + n_anc_syndrome + n_anc_logical_Z

    circuit = vcat(DQC_circuit, syndrome_circ, logical_Z_circ)
    logical_failures = 0
    #println("\n\n\n $circuit \n\n\n")
    for _ in 1:network_specs.num_shots
        initial_state = Register(one(MixedDestabilizer,total_number_qubits),network_specs.num_registers*(network_specs.num_registers-1)+length(syndrome_bits)+length(logical_Z_bits))
        st, stat = mctrajectory!(initial_state, copy(circuit))
        
        matching = compare_states(st.stab, code_params.target_state, network_specs)
        println("Matching States: $matching")
        syndrome = st.bits[syndrome_bits]
        measured_logical_Z_bits = st.bits[logical_Z_bits] 

        println(st.bits)
        println(syndrome)
        println(measured_logical_Z_bits)

        error_guess = decode(css_lut_decoder, syndrome)
        println(error_guess)
        if isnothing(error_guess)
            # Table decoder can't find a matching syndrome bc error weight is too high
            logical_failures += 1
            continue
        end

        #TODO: Apply the correction gate (since there is a correctable syndrome, this measn that the correction will bring us back to the codespace)
            # applying pauli correction: https://github.com/QuantumSavory/QuantumClifford.jl/blob/444f341a50d2926b16b63d98586b8b06a7b6ac10/src/ecc/decoder_correction_gate.jl maps back to the codespace so stabilisers are satisfied. 
        #and see if the logical syndrome is correct now (if everything is zero we are in the codespace and the logicals are satisfied)
        DecoderCorrectionGate(css_lut_decoder, network_specs.data_qubits,syndrome )

        # initial_state = Register(one(MixedDestabilizer,total_number_qubits),network_specs.num_registers*(network_specs.num_registers-1)+length(syndrome_bits)+length(logical_Z_bits))
        # st, stat = mctrajectory!(initial_state, copy(circuit))
        
        # matching = compare_states(st.stab, code_params.target_state, network_specs)
        # println("Matching States: $matching")
        # syndrome = st.bits[syndrome_bits]
        # measured_logical_Z_bits = st.bits[logical_Z_bits] 

        # println(st.bits)
        # println(syndrome)
        # println(measured_logical_Z_bits)
    end

    return nothing
end


end
