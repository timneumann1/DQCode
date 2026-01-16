# simulation.jl

#using .Types
using .Types: SimulationParameters, SimulationFidelity, Circuit, HadamardGate, IdentityGate, PauliXGate, PauliYGate, PauliZGate, CNOT_Gate, Gate

# For convenient graph data structures
using Graphs

# For discrete event simulation
using ResumableFunctions
using ConcurrentSim

# The workhorse for the simulation
using QuantumSavory
using QuantumSavory.ProtocolZoo: EntanglerProt
using QuantumSavory: H, CNOT

const bell = StabilizerState("XX ZZ")
const noisy_pair_func(F) = (1-F)*MixedState(bell) + F*projector(bell)

gate_to_apply(::Type{Types.HadamardGate}) = H


@resumable function single_qubit_gate(sim, network, register, qubit, gate)
    # Accepts an array of single-qubit operations and executes them in parallel (during on time slice)
    #num_ops = length(register_array)
    #@assert num_ops === length(qubit_array) === length(ops_array)
    print("Single qubit gate:")
    print("\n $register, $qubit, $gate")
    @yield request(network[register, qubit])
    print("HEEEERE")
    print("Executing $gate gate on ($register, $qubit)")
    print(network[register,qubit])
    apply!(network[register,qubit], gate)
    @simlog sim "Executing $gate gate on ($register, $qubit)"
    @yield unlock(network[register, qubit])
    print("HEEEERE")
end

@resumable function CNOT_gate(sim, network, register, control_qubit, target_qubit)
    # Accepts an array of control and target qubits for CNOT applications
    # that can be executed in parallel (during on time slice) within one core
    @yield request(network[register, control_qubit])
    @yield request(network[register, target_qubit])
    apply!( ( network[register, control_qubit], network[register, target_qubit]), CNOT)
    @simlog sim "Executing CNOT between ($register, $control_qubit) and ($register, $target_qubit)"
    @yield unlock(network[register, control_qubit])
    @yield unlock(network[register, target_qubit])    
end

@resumable function CNOT_telegate(sim, network, control_register, target_register, control_qubit, target_qubit, params)
    # Accepts control and target qubit for remote CNOT application across cores
    @assert network[control_register, 1] !== nothing "Missing flying qubit at control register"
    @assert network[target_register, 1] !== nothing "Missing flying qubit at target register"
    @assert network[control_register, control_qubit] !== nothing "Missing control qubit"
    @assert network[target_register, target_qubit] !== nothing "Missing target qubit"
    @assert control_qubit !== 1 "Control qubit is not a memory qubit"
    @assert target_qubit !== 1 "Target qubit is not a memory qubit"
        
    # Create Bell Pair between first qubits of registers i and j
    noisy_pair = noisy_pair_func(params.bell_state_fidelity)       # Retrieving a (noisy) Bell pair with fidelity F

    eprot = EntanglerProt(sim, network, control_register, target_register; pairstate=noisy_pair, chooseslotA=1, chooseslotB=1, rounds=-1, attempts = -1, success_prob=params.success_prob, attempt_time = params.attempt_time) 
    @process eprot()
    #@yield timeout(sim, entangler_busy_time)
    
    # Perform CNOT telegate between qubits q1 and q2 (TODO: extend to application of more gates, cf. EJPP protocol)
    @yield request(network[control_register,control_qubit]) 
    @yield request(network[target_register,target_qubit])
    @yield request(network[control_register,1]) 
    @yield request(network[target_register,1])

    @simlog sim "Executing Telegate-CNOT between ($control_register, $control_qubit) and ($target_register, $target_qubit)"
    apply!((network[control_register,control_qubit], network[control_register,1]), CNOT)  
    @yield timeout(sim, params.two_qubit_gate_exec_time)
    m1 = project_traceout!(network[control_register,1], Z)
    @yield timeout(sim, params.projective_measurement_time)
    # If the result is |1> (m1 == 2), flip flying qubit 1 in target_register
    # TODO: Add communication channel
    if m1 == 2
        apply!(network[target_register,1], X)
        @yield timeout(sim, params.classical_communication_time)
        @yield timeout(sim, params.single_qubit_gate_exec_time)
    end
    apply!((network[target_register,1], network[target_register,target_qubit]), CNOT)
    @yield timeout(sim, params.two_qubit_gate_exec_time)
    apply!(network[target_register,1], H)
    @yield timeout(sim, params.single_qubit_gate_exec_time)
    m2 = project_traceout!(network[target_register,1], Z)
    @yield timeout(sim, params.projective_measurement_time)
    # If the result is '1' (m1 == 2), apply Z to q1 in ri
    if m2 == 2
        apply!(network[control_register,control_qubit], Z)
        @yield timeout(sim, params.single_qubit_gate_exec_time)
    end
    @yield unlock(network[control_register,control_qubit])
    @yield unlock(network[target_register,target_qubit])
    @yield unlock(network[control_register,1])
    @yield unlock(network[target_register,1])
end

#function gate_to_apply(::Type{HadamardGate}) = H

# function gate_to_apply(g::Gate)
#     if g isa HadamardGate
#         return H
#     elseif g isa PauliXGate
#         return X
#     elseif g isa PauliYGate
#         return Y
#     elseif g isa PauliZGate
#         return Z
#     # elseif g isa IdentityGate
#     #     return nothing  # identity, skip
#     # elseif g isa CNOTGate
#     #     return CNOT  # or whatever simulator CNOT expects
#     # else
#     #     error("Unknown gate")
#     end
# end


#     elseif g isa PauliXGate
#         return X
#     elseif g isa PauliYGate
#         return Y
#     elseif g isa PauliZGate
#         return Z

@resumable function single_qubit_gates_per_layer(sim, network, single_qubit_gates, params, register_lookup_array)
    print("In the function for single qubit gate*s now.")
    print(single_qubit_gates)
    for gtype in keys(single_qubit_gates)
        indices = get!(single_qubit_gates, gtype, Int[])   
        print(indices) 
        print(gtype)    
        for index in indices
            print(register_lookup_array[index], index, gate_to_apply(gtype))
            print("Initialising a process of single qubit gate for $gate_to_apply(gtype) on index $index")
            @process single_qubit_gate(sim, network, register_lookup_array[index], index, gate_to_apply(gtype))#gate_to_apply(gtype))
            print("Made it here")
        end
    end
    @yield timeout(sim, params.single_qubit_gate_exec_time)
end

@resumable function CNOT_gates_per_layer(sim, network, CNOT_gates, params, register_lookup_array)
    for (control, target) in CNOT_gates
        @process CNOT_gate(sim, network, register_lookup_array[control], control, target)
    end
    @yield timeout(sim, params.single_qubit_gate_exec_time)

end

#TODO: potentially refactor back to simple function (if single telegate is kept as assumption)
@resumable function CNOT_telegates_per_layer(sim, network, CNOT_telegates, params, register_lookup_array)
    for (control, target) in CNOT_telegates
        @process CNOT_telegate(sim, network, register_lookup_array[control], register_lookup_array[target], control, target, params)
    end
end

# @resumable function steane_encoding_circuit(sim, network, params)
#     # Initialize memory qubits in zero state
#     initialize!([network[i,j] for i in 1:length(params.register_sizes) for j in 2:params.register_sizes[i]], SProjector(Z1⊗Z1⊗Z1⊗Z1⊗Z1⊗Z1⊗Z1)  ) 

#     # TODO: For each layer, identify the longest execution time among parallel operations
#     # make adding gates possible
#     @yield @process single_qubit_gate(sim, network, 1, 2, H, params)
#     @yield @process single_qubit_gate(sim, network, 1, 3, H, params)
#     @yield @process single_qubit_gate(sim, network, 1, 5, H, params)

#     @yield @process CNOT_gate(sim, network, 1, 2, 4, params)
#     @yield @process telegate_CNOT(sim, network, 1, 2, 5, 2, params)
#     @yield @process telegate_CNOT(sim, network, 1, 2, 3, 3, params)
#     @yield @process CNOT_gate(sim, network, 1, 2, 6, params)
#     @yield @process telegate_CNOT(sim, network, 1, 2, 5, 3, params)
#     @yield @process telegate_CNOT(sim, network, 1, 2, 3, 2, params)
#     @yield @process telegate_CNOT(sim, network, 1, 2, 2, 3, params)
#     @yield @process CNOT_gate(sim, network, 1, 3, 4, params)
#     @yield @process CNOT_gate(sim, network, 1, 5, 6, params)
# end

function initialise_simulation(params)
    # Initialise the network
    sizes = params.register_sizes
    R = length(sizes) # Number of registers
    # All of the quantum register we will be simulating
    registers = Register[]
    for s in sizes
        traits = [Qubit() for _ in 1:s]
        repr = [QuantumOpticsRepr() for _ in 1:s]  # recast to Clifford
        bg = [T2Dephasing(params.T2_dephasing) for _ in 1:s]          # define other noise on the registers?
        push!(registers, Register(traits,repr, bg)) 
    end

    graph = grid([R])
    network = RegisterNet(graph, registers) # A graphs with extra "meta data"

    # The scheduler datastructure for the discrete event simulation
    sim = get_time_tracker(network)
    initialize!([network[i,j] for i in eachindex(sizes) for j in 2:sizes[i]], SProjector(Z1⊗Z1⊗Z1⊗Z1⊗Z1⊗Z1⊗Z1)  ) 

    return sim, network
end

function run_simulation(params::Types.SimulationParameters, circuit::Types.Circuit, register_lookup_array::Array{Int})::Types.SimulationFidelity
    
    sim, network = initialise_simulation(params)
    
    circuit_matrix = circuit.gates
    for col in axes(circuit_matrix, 2) # each column corresponds to one layer
        single_qubit_gates = Dict{Type, Vector{Int}}() # stores single qubit gates and corresponding qubit indices per layer
        CNOT_gates = Vector{Tuple{Int, Int}}()
        CNOT_telegates = Vector{Tuple{Int, Int}}()
        flags = Set{Int}() # to flag the control/target qubits that can be ignored

        for row in axes(circuit_matrix, 1) # each row corresponds to one qubit
            gate = circuit_matrix[row,col]
            #print(gate)
            if gate isa Union{PauliXGate, PauliYGate, PauliZGate, HadamardGate} #  Gate && !(gate isa CNOT) 
                gate_type = typeof(gate)
                push!(get!(single_qubit_gates, gate_type, Int[]), row)
                #single_qubit_gates add row (value) to circuit_matrix.gates[row, column] (gate) in the dict
            elseif gate isa CNOTGate
                row in flags && continue

                control = gate.control
                target = gate.target
                if register_lookup_array[control] == register_lookup_array[target]
                    push!(CNOT_gates, (control, target) )
                    #add (control, target) to CNOT_gates array
                else 
                    push!(CNOT_telegates, (control, target) )
                    #add (control, target) to CNOT_telegates
                end
                push!(flags, control)
                push!(flags, target)
            end
        end
        println("Constructed my matrix: $circuit_matrix")
        @process single_qubit_gates_per_layer(sim, network, single_qubit_gates, params, register_lookup_array)
        println("Here now")
        #@yield @process CNOT_gates(sim, network, CNOT_gates, params, register_lookup_array)
        #@yield @process CNOT_telegates(sim, network, CNOT_telegates, params, register_lookup_array)

    end

    #TODO: per slice, add check that all cnots are correct

    #@process simulate_circuit(sim, network, circuit, params)

    run(sim, params.simulation_time)
    steane_7_state = StabilizerState("ZIZIZIZ XIXIXIX IZZIIZZ IXXIIXX IIIZZZZ IIIXXXX ZZZZZZZ") 

    fidelity = real(observable(vcat( [network[1,i] for i in 2:6],[network[2,j] for j in 2:3]), SProjector(steane_7_state)))
    return Types.SimulationFidelity(fidelity)
end
