# simulation.jl
module Simulation

using ..Types

using Graphs
using ResumableFunctions
using ConcurrentSim
using QuantumSavory
using QuantumSavory.ProtocolZoo: EntanglerProt
using QuantumSavory: H, CNOT, X, Y, Z

using QuantumClifford: MixedDestabilizer, sHadamard, sCNOT, @S_str

export run_simulation

# For QuantumOpticsRepr
const bell_pair = StabilizerState("XX ZZ")
const noisy_pair_func(F) = (1-F)*MixedState(bell_pair) + F*projector(bell_pair)

# For CliffordRepr
const bell_pair_tableau = S"XX ZZ"
const bell_pair_stabilizer = StabilizerState(bell_pair_tableau)
const bell_pair_stabilizer_density = SProjector(bell_pair_stabilizer)
const noisy_pair_func_stabilizer(F) = (1-F)*MixedState(bell_pair)+ F*bell_pair_stabilizer_density

# gate_to_apply(::Type{HadamardGate}) = H  # QuantumOpticsRepr
gate_to_apply(::Type{HadamardGate}) = sHadamard  # CliffordRepr
gate_to_apply(::Type{PauliXGate}) = X
gate_to_apply(::Type{PauliYGate}) = Y
gate_to_apply(::Type{PauliZGate}) = Z

@resumable function single_qubit_gate(sim, network, register, qubit, gate, single_qubit_gates_layer_executed, num_single_qubit_gates_processes)
    
    @yield request(network[register, qubit])
    apply!(network[register,qubit], gate)
    @simlog sim "Executing $gate gate on ($register, $qubit)"
    @yield unlock(network[register, qubit])

    num_single_qubit_gates_processes[]-=1
    @simlog sim "num_single_qubit_gates_processes remaining: $(num_single_qubit_gates_processes[])"   
    num_single_qubit_gates_processes[]==0 && succeed(single_qubit_gates_layer_executed)

end

@resumable function CNOT_gate(sim, network, register, control_qubit, target_qubit, num_CNOT_gate_processes, CNOT_gate_layer_executed)
   
    @yield request(network[register, control_qubit])
    @yield request(network[register, target_qubit])
    apply!( ( network[register, control_qubit], network[register, target_qubit]), CNOT)
    @simlog sim "Executing CNOT between ($register, $control_qubit) and ($register, $target_qubit)"
    @yield unlock(network[register, control_qubit])
    @yield unlock(network[register, target_qubit])  
    
    num_CNOT_gate_processes[]-=1
    @simlog sim "num CNOT_gate_processes remaining: $(num_CNOT_gate_processes[])"
    num_CNOT_gate_processes[]==0 && succeed(CNOT_gate_layer_executed)
end

@resumable function CNOT_telegate(sim, network, control_register, target_register, control_qubit, target_qubit, params, num_CNOT_telegate_processes, CNOT_telegate_layer_executed)
    
    @assert network[control_register, 1] !== nothing "Missing flying qubit at control register"
    @assert network[target_register, 1] !== nothing "Missing flying qubit at target register"
    @assert network[control_register, control_qubit] !== nothing "Missing control qubit"
    @assert network[target_register, target_qubit] !== nothing "Missing target qubit"
    @assert control_qubit !== 1 "Control qubit is not a memory qubit"
    @assert target_qubit !== 1 "Target qubit is not a memory qubit"
        
    # Create Bell Pair between first qubits of registers i and j
    noisy_pair = noisy_pair_func_stabilizer(params.bell_state_fidelity)  # noisy_pair_func for QuantumOptics     # Retrieving a (noisy) Bell pair with fidelity F
    eprot = EntanglerProt(sim, network, control_register, target_register; pairstate=noisy_pair, chooseslotA=1, chooseslotB=1, rounds=1, attempts = -1, success_prob=params.success_prob, attempt_time = params.attempt_time) 
    @yield @process eprot() # @yield is redundant since the request depends on the process to finish already
    
    # Perform CNOT telegate between qubits q1 and q2 (TODO: extend to application of more gates, cf. EJPP protocol)
    @yield request(network[control_register,control_qubit]) 
    @yield request(network[target_register,target_qubit])
    @yield request(network[control_register,1]) 
    @yield request(network[target_register,1])

    @simlog sim "Executing Telegate-CNOT between ($control_register, $control_qubit) and ($target_register, $target_qubit)"
    @simlog sim network[1]
    apply!((network[control_register,control_qubit], network[control_register,1]), CNOT)  

    @yield timeout(sim, params.two_qubit_gate_exec_time)
    m1 = project_traceout!(network[control_register,1], Z)
    @yield timeout(sim, params.projective_measurement_time)
    # If the result is |1> (m1 == 2), flip flying qubit 1 in target_register
    # TODO: Add communication channel
    if m1 == 2
        apply!(network[target_register,1], X)
        @yield timeout(sim, params.single_qubit_gate_exec_time)
    end
    @simlog sim "Measured m1 = $m1"
    @yield timeout(sim, params.classical_communication_time)
    apply!((network[target_register,1], network[target_register,target_qubit]), CNOT)
    @yield timeout(sim, params.two_qubit_gate_exec_time)
    apply!(network[target_register,1], sHadamard) # only H for QuantumOpticsRepr
    @yield timeout(sim, params.single_qubit_gate_exec_time)
    m2 = project_traceout!(network[target_register,1], Z)
    @yield timeout(sim, params.projective_measurement_time)
    # If the result is '1' (m1 == 2), apply Z to q1 in ri
    if m2 == 2
        apply!(network[control_register,control_qubit], Z)
        @yield timeout(sim, params.single_qubit_gate_exec_time)
    end
    @simlog sim "Measured m2 = $m2"
    @yield timeout(sim, params.classical_communication_time)
    @yield unlock(network[control_register,control_qubit])
    @yield unlock(network[target_register,target_qubit])
    @yield unlock(network[control_register,1])
    @yield unlock(network[target_register,1])
    
    num_CNOT_telegate_processes[]-=1
    @simlog sim "num CNOT_telegate_processes remaining: $(num_CNOT_telegate_processes[])"
    num_CNOT_telegate_processes[]==0 && succeed(CNOT_telegate_layer_executed)
end

@resumable function single_qubit_gates_per_layer(sim, network, single_qubit_gates, params, register_lookup_array, layer_executed, num_processes)
    
    single_qubit_gates_layer_executed = Event(sim)
    num_single_qubit_gates_processes = isempty(single_qubit_gates) ? Ref(0) : Ref(sum(length, values(single_qubit_gates)))#@isdefined(single_qubit_gates) ? Ref() : Ref(0)

    for gtype in keys(single_qubit_gates)
        indices = get!(single_qubit_gates, gtype, Int[])     
        for index in indices
            register = register_lookup_array[index]
            index_in_register = index - sum(params.register_sizes[1:(register-1)])  # need to account for index of qubit within register
            @process single_qubit_gate(sim, network, register, index_in_register, gate_to_apply(gtype), single_qubit_gates_layer_executed, num_single_qubit_gates_processes)#gate_to_apply(gtype))
        end
    end
    if num_single_qubit_gates_processes[] != 0
        @yield single_qubit_gates_layer_executed 
        @yield timeout(sim, params.single_qubit_gate_exec_time)
    end
    
    num_processes[]-=1
    @simlog sim "num layer_executed processes remaining: $(num_processes[])"
    num_processes[]==0 && succeed(layer_executed)
end

@resumable function CNOT_gates_per_layer(sim, network, CNOT_gates, params, register_lookup_array, layer_executed, num_processes)
    
    CNOT_gate_layer_executed = Event(sim)
    num_CNOT_gate_processes = isempty(CNOT_gates) ? Ref(0) : Ref(length(CNOT_gates)) 

    for (control, target) in CNOT_gates
        register = register_lookup_array[control]
        control_index_in_register = control - sum(params.register_sizes[1:(register-1)])
        target_index_in_register = target - sum(params.register_sizes[1:(register-1)])
        @process CNOT_gate(sim, network, register, control_index_in_register, target_index_in_register, num_CNOT_gate_processes, CNOT_gate_layer_executed)
    end

    if num_CNOT_gate_processes[] != 0
        @yield CNOT_gate_layer_executed
        @yield timeout(sim, params.two_qubit_gate_exec_time)
    end
    
    num_processes[]-=1
    @simlog sim "num layer_executed processes remaining: $(num_processes[])"
    num_processes[]==0 && succeed(layer_executed)
end

#TODO: potentially refactor back to simple function (if single telegate is kept as assumption)
@resumable function CNOT_telegates_per_layer(sim, network, CNOT_telegates, params, register_lookup_array, layer_executed, num_processes)

    CNOT_telegate_layer_executed = Event(sim)
    num_CNOT_telegate_processes = isempty(CNOT_telegates) ? Ref(0) : Ref(length(CNOT_telegates))

    for (control, target) in CNOT_telegates
        control_register = register_lookup_array[control]
        target_register = register_lookup_array[target]
        control_index_in_register = control - sum(params.register_sizes[1:(control_register-1)])
        target_index_in_register = target - sum(params.register_sizes[1:(target_register-1)])
        @process CNOT_telegate(sim, network, control_register, target_register, control_index_in_register, target_index_in_register, params, num_CNOT_telegate_processes, CNOT_telegate_layer_executed)
    end

    if num_CNOT_telegate_processes[] != 0
        @yield CNOT_telegate_layer_executed
    end

    num_processes[]-=1
    @simlog sim "num layer_executed processes remaining: $(num_processes[])"
    num_processes[]==0 && succeed(layer_executed)
end

function initialise_simulation(params)

    # Initialise the network
    sizes = params.register_sizes
    R = length(sizes) # Number of registers
    registers = Register[]
    
    for s in sizes
        traits = [Qubit() for _ in 1:s]
        repr = [CliffordRepr() for _ in 1:s]        #QuantumOpticsRepr  #CliffordRepr      # TODO: recast to Clifford
        bg = [T2Dephasing(params.T2_dephasing) for _ in 1:s]    # TODO: define other noise on the registers?
        push!(registers, Register(traits,repr, bg)) 
    end

    # important: Currently, to switch between Clifford and QuantumOptics representation, must only change repr, the entangler prot argument for Bell states, and sHadamard

    graph = grid([R])
    network = RegisterNet(graph, registers)
    sim = get_time_tracker(network) # the scheduler datastructure for the discrete event simulation, creates Simulation() for all registers

    data_qubits = [network[i,j] for i in eachindex(sizes) for j in 2:sizes[i]]
    #initialize!(data_qubits, SProjector(Z1⊗Z1⊗Z1⊗Z1⊗Z1⊗Z1⊗Z1)  )   # QuantumOpticsRepr() formulation
    initialize!(data_qubits, one(MixedDestabilizer,7) )             # CliffordRepr() formulation
    
    return sim, network, data_qubits
end

@resumable function execute_layer(sim, network, single_qubit_gates, CNOT_gates, CNOT_telegates, params, register_lookup_array)
    layer_executed = Event(sim)
    num_processes = Ref(3)
    
    @process single_qubit_gates_per_layer(sim, network, single_qubit_gates, params, register_lookup_array, layer_executed, num_processes)
    @process CNOT_gates_per_layer(sim, network, CNOT_gates, params, register_lookup_array, layer_executed, num_processes)#, layer_executed, num_processes)
    @process CNOT_telegates_per_layer(sim, network, CNOT_telegates, params, register_lookup_array, layer_executed, num_processes)
    
    @yield layer_executed
end

@resumable function build_simulation_process(sim, network, params, circuit_matrix, register_lookup_array)
    
    for col in axes(circuit_matrix, 2) # each column corresponds to one layer
        @simlog sim "Entered the iteration $col in the outer loop"
        single_qubit_gates = Dict{Type, Vector{Int}}() # stores single qubit gates and corresponding qubit indices per layer
        CNOT_gates = Vector{Tuple{Int, Int}}()
        CNOT_telegates = Vector{Tuple{Int, Int}}()
        flags = Set{Int}() # to flag the control/target qubits that can be ignored

        for row in axes(circuit_matrix, 1) # each row corresponds to one qubit
            gate = circuit_matrix[row,col]
            if gate isa Union{PauliXGate, PauliYGate, PauliZGate, HadamardGate} #  Gate && !(gate isa CNOT) 
                gate_type = typeof(gate)
                push!(get!(single_qubit_gates, gate_type, Int[]), row)
            elseif gate isa CNOT_Gate
                row in flags && continue
                control = gate.control
                target = gate.target
                if register_lookup_array[control] == register_lookup_array[target]
                    push!(CNOT_gates, (control, target) )
                else 
                    push!(CNOT_telegates, (control, target) )
                end
                push!(flags, control)
                push!(flags, target)
            end
        end
        #TODO: per slice, add check that all cnots are correct
        @simlog sim "Now yielding the process for column $col: $single_qubit_gates, $CNOT_gates, $CNOT_telegates"
        @yield @process execute_layer(sim, network, single_qubit_gates, CNOT_gates, CNOT_telegates, params, register_lookup_array)

    end
end

function run_simulation(params::SimulationParameters, circuit::Circuit, register_lookup_array::Array{Int})::SimulationFidelity
    
    sim, network, data_qubits = initialise_simulation(params)
    circuit_matrix = circuit.gates

    print("Initial state is $(network[1]) \n $(network[2])")
    print("Second qubit in first register is $(network[1,2])")

    execute = @process build_simulation_process(sim, network, params, circuit_matrix, register_lookup_array)
    run(sim, execute)

    steane_7_state = StabilizerState("ZIZIZIZ XIXIXIX IZZIIZZ IXXIIXX IIIZZZZ IIIXXXX ZZZZZZZ") 
    
    print("Final state is $(network[1,2]) \n $(network[2])")
    
    fidelity = real(observable(data_qubits, SProjector(steane_7_state)))
    # TODO: recast the fidelity as the Tableau distance to the canoncial form of the stabilizer tableau

    return SimulationFidelity(fidelity)
end


end
