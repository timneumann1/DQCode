# simulation.jl

using ..Types

# For convenient graph data structures
using Graphs

# For discrete event simulation
using ResumableFunctions
using ConcurrentSim

# The workhorse for the simulation
using QuantumSavory
using QuantumSavory.ProtocolZoo: EntanglerProt

const bell = StabilizerState("XX ZZ")
const noisy_pair_func(F) = (1-F)*MixedState(bell) + F*projector(bell)

@resumable function single_qubit_gate(sim, network, register, qubit, gate, params)
    # Accepts an array of single-qubit operations and executes them in parallel (during on time slice)
    #num_ops = length(register_array)
    #@assert num_ops === length(qubit_array) === length(ops_array)
    @yield request(network[register, qubit])
    apply!(network[register,qubit], gate)
    @simlog sim "Executing $gate gate on ($register, $qubit)"
    @yield timeout(sim, params.single_qubit_gate_exec_time)
    @yield unlock(network[register, qubit])
    
end

@resumable function CNOT_gate(sim, network, register, control_qubit, target_qubit, params)
    # Accepts an array of control and target qubits for CNOT applications
    # that can be executed in parallel (during on time slice) within one core
    @yield request(network[register, control_qubit])
    @yield request(network[register, target_qubit])
    apply!( ( network[register, control_qubit], network[register, target_qubit]), CNOT)
    @yield timeout(sim, params.single_qubit_gate_exec_time)
    @simlog sim "Executing CNOT between ($register, $control_qubit) and ($register, $target_qubit)"
    @yield unlock(network[register, control_qubit])
    @yield unlock(network[register, target_qubit])    
end

@resumable function telegate_CNOT(sim, network, control_register, target_register, control_qubit, target_qubit, params)
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
    unlock(network[control_register,control_qubit])
    unlock(network[target_register,target_qubit])
    unlock(network[control_register,1])
    unlock(network[target_register,1])
end

@resumable function steane_encoding_circuit(sim, network, params)
    # Initialize memory qubits in zero state
    initialize!([network[i,j] for i in 1:length(params.register_sizes) for j in 2:params.register_sizes[i]], SProjector(Z1⊗Z1⊗Z1⊗Z1⊗Z1⊗Z1⊗Z1)  ) 

    # TODO: For each layer, identify the longest execution time among parallel operations
    # make adding gates possible
    @yield @process single_qubit_gate(sim, network, 1, 2, H, params)
    @yield @process single_qubit_gate(sim, network, 1, 3, H, params)
    @yield @process single_qubit_gate(sim, network, 1, 5, H, params)

    @yield @process CNOT_gate(sim, network, 1, 2, 4, params)
    @yield @process telegate_CNOT(sim, network, 1, 2, 5, 2, params)
    @yield @process telegate_CNOT(sim, network, 1, 2, 3, 3, params)
    @yield @process CNOT_gate(sim, network, 1, 2, 6, params)
    @yield @process telegate_CNOT(sim, network, 1, 2, 5, 3, params)
    @yield @process telegate_CNOT(sim, network, 1, 2, 3, 2, params)
    @yield @process telegate_CNOT(sim, network, 1, 2, 2, 3, params)
    @yield @process CNOT_gate(sim, network, 1, 3, 4, params)
    @yield @process CNOT_gate(sim, network, 1, 5, 6, params)
end

function run_simulation(params::Types.SimulationParameters)::Types.SimulationFidelity
    
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

    #TODO: pass the circuit object here
    @process steane_encoding_circuit(sim, network, params)

    run(sim, params.simulation_time)
    steane_7_state = StabilizerState("ZIZIZIZ XIXIXIX IZZIIZZ IXXIIXX IIIZZZZ IIIXXXX ZZZZZZZ") 

    fidelity = real(observable(vcat( [network[1,i] for i in 2:6],[network[2,j] for j in 2:3]), SProjector(steane_7_state)))
    return Types.SimulationFidelity(fidelity)
end
