module Evolutionary

#include("helper_functions.jl") 
# For convenient graph data structures
using Graphs

# For discrete event simulation
using ResumableFunctions
using ConcurrentSim

# Useful for interactive work
# Enables automatic re-compilation of modified codes
using Revise

# The workhorse for the simulation
using QuantumSavory

using QuantumSavory.ProtocolZoo: EntanglerProt


# System Parameters
sizes = [6,3]        # Number of qubits in each register
T1 = 12000.00         # T1 relaxation time of all qubits   [Wang]
T2 = 4200.0         # T2 dephasing time of all qubits     [Wang]
F = 1              # Fidelity of the raw Bell pairs
t = 75

# TODO: Define units and insert realistic values
#entangler_wait_time = 0.1  # How long to wait if all qubits are busy before retrying entangling
#entangler_busy_time = 1.0  # How long it takes to establish a newly entangled pair
#swapper_wait_time = 0.1    # How long to wait if all qubits are unavailable for swapping
#swapper_busy_time = 0.15   # How long it takes to swap two qubits
#purifier_wait_time = 0.15  # How long to wait if there are no pairs to be purified
#purifier_busy_time = 0.2   # How long the purification circuit takes to execute
single_qubit_gate_exec_time = 20*10^(-6)  # Execution Time of a single qubit gate
two_qubit_gate_exec_time = 200*10^(-6)     # Execution Time of a two qubit gate
projective_measurement_time = 10^(-5)     # Time to peform a measurement
init_time = 1                         # Time to initialise the system (e.g., in the all-zero state)
classical_communication_time = 0.001  # define in relation to speed of photons in fiber cable

#TODO: Add state preparation fidelity and single-shot readout of 99.93% [Harty], single-qubit gate fidelity of 99.99916%, two-qubit fidelity of 99.97% [Löschnauer]


# characteristic_time = 1000
# p = 1-exp(-1/characteristic_time) # define probability for Pauli Noise application (Poisson point process)

#TODO: Define single qubit error rate

bell = StabilizerState("XX ZZ")
noisy_pair_func(F) = (1-F)*MixedState(bell) + F*projector(bell)
noisy_pair = noisy_pair_func(0.9689)       # Retrieving a (noisy) Bell pair with fidelity F

@resumable function single_qubit_gate(sim, network, register, qubit, gate)
    # Accepts an array of single-qubit operations and executes them in parallel (during on time slice)
    #num_ops = length(register_array)
    #@assert num_ops === length(qubit_array) === length(ops_array)
    @yield request(network[register, qubit])
    apply!(network[register,qubit], gate)
    @yield timeout(sim, single_qubit_gate_exec_time)
    @yield unlock(network[register, qubit])
    
end

@resumable function CNOT_gate(sim, network, register, control_qubit, target_qubit)
    # Accepts an array of control and target qubits for CNOT applications
    # that can be executed in parallel (during on time slice) within one core
    @yield request(network[register, control_qubit])
    @yield request(network[register, target_qubit])
    apply!( ( network[register, control_qubit], network[register, target_qubit]), CNOT)
    @yield timeout(sim, single_qubit_gate_exec_time)
    @yield unlock(network[register, control_qubit])
    @yield unlock(network[register, target_qubit])

    #num_ops = length(register_array)
    # @assert num_ops === length(control_qubits) === length(target_qubits)
    # for i in 1:num_ops
    #     r_i = register_array[i]
    #     control_i = control_qubits[i]
    #     target_i = target_qubits[i]
    #     @yield request(network[r_i,control_i]) 
    #     @yield request(network[r_i,target_i]) 
    #     apply!( (network[r_i,control_i], network[r_i,target_i]), CNOT)
    #     unlock(network[r_i, control_i])
    #     unlock(network[r_i, target_i])
    # end
    
end

@resumable function telegate_CNOT(sim, network, control_register, target_register, control_qubit, target_qubit)
    # Accepts control and target qubit for remote CNOT application across cores
    @assert network[control_register, 1] !== nothing "Missing flying qubit at control register"
    @assert network[target_register, 1] !== nothing "Missing flying qubit at target register"
    @assert network[control_register, control_qubit] !== nothing "Missing control qubit"
    @assert network[target_register, target_qubit] !== nothing "Missing target qubit"
    @assert control_qubit !== 1 "Control qubit is not a memory qubit"
    @assert target_qubit !== 1 "Target qubit is not a memory qubit"
        
    # Create Bell Pair between first qubits of registers i and j
    eprot = EntanglerProt(sim, network, control_register, target_register; pairstate=noisy_pair, chooseslotA=1, chooseslotB=1, rounds=-1, attempts = -1, success_prob=1.41*10^(-4), attempt_time = 1.168*10^(-9)) # from [Main, 2025]
    @process eprot()
    #@yield timeout(sim, entangler_busy_time)
    
    # Perform CNOT telegate between qubits q1 and q2 (TODO: extend to application of more gates, cf. EJPP protocol)
    @yield request(network[control_register,control_qubit]) 
    @yield request(network[target_register,target_qubit])
    @yield request(network[control_register,1]) 
    @yield request(network[target_register,1])

    @simlog sim "Executing Telegate-CNOT between ($control_register, $control_qubit) and ($target_register, $target_qubit)"
    apply!((network[control_register,control_qubit], network[control_register,1]), CNOT)  
    @yield timeout(sim, two_qubit_gate_exec_time)
    m1 = project_traceout!(network[control_register,1], Z)
    @yield timeout(sim, projective_measurement_time)
    # If the result is |1> (m1 == 2), flip flying qubit 1 in target_register
    # TODO: Add communication channel
    if m1 == 2
        apply!(network[target_register,1], X)
        @yield timeout(sim, classical_communication_time)
        @yield timeout(sim, single_qubit_gate_exec_time)
    end
    apply!((network[target_register,1], network[target_register,target_qubit]), CNOT)
    @yield timeout(sim, two_qubit_gate_exec_time)
    apply!(network[target_register,1], H)
    @yield timeout(sim, single_qubit_gate_exec_time)
    m2 = project_traceout!(network[target_register,1], Z)
    @yield timeout(sim, projective_measurement_time)
    # If the result is '1' (m1 == 2), apply Z to q1 in ri
    if m2 == 2
        apply!(network[control_register,control_qubit], Z)
        @yield timeout(sim, single_qubit_gate_exec_time)
    end
    unlock(network[control_register,control_qubit])
    unlock(network[target_register,target_qubit])
    unlock(network[control_register,1])
    unlock(network[target_register,1])
end

@resumable function steane_encoding_circuit(sim, network)
    # Initialize memory qubits in zero state
    initialize!([network[i,j] for i in 1:length(sizes) for j in 2:sizes[i]], SProjector(Z1⊗Z1⊗Z1⊗Z1⊗Z1⊗Z1⊗Z1)  ) 

    # TODO: For each layer, identify the longest execution time among parallel operations
    # make adding gates possible
    @yield @process single_qubit_gate(sim, network, 1, 2, H)
    @yield @process single_qubit_gate(sim, network, 1, 3, H)
    @yield @process single_qubit_gate(sim, network, 1, 5, H)

    @yield @process CNOT_gate(sim, network, 1, 2, 4)
    @yield @process telegate_CNOT(sim, network, 1, 2, 5, 2)
    @yield @process telegate_CNOT(sim, network, 1, 2, 3, 3)
    @yield @process CNOT_gate(sim, network, 1, 2, 6)
    @yield @process telegate_CNOT(sim, network, 1, 2, 5, 3)
    @yield @process telegate_CNOT(sim, network, 1, 2, 3, 2)
    @yield @process telegate_CNOT(sim, network, 1, 2, 2, 3)
    @yield @process CNOT_gate(sim, network, 1, 3, 4)
    @yield @process CNOT_gate(sim, network, 1, 5, 6)
end


# Initialise the network

R = length(sizes) # Number of registers
# All of the quantum register we will be simulating
registers = Register[]
for s in sizes
    traits = [Qubit() for _ in 1:s]
    repr = [QuantumOpticsRepr() for _ in 1:s]  # recast to Clifford
    bg = [T2Dephasing(T2) for _ in 1:s]          # define other noise on the registers?
    push!(registers, Register(traits,repr, bg)) 
end


graph = grid([R])
network = RegisterNet(graph, registers) # A graphs with extra "meta data"

# The scheduler datastructure for the discrete event simulation
sim = get_time_tracker(network)

# Add a register datastructures and event locks to each node.
#for v in vertices(network)
#    # Create an array specifying whether a qubit is entangled with another qubit
#    network[v,:enttrackers] = Any[nothing for i in 1:sizes[v]]
#end

#sim, network = simulation_setup(sizes, T2; representation = QuantumOpticsRepr)  
@process steane_encoding_circuit(sim, network)

run(sim, t)

#step_ts = range(0, 50, step=0.2)

#ts = Observable(Float64[0])
#fidelities = Observable(Float64[0])

steane_7_state = StabilizerState("ZIZIZIZ XIXIXIX IZZIIZZ IXXIIXX IIIZZZZ IIIXXXX ZZZZZZZ") 

fidelity = real(observable(vcat( [network[1,i] for i in 2:6],[network[2,j] for j in 2:3]), SProjector(steane_7_state)))
print("Fidelity: $fidelity at time $t \n")

end


