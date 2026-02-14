# types.jl

module Types

using QuantumClifford

# types for simulation
using QuantumClifford: AbstractOperation

export SimulationParameters, SimulationFidelity, Circuit, HadamardGate, IdentityGate, PauliXGate, PauliYGate, PauliZGate, CNOT_Gate, SWAP_Gate, Gate, ConditionalGate


struct ConditionalGate <: AbstractOperation
    truegate::AbstractOperation
    falsegate::AbstractOperation
    controlbit::Int
end

struct SimulationParameters
    register_sizes::Vector{Int}        # Number of qubits in each register
    #T1_relaxation::Float64          # T1 relaxation time of all qubits   [Wang]
    depolarising_noise::Float64       #  T2_dephasing::Float64         # T2 dephasing time of all qubits     [Wang]
    single_qubit_gate_exec_time::Float64   # Execution Time of a single qubit gate
    two_qubit_gate_exec_time::Float64      # Execution Time of a two qubit gate
    projective_measurement_time::Float64      # Time to peform a measurement
    classical_communication_time::Float64   # define in relation to speed of photons in fiber cable
    bell_state_fidelity::Float64
    success_prob::Float64
    attempt_time::Float64
    #TODO: Add measurement fidelity
end

struct SimulationFidelity
    fidelity::Float64
end

### types for circuit representation

abstract type Gate end

# single-qubit gates
struct IdentityGate <: Gate end
struct PauliXGate <: Gate end
struct PauliYGate <: Gate end
struct PauliZGate <: Gate end
struct HadamardGate <: Gate end

# two-qubit gates (store the counterpart of the operation)
struct CNOT_Gate <: Gate
    control::Int
    target::Int
end

struct SWAP_Gate <: Gate
    qubit_1::Int
    qubit_2::Int
end

struct Circuit
    gates::Matrix{Gate}
end

function Circuit(num_qubits::Int, num_layers::Int)
    @assert num_qubits > 0
    @assert num_layers > 0
    Circuit(fill(IdentityGate(), num_qubits, num_layers))
end

end