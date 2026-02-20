# types.jl

module Types

using QuantumClifford

# types for simulation
using QuantumClifford: AbstractOperation

import Quantikz: QuantikzOp, ClassicalDecision

export SimulationParameters, SimulationFidelity, Circuit, HadamardGate, IdentityGate, PauliXGate, PauliYGate, PauliZGate, CNOT_Gate, SWAP_Gate, Gate, ConditionalGate, GeneticParameters, NetworkingParameters


struct ConditionalGate <: AbstractOperation
    truegate::AbstractOperation
    falsegate::AbstractOperation
    controlbit::Int
end

# For Quantikz
function _conditional_gate_label(g::AbstractOperation)
    repr_g = string(g)
    if occursin("sX", repr_g)
        return "X"
    elseif occursin("sY", repr_g)
        return "Y"
    elseif occursin("sZ", repr_g)
        return "Z"
    elseif occursin("sHadamard", repr_g)
        return "H"
    elseif occursin("sId1", repr_g)
        return "I"
    end
    return "U"
end

function QuantikzOp(op::ConditionalGate)
    targets = collect(affectedqubits(op.truegate))
    label = _conditional_gate_label(op.truegate)
    return ClassicalDecision(label, targets, op.controlbit)
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

struct NetworkingParameters
    register_sizes::Vector{Int}        # Number of qubits in each register
    depolarising_noise::Float64        # Circuit Noise probability
    gate_noise::Float64                 # Gate Noise probability
    telegate_noise::Float64
end

struct GeneticParameters
    num_individuals::Int
    num_generations::Int
    num_shots::Int
    mutation_rate::Float64
    tournament_size::Int
    selection_ratio::Float64
    depth::Int
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

mutable struct Circuit
    gates::Matrix{Gate}
end

function Circuit(num_qubits::Int, num_layers::Int)
    @assert num_qubits > 0
    @assert num_layers > 0
    Circuit(fill(IdentityGate(), num_qubits, num_layers))
end

end