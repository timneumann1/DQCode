# types.jl

module Types

using QuantumClifford

# types for simulation
using QuantumClifford: AbstractOperation
using QECCore: AbstractCSSCode

import Quantikz: QuantikzOp, ClassicalDecision

export SimulationParameters, SimulationFidelity, Circuit, CircuitIndividual, HadamardGate, IdentityGate, PauliXGate, PauliYGate, PauliZGate, CNOT_Gate, SGate, SWAP_Gate, Gate, ConditionalGate
export CodeParameters, GeneticParameters, NetworkSpecifications, MCTSParameters, OptimisationParameters

# Parameter structures

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

struct CodeParameters
    qec_code::AbstractCSSCode
    stabilizers::Stabilizer{QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}}}
    target_state::Stabilizer{QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}}}
    target_bit_matrix::Matrix{Int}
    distance::Int
end

struct OptimisationParameters
    tableau_metric::String
end

struct NetworkSpecifications
    register_sizes::Vector{Int}        # Number of qubits in each register
    num_registers::Int
    permutation::Vector{Int}
    mapping::Vector{Tuple{Int, Int}}
    inv_perm::Vector{Int}
    register_lookup_array::Vector{Int}
    data_qubits::Vector{Int}
    comm_qubits::Vector{Int}
    num_data_qubits::Int
    num_comm_qubits_per_register::Int
    num_qubits::Int
    comm_idx::Vector{Int}
    comm_inv_perm_idx::Vector{Int}
    depolarising_noise::Float64        # Circuit Noise probability
    gate_noise::Float64                 # Gate Noise probability
    telegate_noise::Float64
    num_shots::Int
end

struct GeneticParameters
    num_individuals::Int
    num_generations::Int
    max_len::Int
    mutation_rate::Float64
    tournament_size::Int
    selection_ratio::Float64
    num_elite::Int
    warm_start::Bool
    #qec_code::AbstractCSSCode
    #tableau_metric::String
end

struct MCTSParameters
    max_steps::Int
    n_iterations::Int
    exploration_constant::Float64
end

struct SimulationFidelity
    fidelity::Float64
end


# Gate Types

abstract type Gate end

# single-qubit gates
struct IdentityGate <: Gate 
    index::Int
end

struct PauliXGate <: Gate
    index::Int
end

struct PauliYGate <: Gate
    index::Int
end

struct PauliZGate <: Gate 
    index::Int
end

struct HadamardGate <: Gate 
    index::Int
end

struct SGate <: Gate 
    index::Int
end


# two-qubit gates (store the counterpart of the operation)
struct CNOT_Gate <: Gate
    control::Int
    target::Int
end

struct SWAP_Gate <: Gate
    qubit_1::Int
    qubit_2::Int
end

struct ConditionalGate <: AbstractOperation
    truegate::AbstractOperation
    falsegate::AbstractOperation
    controlbit::Int
end

# Genetic Algorithm Types

mutable struct CircuitIndividual
    gates::Vector{Gate}   
end

function CircuitIndividual(num_gates::Int)
    @assert num_gates > 0
    gates = Gate[IdentityGate(1) for _ in 1:num_gates]
    return CircuitIndividual(gates)
end

# mutable struct Circuit
#     gates::Matrix{Gate}
# end

# function Circuit(num_qubits::Int, num_layers::Int)
#     @assert num_qubits > 0
#     @assert num_layers > 0
#     Circuit(fill(IdentityGate(1), num_qubits, num_layers))
# end

# Plotting 

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

end