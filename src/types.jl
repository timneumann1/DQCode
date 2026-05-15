# types.jl

module Types

using QuantumClifford

# types for simulation
using QuantumClifford: AbstractOperation#, AbstractCliffordOperator
using QECCore: AbstractCSSCode

import Quantikz: QuantikzOp, ClassicalDecision

export Circuit, CircuitIndividual #SimulationParameters, SimulationFidelity,
#export GateSet
#export Gate, GateSet, SingleQubitGate, TwoQubitGate
export HadamardGate, IdentityGate, PauliXGate, PauliYGate, PauliZGate, CX_Gate, CZ_Gate, SGate, InvSGate, SqrtXGate, InvSqrtXGate, SWAP_Gate, ConditionalGate
export CodeParameters, GeneticParameters, NetworkSpecifications, MCTSParameters, NoiseSpecs

# Parameter structures

# struct SimulationParameters
#     register_sizes::Vector{Int}        # Number of qubits in each register
#     #T1_relaxation::Float64          # T1 relaxation time of all qubits   [Wang]
#     depolarising_noise::Float64       #  T2_dephasing::Float64         # T2 dephasing time of all qubits     [Wang]
#     single_qubit_gate_exec_time::Float64   # Execution Time of a single qubit gate
#     two_qubit_gate_exec_time::Float64      # Execution Time of a two qubit gate
#     projective_measurement_time::Float64      # Time to peform a measurement
#     classical_communication_time::Float64   # define in relation to speed of photons in fiber cable
#     bell_state_fidelity::Float64
#     success_prob::Float64
#     attempt_time::Float64
#     #TODO: Add measurement fidelity
# end

# struct SimulationFidelity
#     fidelity::Float64
# end


# Gate Types

# struct GateSet
#     single_qubit_gates::Vector{DataType}
#     two_qubit_gates::Vector{DataType}
# end

# abstract type Gate end
# abstract type SingleQubitGate <: Gate end
# abstract type TwoQubitGate <: Gate end

# # single-qubit gates
# struct IdentityGate <: SingleQubitGate 
#     index::Int
# end

# struct PauliXGate <: SingleQubitGate
#     index::Int
# end

# struct PauliYGate <: SingleQubitGate
#     index::Int
# end

# struct PauliZGate <: SingleQubitGate 
#     index::Int
# end

# struct HadamardGate <: SingleQubitGate 
#     index::Int
# end

# struct SGate <: SingleQubitGate 
#     index::Int
# end

# struct InvSGate <: SingleQubitGate 
#     index::Int
# end

# struct SqrtXGate <: SingleQubitGate 
#     index::Int
# end

# struct InvSqrtXGate <: SingleQubitGate 
#     index::Int
# end

# # two-qubit gates (store the counterpart of the operation)
# struct CX_Gate <: TwoQubitGate
#     control::Int
#     target::Int
# end

# struct CZ_Gate <: TwoQubitGate
#     control::Int
#     target::Int
# end

# struct SWAP_Gate <: TwoQubitGate
#     qubit_1::Int
#     qubit_2::Int
# end

struct ConditionalGate <: AbstractOperation
    truegate::AbstractOperation
    falsegate::AbstractOperation
    controlbit::Int
end

# Genetic Algorithm Types

# mutable struct CircuitIndividual
#     #gates::Vector{Gate}   
#     gates::Vector{AbstractOperation}
# end

# function CircuitIndividual(num_gates::Int)
#     @assert num_gates > 0
#     #gates = Gate[IdentityGate(1) for _ in 1:num_gates]
#     gates = [sId1(1) for _ in 1:num_gates]
#     return CircuitIndividual(gates)
# end

# mutable struct Circuit
#     gates::Matrix{Gate}
# end

# function Circuit(num_qubits::Int, num_layers::Int)
#     @assert num_qubits > 0
#     @assert num_layers > 0
#     Circuit(fill(IdentityGate(1), num_qubits, num_layers))
# end


struct CodeParameters
    qec_code::AbstractCSSCode
    #stabilizers::Stabilizer{QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}}}
    num_X_checks::Int
    logical_Zs
    target_state::Stabilizer{QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}}}
    target_bit_matrix::Matrix{Int}
    n::Int
    k::Int
    distance::Int
end

struct NetworkSpecifications
    register_sizes::Vector{Int}        # Number of qubits in each register
    num_registers::Int
    mapping::Vector{Int}
    mapping_transpositions::Vector{Tuple{Int, Int}}
    inv_map::Vector{Int}
    register_lookup_array::Vector{Int}
    data_qubits::Vector{Int}
    comm_qubits::Vector{Int}
    num_data_qubits::Int
    num_comm_qubits::Int
    num_comm_qubits_per_register::Int
    num_data_and_comm_qubits::Int
    #comm_idx::Vector{Int}
    #comm_inv_perm_idx::Vector{Int}
end


# struct OptimisationParameters
#     tableau_metric::String
#     gate_set::GateSet
# end

struct GeneticParameters
    #gate_set::GateSet
    num_individuals::Int
    num_generations::Int
    max_len::Int
    mutation_rate::Float64
    tournament_size::Int
    selection_ratio::Float64
    num_elite::Int
    #standard_encoding::Bool
    #warm_start::Bool
    fitness_weights::Vector{Float64}
    #qec_code::AbstractCSSCode
    tableau_metric::String
end

struct MCTSParameters
    #gate_set::GateSet
    max_steps::Int # Maximum number of steps (actions) that the solver takes, is equivalent to maximum circuit size
    fitness_weights::Vector{Float64}
    discount_factor::Float64 
    tableau_metric::String
    reuse_tree::Bool
    depth::Int # depth that the solver traverses to maximally in each rollout
    n_iterations::Int # number of iterations the solver rolls out to choose the best next action
    exploration_constant::Float64
end


 struct NoiseSpecs # Defining circuit-level noise
    n_samples::Int64
    p::Float64 # Captures memory initialisation, memory/comm single- and two-qubit gates, depth-1 decoherence and measurement noise
    #p_idle_telegate_layer::Float64 # Idle error probability for layer that includes a telegate (accounting for longer/probabilistic creation time)
    p_mixed::Float64 # Captures two-qubit gate noise between species
    p_bell::Float64 # Captures Bell states initialisation of communication qubits via photonic interconnects
 end


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