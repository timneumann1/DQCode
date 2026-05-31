# types.jl 

"""
Type definitions used in DQCode for optimisation and DQC simulation of fault-tolerant encoding circuits.
"""
module Types

using QuantumClifford
using QECCore: AbstractCSSCode
export CodeParameters, GeneticParameters, NetworkSpecifications, MCTSParameters, NoiseSpecs


"""
    CodeParameters

Type that collects relevant code parameters for a QEC code. Used in the optimisation of 
logical state encoding as well as in the DQC simulation.

### Fields

- `qec_code` -- QEC code for which the logical zero state is to be optimised
- `num_X_checks` -- number of X checks in the canonical tableau
- `logical_Zs` -- canonical logical Z operators
- `target_state` -- logical zero target state in canonical form
- `target_bit_matrix` -- target state tableau converted to a bit matrix
- `n` -- number of qubits
- `k` -- dimensions of logical space
- `distance` -- QEC code distance
"""
struct CodeParameters
    qec_code::AbstractCSSCode
    num_X_checks::Int
    logical_Zs::QuantumClifford.AbstractStabilizer
    target_state::QuantumClifford.Stabilizer{QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}}}
    target_bit_matrix::Matrix{Int}
    n::Int
    k::Int
    distance::Int
end


"""
    NetworkSpecifications

Type that collects relevant networking specifications for the Type-II architecture.
    
### Fields

- `register_sizes` -- sizes (number of memory qubits) of each QPU register
- `num_registers` -- total number of registers in the architecture
- `mapping` -- initial placement of data qubits extracted from hypergraph partitioning.
- `mapping_transpositions` -- sequence of transpositions corresponding to the initial mapping (applied right-to-left)
- `inv_map` -- inverse mapping used during circuit execution
- `register_lookup_array` -- lookup array indicating which register core each register slot index belongs to
- `data_qubits` -- list of indices corresponding to the data qubits
- `comm_qubits` -- list of indices corresponding to the communication qubits
- `num_data_qubits` -- total number of data qubits (needs to be equal to `n` in CodeParameters)
- `num_comm_qubits` -- total number of communication qubits across all registers
- `num_comm_qubits_per_register` -- number of communication qubits allocated to each register
- `num_data_and_comm_qubits` -- total combined number of data + communication qubits in the distributed network
"""
struct NetworkSpecifications
    register_sizes::Vector{Int}        
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
end


"""
    GeneticParameters

Type that collects optimisation parameters for the genetic search.

### Fields
- `num_individuals` -- size of the population per generation
- `num_generations` -- total number of generations in the genetic algorithm
- `max_len` -- maximum allowed number of gates in each circuit individual
- `mutation_rate` -- probability of an offspring individual undergoing mutation after crossover
- `tournament_size` -- number of candidate individuals participating in tournament selection
- `selection_ratio` -- fraction of the population selected as parents for the next generation
- `num_elite` -- number of best-performing individuals passed directly into the next generation unmodified
- `fitness_weights` -- vector of weight coefficients for the multi-objective fitness function, weighing single-, two-qubit and telegates
- `tableau_metric` -- string identifier for the metric used to evaluate proximity of the current tableau to the target state
"""
struct GeneticParameters
    num_individuals::Int
    num_generations::Int
    max_len::Int
    mutation_rate::Float64
    tournament_size::Int
    selection_ratio::Float64
    num_elite::Int
    fitness_weights::Vector{Float64}
    tableau_metric::String
end


"""
    MCTSParameters

Type that collects optimisation parameters for the Monte Carlo Tree Search (MCTS) algorithm.

### Fields

- `max_steps` -- maximum number of steps (actions) that the solver takes, capping the maximum circuit size
- `fitness_weights` -- vector of weight coefficients for the multi-objective reward/cost function
- `discount_factor` -- decay rate applied to future rewards during backpropagation 
- `tableau_metric` -- string identifier for the metric used to evaluate proximity of the current tableau to the target state
- `reuse_tree` -- whether the solver retains and reuses the generated search tree across sequential steps
- `depth` -- maximum rollout depth the solver traverses when evaluating a node
- `n_iterations` -- number of tree search iterations the solver performs in expansion to choose the best next action
- `exploration_constant` -- constant balancing exploration vs. exploitation in expansion
"""
struct MCTSParameters
    max_steps::Int 
    fitness_weights::Vector{Float64}
    discount_factor::Float64 
    tableau_metric::String
    reuse_tree::Bool
    depth::Int 
    n_iterations::Int
    exploration_constant::Float64
end

"""
    NoiseSpecs

Type that collects circuit-level and Bell pair initialisation noise definitions for DQC simulations.

### Fields

- `n_samples` -- total number of Monte Carlo trajectories to simulate during the simulation
- `p` -- base probability for memory/communication qubit initialisation, measurement, and two-qubit errors
- `p_idle` -- probability of an idle memory error occurring over a standard depth-1 layer
- `p_idle_telegate_layer` -- idle error probability for qubits waiting during a telegate interaction (accounts for longer/probabilistic creation time)
- `p_single` -- error probability for local single-qubit gates
- `p_bell` -- error probability or infidelity of Bell-pair generation/initialisation for photonic interconnect
"""
struct NoiseSpecs 
    n_samples::Int64
    p::Float64 
    p_idle::Float64
    p_idle_telegate_layer::Float64 
    p_single::Float64
    p_bell::Float64 
end


end