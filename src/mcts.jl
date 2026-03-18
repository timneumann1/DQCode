module MonteCarloTreeSearch

using ..Types
using ..CircuitSimulator
using ..Helper
using ..LogicalEnc
using ..TrivariateBicycleCode
using ..Genetic: data_qubit_partitioning, tableau_to_bitmatrix, tableau_distance

using POMDPs, POMDPTools
using MCTS

using Random
#rng = MersenneTwister(42)
# using Quantikz: savecircuit, @with, classicalbitslayout
using QECCore
using QECCore: Steane7, QuantumTannerGraphProduct, CyclicQuantumTannerGraphProduct, Triangular488
# using QuantumClifford.ECC: DistanceMIPAlgorithm
# using HiGHS
using QuantumClifford
using QuantumClifford: MixedDestabilizer, sHadamard, sCNOT, sSWAP, @S_str, true_success_stat, false_success_stat, continue_stat, failure_stat, PauliMeasurement, VerifyOp
# using BenchmarkTools
# using CairoMakie
# using KaHyPar
# using SparseArrays
using D3Trees

export run_MCTS

function define_parameters()

    #TODO: Rename or introduce a SimulationParameters type for this sim as well (in addition to DTS)
    networking_params = NetworkingParameters(
        [3,3,3,3], #register sizes of Type-II architecture (here: only fewn memory qubits per core), CircuitSim automatically adds comm. qubits (ancillas are only added in DTS)
        0.0, # depolarising_prob 
        0.0, # gate_noise_prob 
        0.0, # Telegate noise (depolarising channel)
    )

    mcts_params = MCTSParameters(
        250, # steps before termination
        30, # max length
        1, #shots
        TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]),# qec code
        "hamming", # tableau distance metric
        2500, #iterations
        1.5 # exploration constant
        )
    return networking_params, mcts_params
end


struct EncodingState
    gates::Vector{Gate}   # raw Gate[] as in CircuitIndividual, no comm qubit indexing
end

function Base.hash(s::EncodingState, h::UInt)
    hash(s.gates, h)
end
function Base.:(==)(a::EncodingState, b::EncodingState)
    length(a.gates) == length(b.gates) && all(a.gates .== b.gates)
end

Base.hash(g::HadamardGate, h::UInt) = hash((:H, g.index), h)
Base.hash(g::SGate, h::UInt) = hash((:S, g.index), h)

Base.hash(g::CNOT_Gate, h::UInt)    = hash((:CNOT, g.control, g.target), h)
Base.:(==)(a::CNOT_Gate, b::CNOT_Gate) = a.control == b.control && a.target == b.target


struct EncodingMDP <: MDP{EncodingState, Gate}
    # From your run_genetic_search setup:
    num_data_qubits::Int
    max_length::Int
    target_bit_matrix::Matrix{Int}        # from tableau_to_bitmatrix
    tableau_metric::String                 # "hamming" or "jaccard"

    # DQC architecture context (passed through to construct_executable_circuit):
    networking_params::NetworkingParameters
    mapping::Vector{Tuple{Int,Int}}
    inv_perm::Vector{Int}
    register_lookup_array::Vector{Int}
    data_qubits::Vector{Int}
    comm_qubits::Vector{Int}
    num_comm_qubits_per_register::Int
    num_qubits::Int
    data_qubit_capacities::Vector{Int}
    num_registers::Int
    
    discount_factor::Float64
end

function POMDPs.actions(mdp::EncodingMDP)
    n = mdp.num_data_qubits
    actions = Gate[]
    for i in 1:n
        push!(actions, HadamardGate(i))
        push!(actions, SGate(i))
    end
    for c in 1:n, t in 1:n
        c == t && continue
        push!(actions, CNOT_Gate(c, t))
    end
    return actions
end

function POMDPs.gen(mdp::EncodingMDP, s::EncodingState, a::Gate, rng)
    # Append the new gate to the current circuit
    new_gates = vcat(s.gates, [a])
    sp = EncodingState(new_gates)

    # Build the executable DQC circuit (handles telegates, comm qubits, etc.)
    exec_circuit = construct_executable_circuit(
        mdp.networking_params.depolarising_noise,
        mdp.networking_params.gate_noise,
        mdp.networking_params.telegate_noise,
        new_gates,
        mdp.mapping, mdp.inv_perm, mdp.register_lookup_array,
        mdp.data_qubits, mdp.num_comm_qubits_per_register,
        mdp.num_qubits, mdp.data_qubit_capacities
    )

    # Single-shot tableau evaluation (num_shots=1 for MCTS speed; see note below)
    mc_result = execute_circuit(exec_circuit, mdp.num_qubits, mdp.num_registers; num_traj=1)
    
    # Compute reward from tableau distance — mirrors evaluate_population exactly
    stab = mc_result[1]
    stab_view = stabilizerview(stab)
    stab_view = traceout!(copy(stab_view), mdp.comm_qubits)
    stab_canon = canonicalize_rref!(stab_view)
    tableau = tab(stab_canon[1])
    current_bit_matrix = tableau_to_bitmatrix(tableau)
    
    dist = tableau_distance(current_bit_matrix, mdp.target_bit_matrix,
                            mdp.data_qubits, mdp.comm_qubits, mdp.tableau_metric)
    
    fidelity = 1.0-dist
    # Potential-based shaping: γΦ(s') - Φ(s), where Φ = -distance
    # This preserves the optimal policy while providing dense signal
    r = 1e6*fidelity+ 1e-6*rand(rng)-length(exec_circuit)  # via the discount factor, large depth will be penalised
    
    return (sp=sp, r=r, fidelity=1.0-dist, gate_count=length(exec_circuit) )
end

POMDPs.discount(mdp::EncodingMDP) = mdp.discount_factor

function POMDPs.isterminal(mdp::EncodingMDP, s::EncodingState)
    # Terminate at max depth or if already perfectly encoding
    if length(s.gates) >= mdp.max_length
        return true
    end
    # Optional: early termination on perfect fidelity (distance == 0)
    # This requires re-evaluating the tableau, so only worth doing if 
    # you cache the last evaluated distance in state (see extension below)
    return false
end

function POMDPs.initialstate(mdp::EncodingMDP)
    # Warm start: seed with first half of compiler circuit, matching your GA
    # Or start from empty circuit:
    return Deterministic(EncodingState(Gate[]))
end

function build_encoding_mdp(networking_params, mcts_params)
    # Replicate the setup from run_genetic_search
    permutation = data_qubit_partitioning(networking_params, Stabilizer(mcts_params.qec_code))
    inv_perm = invperm(permutation)
    mapping = perm_to_transpositions(deepcopy(permutation))
    
    qec_code_required_qubits = code_n(mcts_params.qec_code)
    data_qubit_capacities = networking_params.register_sizes
    num_registers = length(data_qubit_capacities)
    register_lookup_array, data_qubits, num_data_qubits, num_comm_qubits_per_register = 
        create_lookup_array_cliff(data_qubit_capacities)
    num_qubits = num_data_qubits + num_comm_qubits_per_register * num_registers
    all_qubits = collect(1:num_qubits)
    comm_qubits = setdiff(all_qubits, data_qubits)

    code = MixedDestabilizer(mcts_params.qec_code)
    target_state = vcat(stabilizerview(code), logicalzview(code))
    target_canon = canonicalize_rref!(target_state)
    target_tableau = tab(target_canon[1])
    target_bit_matrix = tableau_to_bitmatrix(target_tableau)

    return EncodingMDP(
        num_data_qubits, mcts_params.max_length,
        target_bit_matrix, mcts_params.tableau_metric,
        networking_params, mapping, inv_perm,
        register_lookup_array, data_qubits, comm_qubits,
        num_comm_qubits_per_register, num_qubits,
        data_qubit_capacities, num_registers,
        0.999
    )
end


function make_mcts_planner(networking_params, mcts_params)
    mdp = build_encoding_mdp(networking_params, mcts_params)
    rng = MersenneTwister(rand(UInt))

    solver = MCTSSolver(
        n_iterations = mcts_params.n_iterations,
        depth = mcts_params.max_length,
        exploration_constant = mcts_params.exploration_constant,
        rng = rng,
        reuse_tree = false,
        enable_tree_vis = true,
        estimate_value = 0.0,
    )
    planner = solve(solver, mdp)
    return mdp, planner
end


# To extract the full optimised circuit:
function run_mcts_search(mdp, planner, max_steps)
    s = EncodingState(Gate[])
    final_gate_count = typemax(Int)
    for _ in 1:max_steps
        POMDPs.isterminal(mdp, s) && break
        a = action(planner, s)
        s, r, fidelity, gate_count = POMDPs.gen(mdp, s, a, Random.GLOBAL_RNG)
        println("Applied $(typeof(a)), reward=$r for fidelity=$fidelity and DQC_depth=$gate_count, MC tree depth=$(length(s.gates))")
        if fidelity >= 1.0
            println("Stopping early: fidelity reached $(fidelity).")
            final_gate_count = gate_count
            break
        end
    end
    return (final_gates=s.gates, DQC_gate_count = final_gate_count)
end

function run_MCTS()
    networking_params, mcts_params = define_parameters()
    mdp, planner = make_mcts_planner(networking_params, mcts_params)
    result = run_mcts_search(mdp, planner, mcts_params.max_steps)
    best_circuit = result.final_gates
    #a, info = action_info(planner, EncodingState(Gate[]))
    #inchrome(D3Tree(info[:tree]))
    print("\n\n\nBest circuit: $best_circuit, \nconsisting of $(result.DQC_gate_count) DQC gates.\n\n\n")
    return best_circuit, mdp, planner
end

#mdp = MyMDP() # initializes the MDP
#solver = MCTSSolver(n_iterations=50, depth=20, exploration_constant=5.0) # initializes the Solver type
#planner = solve(solver, mdp) # initializes the planner

end