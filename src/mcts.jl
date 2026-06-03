# mcts.jl

"""
Monte Carlo Tree search for efficient encoding circuits of CSS QEC codes.
"""
module MonteCarloTreeSearch

export monte_carlo_tree_search

using ..Types
using ..Helper

using QuantumClifford: MixedDestabilizer, Tableau, AbstractOperation, stabilizerview, AbstractSingleQubitOperator,
                        AbstractTwoQubitOperator, sHadamard, sCNOT, tab, canonicalize_rref!, Register, mctrajectory!
using POMDPs, POMDPTools
using MCTS
using Random
using Serialization, CSV, DataFrames
using ProgressMeter
                        

"""
    monte_carlo_tree_search(code_params::CodeParameters, 
                            network_specs::NetworkSpecifications, 
                            mcts_params::MCTSParameters)::Tuple{Vector{AbstractOperation}, Bool, Vector{Int}, 
                                                                Vector{Float64}, Vector{Vector{Int}}, Vector{Float64}}

Execute a Monte Carlo Tree Search (MCTS) to find efficient logical zero state encoding circuits for CSS codes.

### Input

- `code_params` -- parameters defining the target QEC code and its logical zero state
- `network_specs` -- hardware networking specifications for the underlying DQC architecture
- `mcts_params` -- Monte Carlo Tree Search hyper- and configuration parameters

### Output

Returns a 6-element Tuple containing the optimised MCTS circuit, its verification and gate count as well as the
vectors capturing the evoluation of fidelity, gate counts and reward throughout the search. `POMDPs.gen`

### Notes

We assume a H-CNOT template, where we initially populate a circuit with `num_X_checks` Hadamards. Large circuits 
are penalised via `fitness_weights` and `discount_factor`.
"""
function monte_carlo_tree_search(code_params::CodeParameters, 
                                network_specs::NetworkSpecifications, 
                                mcts_params::MCTSParameters)::Tuple{Vector{AbstractOperation}, Bool, Vector{Int},
                                                                    Vector{Float64}, Vector{Vector{Int}}, Vector{Float64}}
    @info "Initialised Monte Carlo Tree Search with depth $(mcts_params.depth) and $(mcts_params.n_iterations) iterations ..."
    mdp = EncodingMDP(code_params, network_specs, mcts_params)
    solver = MCTSSolver(
        n_iterations = mdp.mcts_params.n_iterations, # if no solution is found, this `n_iterations` be increased
        depth = mcts_params.depth, # "
        exploration_constant = mcts_params.exploration_constant, # "
        rng = Random.GLOBAL_RNG,
        reuse_tree = mcts_params.reuse_tree,
        enable_tree_vis = false,
        estimate_value = 0.0
    )
    policy = solve(solver, mdp)
    s = initialstate(mdp)
    MCTS_gate_counts = (typemax(Int), typemax(Int), typemax(Int))
    MCTS_fidelity = typemin(Float16)
    MCTS_circuit_state = s
    gate_count_evolution = Vector{Vector{Int64}}()
    fidelity_evolution = Float64[]
    reward_evolution = Float64[]
    p = Progress(mdp.mcts_params.max_steps; desc = "Monte Carlo Tree search", showspeed = true)
    #------------------- MCTS Steps --------------------------------    
    for step in 1:mdp.mcts_params.max_steps 
        a = action(policy, s)    
        s, r = POMDPs.gen(mdp, s, a, Random.GLOBAL_RNG)
        MCTS_circuit_state = s
        MCTS_gate_counts = s.gate_counts
        MCTS_fidelity = s.fidelity
        next!(p; showvalues = [
            (:step, step),
            (:fidelity, round(MCTS_fidelity, digits=5)),
            (:reward, round(r, digits=4)),
            (:gates, MCTS_gate_counts)
        ])
        push!(gate_count_evolution, MCTS_gate_counts)
        push!(fidelity_evolution, MCTS_fidelity)
        push!(reward_evolution, r)
        if POMDPs.isterminal(mdp, s)#, fidelity=fidelity) 
            @info "Terminal condition reached after $step steps. Final fidelity: $(s.fidelity)"
            finish!(p)
            break
        end
    end
    MCTS_circuit = MCTS_circuit_state.circuit    
    verification_logical_state = verify_success(copy(MCTS_circuit), mdp.code_params.target_state, mdp.network_specs)
    @info "Verification of MCTS circuit successful: $verification_logical_state"
    return MCTS_circuit, verification_logical_state, MCTS_gate_counts, 
            fidelity_evolution, gate_count_evolution, reward_evolution
end


"""
    POMDPs.initialstate(mdp::EncodingMDP)::CircuitState

Initialise the tree with an empty circuit state.

### Input

- `mdp` -- Markov Decision Process environment

### Output

Returns a `CircuitState` instantiated with an empty circuit array, a `MixedDestabilizer` 
tableau respresenting the all-zero state, and a fidelity placeholder of `0.0`. 
"""
function POMDPs.initialstate(mdp::EncodingMDP)::CircuitState
    init_state = one(MixedDestabilizer,mdp.network_specs.num_data_qubits)
    init_bit_matrix = tableau_to_bitmatrix( tab(canonicalize_rref!( stabilizerview(init_state) )[1])  ) 
    return CircuitState([], init_state, init_bit_matrix, [0,0,0], 0.0)
end

"""
    POMDPs.discount(mdp::EncodingMDP)::Float64

Expose the discount factor, determining the weight of future reward in the MCTS.
"""
POMDPs.discount(mdp::EncodingMDP) = mdp.mcts_params.discount_factor


"""
    POMDPs.isterminal(mdp::EncodingMDP, s::CircuitState)::Bool

Evaluate the termination condition for a specific node in the MCTS tree

### Input

- `mdp` -- The encoding Markov Decision Process environment (unused but required in POMPDs framework)
- `s` -- Current circuit state to evaluate

### Output

Only returns `true` if the state's simulation fidelity evaluates to `1.0`, i.e., if the circuit
encodes the desired logical zero state, leading to early termination of the MCTS.
"""
function POMDPs.isterminal(mdp::EncodingMDP, s::CircuitState)::Bool
    if s.fidelity >= 1.0 
        return true
    end
    return false
end


"""
    POMDPs.actions(mdp::EncodingMDP, s::CircuitState)::Vector{AbstractOperation}

Define the action space available from the current MCTS circuit state node. 

### Input

- `mdp` -- The Markov Decision Process environment
- `s` -- Current circuit state to take next action from

### Output

Returns an array of candidate `AbstractOperation` implementations that can be appended to the current circuit. 

### Notes

In the H-CNOT template, we restrict to (i) unused qubits in the Hadamard layer, 
and (ii) active qubits in the CNOT layer.
"""
function POMDPs.actions(mdp::EncodingMDP, s::CircuitState)::Vector{AbstractOperation}
    n = mdp.network_specs.num_data_qubits
    actions = Vector{AbstractOperation}()
    affected_qubits = [(g isa AbstractTwoQubitOperator ? g.q2 : g.q) for g in s.circuit]
    if length(s.circuit) < mdp.code_params.num_X_checks # Hadamards == # of X stabilisers 
        for i in 1:n
            if length(s.circuit)==0
                push!(actions, sHadamard(i))
            elseif i ∉ affected_qubits
                push!(actions, sHadamard(i))
            end
        end
    else 
        for c in 1:n, t in 1:n
            c == t && continue
            if last(s.circuit) isa AbstractSingleQubitOperator
                if c ∈ affected_qubits
                    push!(actions, sCNOT(c, t))
                end
            elseif !( (last(s.circuit).q1 == c && last(s.circuit).q2 ==t) || (c ∉ affected_qubits)) 
                push!(actions, sCNOT(c, t))
            end
        end
    end
    return actions
end


"""
    POMDPs.gen(mdp::EncodingMDP, state::CircuitState, action::AbstractOperation, rng)::NamedTuple{(:sp, :r), Tuple{CircuitState, Float64}}

Generative interface for Monte Carlo Tress search, defining the generation and evaluation of a new 
action as well as the new state including its reward.

### Input

- `mdp` -- The encoding Markov Decision Process environment
- `state` -- The current `CircuitState` state
- `action` -- The `AbstractOperation` action selected to be appended to the circuit state
- `rng` -- The random number generator instance (unused but required in POMPDs framework)

### Output

Returns a named tuple `(sp, r)` defining the new state (`sp` stands for `s prime`) and its reward (`r`).

### Notes

The transition from the state-action pair `(s,a)` to the new state `sp`` is deterministic. The reward is updated based
on the nature of the appended action and the marginal fidelity gain/loss.
"""
function POMDPs.gen(mdp::EncodingMDP, state::CircuitState, action::AbstractOperation, rng::Random.TaskLocalRNG)::NamedTuple{(:sp, :r), Tuple{CircuitState, Float64}}
    new_circuit = vcat(state.circuit, [action])
    reward = 0
    initial_quantum_state = copy(state.quantum_state)
    gate_counts = copy(state.gate_counts)
    if action isa AbstractSingleQubitOperator
        reward -= mdp.mcts_params.fitness_weights[2]
        gate_counts[1] += 1
        qubit = action.q
        new_quantum_state = execute_circuit(AbstractOperation[sHadamard(qubit)], initial_quantum_state)
    elseif action isa AbstractTwoQubitOperator
        control = action.q1
        target = action.q2
        control_register = mdp.network_specs.register_lookup_array[mdp.network_specs.inv_map[control]] 
        target_register = mdp.network_specs.register_lookup_array[mdp.network_specs.inv_map[target]] 
        if control_register == target_register 
            reward -= mdp.mcts_params.fitness_weights[3]
            gate_counts[2] += 1
        else
            reward -= mdp.mcts_params.fitness_weights[4]
            gate_counts[3] += 1
        end
        new_quantum_state = execute_circuit(AbstractOperation[sCNOT(control, target)], initial_quantum_state)
    end
    new_quantum_state_tab = tab(canonicalize_rref!( stabilizerview(new_quantum_state) )[1])
    new_quantum_state_bit_matrix = tableau_to_bitmatrix(new_quantum_state_tab) 
    tab_distance = tableau_distance(new_quantum_state_bit_matrix, mdp.code_params.target_bit_matrix, metric = mdp.mcts_params.tableau_metric)
    fidelity = 1 - tab_distance 
    reward += mdp.mcts_params.fitness_weights[1]*(fidelity-state.fidelity)      
    new_state = CircuitState(new_circuit, copy(new_quantum_state), new_quantum_state_bit_matrix, gate_counts, fidelity) 
    return (sp=new_state, r=reward) 
end


"""
    execute_circuit(circuit::Vector{AbstractOperation}, 
                        initial_state::MixedDestabilizer{Tableau{Vector{UInt8}, Matrix{UInt64}}})::MixedDestabilizer

Execute the given sequence of gates applied to a specified `initial_state`.

### Input
- `circuit` -- the specified sequence of quantum gates to apply
- `initial_state` -- the state onto which the circuit is applied

### Output
The final state of the quantum system, captured as `MixedDestabilizer` object.

### Notes
For circuit execution, we leverage the Monte Carlo trajectory simulation
function `mctrajectory!` provided by `QuantumClifford`.
"""
function execute_circuit(circuit::Vector{AbstractOperation}, 
                        initial_state::MixedDestabilizer{Tableau{Vector{UInt8}, Matrix{UInt64}}})::MixedDestabilizer
    initial_state = Register(initial_state, 0)
    state, stat = mctrajectory!(copy(initial_state), circuit)
    return state.stab
end


``` Hash overwrites to discern equality of actions and states```
function Base.hash(s::CircuitState, h::UInt)
    hash(s.circuit, h)
end
function Base.:(==)(s1::CircuitState, s2::CircuitState)
    length(s1.circuit) == length(s2.circuit) && all(s1.circuit .== s2.circuit)
end
Base.hash(g::sHadamard, h::UInt) = hash((:H, g.q), h)
Base.hash(g::sCNOT, h::UInt)    = hash((:CX, g.q1, g.q2), h)
Base.:(==)(g1::sCNOT, g2::sCNOT) = g1.q1 == g2.q1 && g1.q2 == g2.q2



end