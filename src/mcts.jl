module MonteCarloTreeSearch

using ..Types
using ..Helper

using POMDPs, POMDPTools
using MCTS
using Random
using Serialization, CSV, DataFrames
using ProgressMeter
using QuantumClifford
using QuantumClifford: MixedDestabilizer, Tableau, AbstractOperation, sHadamard, sCNOT, sSWAP, @S_str, true_success_stat, false_success_stat, continue_stat, failure_stat, PauliMeasurement, VerifyOp

export monte_carlo_tree_search

struct CircuitState
    circuit::Vector{AbstractOperation}   # raw Gate[] as in CircuitIndividual, no comm qubit indexing
    quantum_state::MixedDestabilizer{Tableau{Vector{UInt8}, Matrix{UInt64}}}   # the corresponding state after circuit execution, starting from the all zero state
    bit_matrix::Matrix{Int}
    gate_counts::Vector{Int}
    fidelity::Float64
end

function Base.hash(s::CircuitState, h::UInt)
    hash(s.circuit, h)
end
function Base.:(==)(s1::CircuitState, s2::CircuitState)
    # only comparing bit matrices is not sufficient since two circuits can have the same action on the zero state, but not the same when another gate is appended!
    length(s1.circuit) == length(s2.circuit) && all(s1.circuit .== s2.circuit)
end


Base.hash(g::sHadamard, h::UInt) = hash((:H, g.q), h)

Base.hash(g::sCNOT, h::UInt)    = hash((:CX, g.q1, g.q2), h)
Base.:(==)(g1::sCNOT, g2::sCNOT) = g1.q1 == g2.q1 && g1.q2 == g2.q2


struct EncodingMDP <: MDP{CircuitState, AbstractOperation}  # Abstract Type of state S and action A (circuit and gate, respectively)
    code_params:: CodeParameters
    network_specs:: NetworkSpecifications
    mcts_params:: MCTSParameters
end

function POMDPs.actions(mdp::EncodingMDP, s::CircuitState)
   
    n = mdp.network_specs.num_data_qubits
    actions = Vector{AbstractOperation}()
    affected_qubits = [(g isa AbstractTwoQubitOperator ? g.q2 : g.q) for g in s.circuit]
    if length(s.circuit) < mdp.code_params.num_X_checks # we add as many H as the number of X stabilisers (assuming blank start)
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

function POMDPs.gen(mdp::EncodingMDP, state::CircuitState, action::AbstractOperation, rng)
    # We use the generative interface since the transition from current state and action to next
    # state is trivial, even though this technically does not match the recommended use case for gen

    new_circuit = vcat(state.circuit, [action])
    # sp stands for s', the next state; it is the circuit obtained by appending the gate (= action) to the current state s ( = circuit)

    reward = 0

    initial_quantum_state = copy(state.quantum_state)
    gate_counts = copy(state.gate_counts)

    if action isa AbstractSingleQubitOperator
        reward -= mdp.mcts_params.fitness_weights[2]
        gate_counts[1] += 1
        qubit = action.q

        new_quantum_state = execute_circuit([sHadamard(qubit)], initial_quantum_state)

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
        
        new_quantum_state = execute_circuit([sCNOT(control, target)], initial_quantum_state)
    end
    
    new_quantum_state_tab = tab(canonicalize_rref!( stabilizerview(new_quantum_state) )[1])
    new_quantum_state_bit_matrix = tableau_to_bitmatrix(new_quantum_state_tab) # extract the stabiliser tableau from MixedDestabilizer object
    tab_distance = tableau_distance(new_quantum_state_bit_matrix, mdp.code_params.target_bit_matrix, metric = mdp.mcts_params.tableau_metric)
    fidelity = 1 - tab_distance # 1 is perfect alignment, here we are in the noiseless setting (one shot)
    reward += mdp.mcts_params.fitness_weights[1]*(fidelity-state.fidelity)  # via the discount factor, large depth will be penalised       
    new_state = CircuitState(new_circuit, copy(new_quantum_state), new_quantum_state_bit_matrix, gate_counts, fidelity) 

    return (sp=new_state, r=reward) 
end

POMDPs.discount(mdp::EncodingMDP) = mdp.mcts_params.discount_factor

function POMDPs.isterminal(mdp::EncodingMDP, s::CircuitState)
    # Terminate early if fidelity = 1
    if s.fidelity >= 1.0 
        return true
    end
    return false
end

function POMDPs.initialstate(mdp::EncodingMDP)
    init_state = one(MixedDestabilizer,mdp.network_specs.num_data_qubits)
    init_bit_matrix = tableau_to_bitmatrix( tab(canonicalize_rref!( stabilizerview(init_state) )[1])  ) 
    return CircuitState([], init_state, init_bit_matrix, [0,0,0], 0.0)
end


function monte_carlo_tree_search(code_params, network_specs, mcts_params)
    
    @info "Initialised Monte Carlo Tree Search with depth $(mcts_params.depth) and $(mcts_params.n_iterations) iterations ..."

    mdp = EncodingMDP(code_params, network_specs, mcts_params)
        
    solver = MCTSSolver(
        n_iterations = mdp.mcts_params.n_iterations,# if no solution is found, this should be increased
        depth = mcts_params.depth, # "
        exploration_constant = mcts_params.exploration_constant,# "
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

    #----------------------------- MCTS steps --------------------------------
    p = Progress(mdp.mcts_params.max_steps; desc = "Monte Carlo Tree search", showspeed = true)
    
    for step in 1:mdp.mcts_params.max_steps 
        a = action(policy, s)    
        s, r = POMDPs.gen(mdp, s, a, Random.GLOBAL_RNG)

        #@info "Applied action $(a)"
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
    
    # ----- Verification ------
    verification_logical_state = verify_success(copy(MCTS_circuit), mdp.code_params.target_state, mdp.network_specs)
    @info "Verification of MCTS circuit successful: $verification_logical_state"

    # ----- Tree Display ----------
    #a, info = action_info(policy, MCTS_circuit_state)
    #inchrome(D3Tree(info[:tree]))

    return MCTS_circuit, verification_logical_state, MCTS_gate_counts, fidelity_evolution, gate_count_evolution, reward_evolution
end


end