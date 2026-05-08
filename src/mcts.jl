module MonteCarloTreeSearch

using ..Types
#using ..CircuitSimulator
using ..Helper
#using ..BaselineEncoding

using POMDPs, POMDPTools
using MCTS

using Random
using Serialization, CSV, DataFrames

#rng = MersenneTwister(42)
#using Quantikz: savecircuit, @with, classicalbitslayout
#using QECCore
using QuantumClifford
using QuantumClifford: MixedDestabilizer, Tableau, AbstractOperation, sHadamard, sCNOT, sSWAP, @S_str, true_success_stat, false_success_stat, continue_stat, failure_stat, PauliMeasurement, VerifyOp
# using BenchmarkTools
using D3Trees

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
    #s1.bit_matrix == s2.bit_matrix # only comparing bit matrices is not sufficient since two circuits can have the same action on the zero state, but not the same when another gate is appended!
    #length(s1.gates) == length(s2.gates) && all(s1.gates .== s2.gates)
    s1.circuit == s2.circuit
end

#Base.hash(g::HadamardGate, h::UInt) = hash((:H, g.index), h)

Base.hash(g::sHadamard, h::UInt) = hash((:H, g.q), h)

# Base.hash(g::SGate, h::UInt) = hash((:S, g.index), h)
# Base.hash(g::InvSGate, h::UInt) = hash((:iS, g.index), h)
# Base.hash(g::SqrtXGate, h::UInt) = hash((:sqX, g.index), h)
# Base.hash(g::InvSqrtXGate, h::UInt) = hash((:isqX, g.index), h)
# Base.hash(g::PauliXGate, h::UInt) = hash((:X, g.index), h)
# Base.hash(g::PauliYGate, h::UInt) = hash((:Y, g.index), h)
# Base.hash(g::PauliZGate, h::UInt) = hash((:Z, g.index), h)

# Base.hash(g::CX_Gate, h::UInt)    = hash((:CX, g.control, g.target), h)
# Base.:(==)(a::CX_Gate, b::CX_Gate) = a.control == b.control && a.target == b.target

Base.hash(g::sCNOT, h::UInt)    = hash((:CX, g.q1, g.q2), h)
Base.:(==)(g1::sCNOT, g2::sCNOT) = g1.q1 == g2.q1 && g1.q2 == g2.q2

# Base.hash(g::CZ_Gate, h::UInt)    = hash((:CZ, g.control, g.target), h)
# Base.:(==)(a::CZ_Gate, b::CZ_Gate) = a.control == b.control && a.target == b.target


struct EncodingMDP <: MDP{CircuitState, AbstractOperation}  # Abstract Type of state S and action A (circuit and gate, respectively)
    code_params:: CodeParameters
    network_specs:: NetworkSpecifications
    mcts_params:: MCTSParameters
end

function POMDPs.actions(mdp::EncodingMDP, s::CircuitState)
   
    n = mdp.network_specs.num_data_qubits
    actions = Vector{AbstractOperation}()
    #single_qubit_gates = mdp.mcts_params.gate_set.single_qubit_gates
    #two_qubit_gates = mdp.mcts_params.gate_set.two_qubit_gates
    affected_qubits = [(g isa AbstractTwoQubitOperator ? g.q2 : g.q) for g in s.circuit]
    #print(affected_qubits)
    #print(length(s.gates), mdp.code_params.num_X_checks )
    if length(s.circuit) < mdp.code_params.num_X_checks # we add as many H as the number of X stabilisers (assuming blank start)
        #push!(actions, HadamardGate(1))
        #push!(actions, HadamardGate(4))
        for i in 1:n
            if length(s.circuit)==0
                push!(actions, sHadamard(i))
            elseif i ∉ affected_qubits
                push!(actions, sHadamard(i))
            end
            # for gate in single_qubit_gates # in this case only Hadamard
            #     if length(s.gates)==0
            #         push!(actions, gate(i))
            #     elseif i ∉ affected_qubits
            #         push!(actions, gate(i))
            #     end
            # end 
        end
    else 
        for c in 1:n, t in 1:n
            c == t && continue
            #for gate in two_qubit_gates # in this case only CNOT
                #println(last(s.gates))
            if last(s.circuit) isa AbstractSingleQubitOperator
                if c ∈ affected_qubits
                    push!(actions, sCNOT(c, t))
                end
            elseif !( (last(s.circuit).q1 == c && last(s.circuit).q2 ==t) || (c ∉ affected_qubits)) 
                # if (typeof(last(s.gates)) <: TwoQubitGate) && (last(s.gates).control == 4 && last(s.gates).target == 5) && (length(s.gates)>6) && (s.fidelity >0.93)
                #     println("Current gates: $(s.gates)")
                #     println("Possible action: $(gate(c,t))")
                # end
                #println("Current gates: $(s.gates)")
                #println("Possible action: $(gate(c,t))")
                push!(actions, sCNOT(c, t))
            end
            #end
        end
    end
    # if length(s.gates) ==1
    #     print("Actions at step 1 for $(s.gates): $actions \n")
    # end
    return actions
end


# TODO: 
# - Apply penalty for appling same gate again,
# - then input MCTS > GA

#  apply reward for attaining corect stabiliser or logicla operator, normalise reward terms 

function POMDPs.gen(mdp::EncodingMDP, state::CircuitState, action::AbstractOperation, rng)
    # We use the generative interface since the transition from current state and action to next
    # state is trivial, even though this technically does not match the recommended use case for gen

    new_circuit = vcat(state.circuit, [action])
    #print("Action is $action")
    # sp stands for s', the next state; it is the circuit obtained by appending the gate (= action) to the current state s ( = circuit)

    reward = 0 # (1e-8)*rand(rng)

    initial_quantum_state = copy(state.quantum_state)
    gate_counts = copy(state.gate_counts)
    #println(typeof(initial_quantum_state))
    #T = typeof(action)
    if action isa AbstractSingleQubitOperator# T <: SingleQubitGate #in mdp.gate_set.single_qubit_gates# isa Union{PauliXGate, PauliYGate, PauliZGate, HadamardGate, SGate} 
        reward -= mdp.mcts_params.fitness_weights[2]
        gate_counts[1] += 1#[1,0,0]
        qubit = action.q

        new_quantum_state = execute_circuit([sHadamard(qubit)], initial_quantum_state)
        #push!(circuit, gate_to_apply(T, n.inv_perm[qubit]) ) 

    elseif action isa AbstractTwoQubitOperator#T <: TwoQubitGate#  in mdp.gate_set.two_qubit_gates 
        control = action.q1
        target = action.q2
        control_register = mdp.network_specs.register_lookup_array[mdp.network_specs.inv_map[control]] 
        target_register = mdp.network_specs.register_lookup_array[mdp.network_specs.inv_map[target]] 
        
        if control_register == target_register # the lookup array does not account for the communication qubits
            reward -= mdp.mcts_params.fitness_weights[3]
            gate_counts[2] += 1#[0,1,0]
            #push!(circuit, gate_to_apply(T, n.comm_inv_perm_idx[control], n.comm_inv_perm_idx[target] ))
        else
            reward -= mdp.mcts_params.fitness_weights[4]
            gate_counts[3] += 1#[0,0,1]
        end
        
        new_quantum_state = execute_circuit([sCNOT(control, target)], initial_quantum_state)
    end
    # Build the executable DQC circuit (handles telegates, comm qubits, etc.)
    #quantum_clifford_circuit, num_single_qubit_gates, num_two_qubit_gates, num_telegates = construct_executable_circuit(new_circuit, mdp.gate_set, mdp.network_specs)        
    #gate_counts = (num_single_qubit_gates, num_two_qubit_gates, num_telegates)

    # Fidelity in (0,1)
    # Reward 1/n in (0,1) for each correct mapping of stabilisers or logicals
    # Circuit Sizes each normalised to (0,1)
    # very large reward for correct state
    #quantum_clifford_circuit, num_single_qubit_gates, num_two_qubit_gates, num_telegates = construct_executable_circuit(new_state.gates, mdp.gate_set, mdp.network_specs)        
    #mc_result = execute_circuit(quantum_clifford_circuit, mdp.network_specs.num_qubits, mdp.network_specs.num_registers; num_traj= mdp.network_specs.num_shots)#, keepstates = true) # if specifying num_traj, we use MC sampling, otherwise perturbation.
    #println(new_quantum_state)
    #new_quantum_state = only(new_quantum_state)
    new_quantum_state_tab = tab(canonicalize_rref!( stabilizerview(new_quantum_state) )[1])
    #print(new_quantum_state)
    #stab_view = traceout!(copy(stab_view), mdp.network_specs.comm_qubits) # TODO: This can be refactored to ptrace upon stable QS release
    # NOTE: if we swtich to ptrace, then also tableau_distance in the helper.jl needs to be adapted!
    #stab_canon = canonicalize_rref!( new_quantum_state )
    #tableau = tab(stab_canon[1])
    #println("Tableau: $tableau")
    new_quantum_state_bit_matrix = tableau_to_bitmatrix(new_quantum_state_tab) # extract the stabiliser tableau from MixedDestabilizer object
    tab_distance = tableau_distance(new_quantum_state_bit_matrix, mdp.code_params.target_bit_matrix, metric = mdp.mcts_params.tableau_metric)#, mdp.network_specs.data_qubits, mdp.network_specs.comm_qubits, mdp.opt_params.tableau_metric)
    #println(tab_distance)
    fidelity = 1 - tab_distance # 1 is perfect alignment, here we are in the noiseless setting (one shot)
    #circuit_size(quantum_clifford_circuit) #  length(quantum_clifford_circuit)
    reward += mdp.mcts_params.fitness_weights[1]*(fidelity-state.fidelity)# - sum(mdp.mcts_params.fitness_weights[2:4] .* gate_counts)    # via the discount factor, large depth will be penalised

    #reward += mdp.mcts_params.fitness_weights[1]*(fidelity-state.fidelity)

    # q_state = one(Stabilizer, mdp.code_params.n)
    # #     q_state = copy(init_state)
    # #     for gate in new_circuit
    # #         T = typeof(gate)
    # #         apply!(q_state, gate_to_apply(T, gate.index))
    # #     end
    # #     if q_state in mdp.code_params.stabilizer_group,
    # # end 
    # correct_stabs = 0
    # #println("Final Tablea: $q_state")
    # #println("Final Tablea: $(q_state[1])")
    # for i in eachindex(q_state)
    #     if q_state[i] in vcat(mdp.code_params.stabilizer_group, mdp.code_params.logical_Zs)
    #         #reward += 1/mdp.code_params.n
    #         correct_stabs += 1
    #     end
    # end
    
    # if correct_stabs == mdp.code_params.n
    #     reward += 1e2
    # end

    #gate_counts_norm = gate_counts ./ sum(gate_counts)

    #weights_norm = mdp.mcts_params.fitness_weights[2:4] / sum(mdp.mcts_params.fitness_weights[2:4])

    #println("Gate counts norm: $gate_counts_norm \n Weights norm: $weights_norm")

    #circuit_size = sum(weights_norm .* gate_counts_norm )

    #reward -= circuit_size
    #reward += 1e-6*rand(rng)
    #     if stabilizerview(q_state)[1] in mdp.code_params.stabilizer_group
    #         reward += 0.25
    #     end
    # end
    # for i in eachindex(mdp.code_params.logical_Zs)
    #     q_state = Stabilizer(QuantumClifford.Tableau([mdp.code_params.logical_Zs[i]]))
    #     for gate in new_circuit
    #         T = typeof(gate)
    #         apply!(q_state, gate_to_apply(T, gate.index))
    #     end
    #     if stabilizerview(q_state)[1] in mdp.code_params.target_logical_Zs[i]
    #         reward += 0.25
    #     end
    # end
    # Potential-based shaping: γΦ(s') - Φ(s), where Φ = -distance
    # This preserves the optimal policy while providing dense signal
    
    new_state = CircuitState(new_circuit, copy(new_quantum_state), new_quantum_state_bit_matrix, gate_counts, fidelity) 
    return (sp=new_state, r=reward) # state.fidelity  gc for logging purposes

    
end

POMDPs.discount(mdp::EncodingMDP) = mdp.mcts_params.discount_factor

function POMDPs.isterminal(mdp::EncodingMDP, s::CircuitState)
    # Terminate at max depth or if already perfectly encoding

    #if mdp.code_params isa CodeParameters
    if s.fidelity >= 1.0 # fidelity
        return true
    end
    # elseif mdp.code_params isa CodeParametersLogical
    #     if s.metric == length(mdp.code_params.stabilizer_generators) + 20*length(mdp.code_params.logical_Xs)# shaped reward,  mdp.code_params.n + mdp.code_params.k # n+k checks were successful 
    #         return true
    #     end
    #end
        # Optional: early termination on perfect fidelity (distance == 0)
    # This requires re-evaluating the tableau, so only worth doing if 
    # you cache the last evaluated distance in state (see extension below)
    return false
end

function POMDPs.initialstate(mdp::EncodingMDP)
    # Warm start: seed with first half of compiler circuit, matching your GA
    # Or start from empty circuit:
    # can encode some warm start circuit here (for example, if run terminated too early)

    init_state = one(MixedDestabilizer,mdp.network_specs.num_data_qubits)
    init_bit_matrix = tableau_to_bitmatrix( tab(canonicalize_rref!( stabilizerview(init_state) )[1])  ) 
    return CircuitState([], init_state, init_bit_matrix, [0,0,0], 0.0)
end



# function estimate_value_f(mdp, s, remaining_depth) #POMPDs typing requirement
#     return s.fidelity
# end

function monte_carlo_tree_search(code_params, network_specs, mcts_params, folder)
    
    @info "Monte Carlo Tree Search with depth $(mcts_params.depth) and $(mcts_params.n_iterations) iterations started..."

    #rng = MersenneTwister(rand(UInt))

    mdp = EncodingMDP(
        code_params, network_specs, mcts_params
    )
    
    
    
    #solver = DPWSolver(n_iterations=mdp.mcts_params.n_iterations, depth=mdp.mcts_params.depth, exploration_constant=mdp.mcts_params.exploration_constant, alpha_state=1/8, tree_in_info=false)
    solver = MCTSSolver(
        n_iterations = mdp.mcts_params.n_iterations,#mdp.mcts_params.n_iterations, # if no solution is found, this should be increased
        depth = mcts_params.depth,# mdp.mcts_params.depth, # similar to above
        exploration_constant = mcts_params.exploration_constant,#mdp.mcts_params.exploration_constant, # should be increased if search space is not explored well
        rng = Random.GLOBAL_RNG,
        reuse_tree = false,
        enable_tree_vis = false,
        estimate_value = 0.0# estimate_value_f # RolloutEstimator(RandomSolver(Random.GLOBAL_RNG,))#0.0
    )
    
    policy = solve(solver, mdp)

    s = initialstate(mdp)# CircuitState(Gate[])
    MCTS_gate_counts = (typemax(Int), typemax(Int), typemax(Int))
    MCTS_fidelity = typemin(Float16)
    MCTS_circuit_state = s

    gate_count_evolution = Vector{Vector{Int64}}()
    fidelity_evolution = Float64[]
    reward_evolution = Float64[]

    steps = 0
    while steps < mdp.mcts_params.max_steps 
       
        a = action(policy, s)
        # tree = policy.tree
        # if steps ==1
        #     # sns = collect(MCTS.state_nodes(tree))
        #     #     # Root is the state node whose outgoing action nodes have the largest total visits
        #     #     root_sn = sns[argmax([sum(MCTS.n(san) for san in MCTS.children(sn)) for sn in sns])]
        #     #     root_children = collect(MCTS.children(root_sn))
        #     #     sort!(root_children, by = san -> MCTS.n(san), rev = true)

        #     #     println("Root action stats:")
        #     #     for san in root_children
        #     #         println("a:$(action(san)) Q:$(MCTS.q(san)) N:$(MCTS.n(san))")
        #     #     end
        #     for sn in MCTS.state_nodes(tree)
        #         #@info " State node is $sn"
        #         for san in MCTS.children(sn)
        #             @info "s:$(action(san)) Q:$(MCTS.q(san)), N:$(MCTS.n(san))"
        #         end
        #     end
        # end
        
        #if mdp.code_params isa CodeParameters
        s, r = POMDPs.gen(mdp, s, a, Random.GLOBAL_RNG)

        @info "Applied $(a), reward=$r for fidelity=$(s.fidelity) and gate_counts=$(s.gate_counts), DQC circuit length=$(length(s.circuit))"
        MCTS_circuit_state = s
        MCTS_gate_counts = s.gate_counts
        MCTS_fidelity = s.fidelity

        push!(gate_count_evolution, MCTS_gate_counts)
        push!(fidelity_evolution, MCTS_fidelity)
        push!(reward_evolution, r)

        if POMDPs.isterminal(mdp, s)#, fidelity=fidelity) 
            println("Terminal condition reached after $(steps+1) steps. Final fidelity: $(s.fidelity).")
            break
        end
        steps += 1
    end
    

    MCTS_circuit = MCTS_circuit_state.circuit
    @info "MCTS optimised circuit length in DQC setting:: Single-qubit gates: $(MCTS_gate_counts[1]), Two-qubit gates: $(MCTS_gate_counts[2]), Telegates: $(MCTS_gate_counts[3])"
    
    # ----- Verification ------

    verification_logical_state = verify_success(copy(MCTS_circuit), mdp.code_params.target_state, mdp.network_specs)
    # ^NOTE: this appends a verifyop operation, but we count before so this is irrelevant
    #verification_logical_state = verification_logical_state == 1.0 ? true : false
    @info "Verification of MCTS circuit successful: $verification_logical_state"


    # ----- Data Storage ----------

    #a, info = action_info(policy, MCTS_circuit_state)
    #print(info[:tree])
    #inchrome(D3Tree(info[:tree]))

    #dir = joinpath(folder, "MCTS")
    #MCTS_dir = next_run_dir(folder) # accounting for multiple MCTS runs in the same folder (probabilistic algorithm yields dffernet outcomes for same hyperparameters)
    MCTS_dir = folder

    @info "Saving results of MCTS run to $(MCTS_dir)"

    serialize( joinpath(MCTS_dir, "MCTS_gates.jls"), MCTS_circuit )

    # open(joinpath(MCTS_dir, "MCTS_circuit.txt"), "w") do io
    #     println(io, "# Raw gate sequence of size $(sum(MCTS_gate_counts))")
    #     for (i, g) in enumerate(MCTS_circuit)
    #         println(io, i, "\t", repr(g))
    #     end
    # end
    
    save_circuit_diagram(copy(MCTS_circuit), MCTS_dir, "MCTS_circuit.png")

    open(joinpath(MCTS_dir, "MCTS_params.txt"), "w") do io
        println(io, "MCTS parameters")
        for fname in fieldnames(Types.MCTSParameters)
            println(io, fname, " = ", repr(getfield(mcts_params, fname)))
        end
    end
            
    open(joinpath(MCTS_dir, "summary.txt"), "w") do io
        println(io, "# Encoding successful: $verification_logical_state")
        println(io, "# Raw gate sequence of size $(sum(MCTS_gate_counts))")
        println(io, "# Executable DQC circuit with $(MCTS_gate_counts[1]) single qubit gates, $(MCTS_gate_counts[2]) two qubit gates and $(MCTS_gate_counts[3]) telegates ")
    end
    
    data = DataFrame(fidelity_evolution = fidelity_evolution, 
                   gate_count_evolution = gate_count_evolution,
                   reward_evolution = reward_evolution)

    CSV.write(joinpath(MCTS_dir, "evolution.csv"), data)
    plot_evolution(MCTS_dir, "MCTS", fidelity_evolution, gate_count_evolution, mcts_params, verification_logical_state)

    return verification_logical_state, MCTS_gate_counts#, MCTS_gates, MCTS_dir
end


end