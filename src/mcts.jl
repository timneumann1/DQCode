module MonteCarloTreeSearch

using ..Types
using ..CircuitSimulator
using ..Helper
using ..LogicalEnc

using POMDPs, POMDPTools
using MCTS

using Random
#rng = MersenneTwister(42)
#using Quantikz: savecircuit, @with, classicalbitslayout
#using QECCore
#using QECCore: Steane7, QuantumTannerGraphProduct, CyclicQuantumTannerGraphProduct, Triangular488
# using QuantumClifford.ECC: DistanceMIPAlgorithm
# using HiGHS
using QuantumClifford
using QuantumClifford: MixedDestabilizer, sHadamard, sCNOT, sSWAP, @S_str, true_success_stat, false_success_stat, continue_stat, failure_stat, PauliMeasurement, VerifyOp
# using BenchmarkTools
# using CairoMakie
# using KaHyPar
# using SparseArrays
using D3Trees

export monte_carlo_tree_search

struct CircuitState
    gates::Vector{Gate}   # raw Gate[] as in CircuitIndividual, no comm qubit indexing
    fidelity::Float64
end

function Base.hash(s::CircuitState, h::UInt)
    hash(s.gates, h)
end
function Base.:(==)(a::CircuitState, b::CircuitState)
    length(a.gates) == length(b.gates) && all(a.gates .== b.gates)
end

Base.hash(g::HadamardGate, h::UInt) = hash((:H, g.index), h)
Base.hash(g::SGate, h::UInt) = hash((:S, g.index), h)
Base.hash(g::InvSGate, h::UInt) = hash((:iS, g.index), h)
Base.hash(g::SqrtXGate, h::UInt) = hash((:sqX, g.index), h)
Base.hash(g::InvSqrtXGate, h::UInt) = hash((:isqX, g.index), h)
Base.hash(g::PauliXGate, h::UInt) = hash((:X, g.index), h)
Base.hash(g::PauliYGate, h::UInt) = hash((:Y, g.index), h)
Base.hash(g::PauliZGate, h::UInt) = hash((:Z, g.index), h)

Base.hash(g::CX_Gate, h::UInt)    = hash((:CX, g.control, g.target), h)
Base.:(==)(a::CX_Gate, b::CX_Gate) = a.control == b.control && a.target == b.target

Base.hash(g::CZ_Gate, h::UInt)    = hash((:CZ, g.control, g.target), h)
Base.:(==)(a::CZ_Gate, b::CZ_Gate) = a.control == b.control && a.target == b.target


struct EncodingMDP <: MDP{CircuitState, Gate}  # Abstract Type of state S and action A (circuit and gate, respectively)
    discount_factor::Float64
    code_params:: Union{CodeParameters, CodeParametersLogical}
    network_specs:: NetworkSpecifications
    opt_params:: OptimisationParameters
    mcts_params:: MCTSParameters
    gate_set:: GateSet
end

function POMDPs.actions(mdp::EncodingMDP, s::CircuitState)
    n = mdp.network_specs.num_data_qubits
    actions = Gate[]
    single_qubit_gates = mdp.gate_set.single_qubit_gates
    two_qubit_gates = mdp.gate_set.two_qubit_gates
    affected_qubits = affected_qubits = [q for g in s.gates for q in (g isa CX_Gate || g isa CZ_Gate ? (g.control, g.target) : (g.index,))]
    #print(affected_qubits)
    #print(length(s.gates), mdp.code_params.num_X_checks )
    if length(s.gates) < mdp.code_params.num_X_checks # we add as many H as the number of X stabilisers (assuming blank start)
        for i in 1:n
            for gate in single_qubit_gates # in this case only Hadamard
                if length(s.gates)==0
                    push!(actions, gate(i))
                elseif last(s.gates).index != i
                    push!(actions, gate(i))
                end
            end 
        end
    else 
        for c in 1:n, t in 1:n
            c == t && continue
            for gate in two_qubit_gates # in this case only CNOT
                #println(last(s.gates))
                if typeof(last(s.gates)) in single_qubit_gates
                    push!(actions, gate(c, t))
                elseif !( (last(s.gates).control == c && last(s.gates).target ==t) || (c ∉ affected_qubits && t ∉ affected_qubits)) 
                    push!(actions, gate(c, t))
                end
            end
        end
    end
    return actions
end


# TODO: 
# - Apply penalty for appling same gate again,
# - then input MCTS > GA

#  apply reward for attaining corect stabiliser or logicla operator, normalise reward terms 

function POMDPs.gen(mdp::EncodingMDP, state::CircuitState, action::Gate, rng)
    # We use the generative interface since the transition from current state and action to next
    # state is trivial, even though this technically does not match the recommended use case for gen
    
    new_circuit = vcat(state.gates, [action])
    print("Action is $action")
    # sp stands for s', the next state; it is the circuit obtained by appending the gate (= action) to the current state s ( = circuit)

    # Build the executable DQC circuit (handles telegates, comm qubits, etc.)
    quantum_clifford_circuit, num_single_qubit_gates, num_two_qubit_gates, num_telegates = construct_executable_circuit(new_circuit, mdp.gate_set, mdp.network_specs)        
    gate_counts = (num_single_qubit_gates, num_two_qubit_gates, num_telegates)

    if mdp.code_params isa CodeParameters


        # Fidelity in (0,1)
        # Reward 1/n in (0,1) for each correct mapping of stabilisers or logicals
        # Circuit Sizes each normalised to (0,1)
        # very large reward for correct state
        #quantum_clifford_circuit, num_single_qubit_gates, num_two_qubit_gates, num_telegates = construct_executable_circuit(new_state.gates, mdp.gate_set, mdp.network_specs)        
        mc_result = execute_circuit(quantum_clifford_circuit, mdp.network_specs.num_qubits, mdp.network_specs.num_registers; num_traj= mdp.network_specs.num_shots)#, keepstates = true) # if specifying num_traj, we use MC sampling, otherwise perturbation.
        
        stab_view = stabilizerview(only(mc_result))
        stab_view = traceout!(copy(stab_view), mdp.network_specs.comm_qubits) # TODO: This can be refactored to ptrace upon stable QS release
        # NOTE: if we swtich to ptrace, then also tableau_distance in the helper.jl needs to be adapted!
        stab_canon = canonicalize_rref!( stab_view )
        tableau = tab(stab_canon[1])
        #println("Tableau: $tableau")
        current_bit_matrix = tableau_to_bitmatrix(tableau) # extract the stabiliser tableau from MixedDestabilizer object
        tab_distance = tableau_distance(current_bit_matrix, mdp.code_params.target_bit_matrix, mdp.network_specs.data_qubits, mdp.network_specs.comm_qubits, mdp.opt_params.tableau_metric)
        fidelity = 1 - tab_distance # 1 is perfect alignment, here we are in the noiseless setting (one shot)
         #circuit_size(quantum_clifford_circuit) #  length(quantum_clifford_circuit)
        
        reward = mdp.mcts_params.fitness_weights[1]*fidelity

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

        gate_counts_norm = gate_counts ./ sum(gate_counts)

        weights_norm = mdp.mcts_params.fitness_weights[2:4] / sum(mdp.mcts_params.fitness_weights[2:4])
    
        #println("Gate counts norm: $gate_counts_norm \n Weights norm: $weights_norm")

        circuit_size = sum(weights_norm .* gate_counts_norm )

        reward -= circuit_size
        reward += 1e-6*rand(rng)
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
        #reward = mdp.mcts_params.fitness_weights[1]*fidelity - sum(mdp.mcts_params.fitness_weights[2:4] .* gate_counts) + 1e-6*rand(rng)   # via the discount factor, large depth will be penalised
        
        new_state = CircuitState(new_circuit, fidelity) 
        return (sp=new_state, r=reward, gc=gate_counts ) # state.fidelity


    # elseif mdp.code_params isa CodeParametersLogical
    #     #println("CIRCUIT: $quantum_clifford_circuit")
    #     #println("new_circuit:$new_circuit")
    #     reward = 0.0

    #     # TODO: Parellelis the below, find intermediate reward?
    #     for g in mdp.code_params.stabilizer_generators
    #         q_state = Stabilizer(QuantumClifford.Tableau([g]))
    #         for gate in new_circuit
    #             T = typeof(gate)
    #             apply!(q_state, gate_to_apply(T, gate.index))
    #         end
    #         if stabilizerview(q_state)[1] in mdp.code_params.stabilizer_group
    #             reward += 0.1
    #         end
    #     end
    #     for i in eachindex(mdp.code_params.logical_Xs)
    #         q_state = Stabilizer(QuantumClifford.Tableau([mdp.code_params.logical_Xs[i]]))
    #         for gate in new_circuit
    #             T = typeof(gate)
    #             apply!(q_state, gate_to_apply(T, gate.index))
    #         end
    #         if stabilizerview(q_state)[1] in mdp.code_params.target_logical_Xs[i]
    #             reward += 5.0
    #         end
    #     end
    #     for i in eachindex(mdp.code_params.logical_Zs)
    #         q_state = Stabilizer(QuantumClifford.Tableau([mdp.code_params.logical_Zs[i]]))
    #         for gate in new_circuit
    #             T = typeof(gate)
    #             apply!(q_state, gate_to_apply(T, gate.index))
    #         end
    #         if stabilizerview(q_state)[1] in mdp.code_params.target_logical_Zs[i]
    #             reward += 5.0
    #         end
    #     end
    #     #reward = mdp.mcts_params.fitness_weights[1]*reward - sum(mdp.mcts_params.fitness_weights[2:4] .* gate_counts) + 1e-6*rand(rng)
    #     new_state = CircuitState(new_circuit, reward) 
    #     return (sp=new_state, r=reward-state.metric, quality = reward, gc=gate_counts )
    end
    
end

POMDPs.discount(mdp::EncodingMDP) = mdp.discount_factor

function POMDPs.isterminal(mdp::EncodingMDP, s::CircuitState)
    # Terminate at max depth or if already perfectly encoding

    if mdp.code_params isa CodeParameters
        if s.fidelity >= 1.0 # fidelity
            return true
        end
    # elseif mdp.code_params isa CodeParametersLogical
    #     if s.metric == length(mdp.code_params.stabilizer_generators) + 20*length(mdp.code_params.logical_Xs)# shaped reward,  mdp.code_params.n + mdp.code_params.k # n+k checks were successful 
    #         return true
    #     end
    end
        # Optional: early termination on perfect fidelity (distance == 0)
    # This requires re-evaluating the tableau, so only worth doing if 
    # you cache the last evaluated distance in state (see extension below)
    return false
end

function POMDPs.initialstate()
    # Warm start: seed with first half of compiler circuit, matching your GA
    # Or start from empty circuit:
    # can encode some warm start circuit here (for example, if run terminated too early)
    return CircuitState([], 0.0)
end


function monte_carlo_tree_search(code_params, network_specs, opt_params, mcts_params, gate_set)
    
    mdp = EncodingMDP(
        mcts_params.discount_factor,
        code_params, network_specs, opt_params, mcts_params, gate_set
    )
    rng = MersenneTwister(rand(UInt))

    solver = MCTSSolver(
        n_iterations = mdp.mcts_params.n_iterations, # if no solution is found, this should be increased
        depth = mdp.mcts_params.depth, # similar to above
        exploration_constant = mdp.mcts_params.exploration_constant, # should be increased if search space is not explored well
        rng = rng,
        reuse_tree = true,
        enable_tree_vis = false,
        #estimate_value = 0.0,
    )
    # solver = DPWSolver(n_iterations=n_iter, depth=depth, exploration_constant=ec, alpha_state=1/8, tree_in_info=true)


    policy = solve(solver, mdp)

    s = initialstate()# CircuitState(Gate[])
    MCTS_gate_counts = (typemax(Float16), typemax(Float16), typemax(Float16))
    MCTS_fidelity = typemin(Float16)
    MCTS_circuit = s
    
    steps = 0
    while steps < mdp.mcts_params.max_steps 
        a = action(policy, s)
        if mdp.code_params isa CodeParameters
            s, r, gate_counts = POMDPs.gen(mdp, s, a, Random.GLOBAL_RNG)
            println("Applied $(a), reward=$r for fidelity=$(s.fidelity) and DQC_depth=$gate_counts, MC tree depth=$(length(s.gates))")
            MCTS_circuit = s
            MCTS_gate_counts = gate_counts
            MCTS_fidelity = s.fidelity
            if POMDPs.isterminal(mdp, s)#, fidelity=fidelity) 
                println("Terminal condition reached after $(steps+1) steps. Final fidelity: $(s.fidelity).")
                break
            end
            
        # elseif mdp.code_params isa CodeParametersLogical
        #     s, r, q, gate_counts = POMDPs.gen(mdp, s, a, Random.GLOBAL_RNG)
        #     println("Applied $(a), reward=$r for quality=$q and DQC_depth=$gate_counts, MC tree depth=$(length(s.gates))")

        #     if POMDPs.isterminal(mdp, s) 
        #         println("Terminal condition reached after $(steps+1) steps. Final reward: $(r).")
        #         MCTS_gate_counts = gate_counts
        #         MCTS_fidelity = r
        #         MCTS_circuit = s
        #         break
        #     end
        end
        
        steps += 1
    end
    
    #### Verification

    MCTS_exec_circuit, _,_,_= construct_executable_circuit(MCTS_circuit.gates, mdp.gate_set, mdp.network_specs, telegate_overhead = true)

    println("\n Optimised circuit length (DQC setting):: Single-qubit gates: $(MCTS_gate_counts[1]), Two-qubit gates: $(MCTS_gate_counts[2]), Telegates: $(MCTS_gate_counts[3])") #  vs. $baseline_exec_circuit_size in baseline")
    verification_logical_state = nothing
    if mdp.code_params isa CodeParameters
        verification_logical_state = verify_success(MCTS_exec_circuit, mdp.code_params.target_state, mdp.network_specs)
        # ^NOTE: this appends a verifyop operation, but we count before so this is irrelevant
        verification_logical_state = verification_logical_state == 1.0 ? true : false
        println("\nVerification MCTS Algorithm successful (target state fidelity; only expressive (binary) in noiseless setting): $verification_logical_state")
    end

    ###########################################
    ############# DATA STORAGE ################
    ###########################################

    #a, info = action_info(policy, MCTS_circuit)
    #inchrome(D3Tree(info[:tree]))

    base_mcts_dir = joinpath(@__DIR__, "results", string(code_dirname(mdp.code_params.qec_code)), "MCTS")
    MCTS_dir = next_run_dir(base_mcts_dir)

    println("Saving results to $(MCTS_dir)")
    save_circuit_diagram(MCTS_circuit.gates, mdp.gate_set, MCTS_dir, "MCTS_raw_circuit__size_$(sum(MCTS_gate_counts)).png")
    #save_circuit_diagram(MCTS_exec_circuit, MCTS_dir, "MCTS_exec_circuit__size_$(MCTS_gate_counts).png")

    open(joinpath(MCTS_dir, "network_specs.txt"), "w") do io
        println(io, "Network Specifications")
        for fname in fieldnames(Types.NetworkSpecifications)
            println(io, fname, " = ", repr(getfield(network_specs, fname)))
        end
    end
    
    if code_params isa CodeParameters
        open(joinpath(MCTS_dir, "code_params.txt"), "w") do io
            println(io, "Code parameters")
            for fname in fieldnames(Types.CodeParameters)
                println(io, fname, " = ", repr(getfield(code_params, fname)))
            end
        end
    # elseif code_params isa CodeParametersLogical
    #     open(joinpath(MCTS_dir, "code_params.txt"), "w") do io
    #         println(io, "Code parameters")
    #         for fname in fieldnames(Types.CodeParametersLogical)
    #             println(io, fname, " = ", repr(getfield(code_params, fname)))
    #         end
    end
            

    open(joinpath(MCTS_dir, "MCTS_params.txt"), "w") do io
        println(io, "MCTS parameters")
        for fname in fieldnames(Types.MCTSParameters)
            println(io, fname, " = ", repr(getfield(mcts_params, fname)))
        end
    end

    open(joinpath(MCTS_dir, "MCTS_raw_circuit.txt"), "w") do io
        println(io, "# Raw gate sequence of size $(sum(MCTS_gate_counts))")
        for (i, g) in enumerate(MCTS_circuit.gates)
            println(io, i, "\t", repr(g))
        end
    end

    open(joinpath(MCTS_dir, "MCTS_exec_circuit.txt"), "w") do io
        println(io, "# Executable (DQC) circuit operations of size:: $MCTS_gate_counts")
        for (i, op) in enumerate(MCTS_exec_circuit)
            println(io, i, "\t", repr(op))
        end
    end

    open(joinpath(MCTS_dir, "summary.txt"), "w") do io
        println(io, "# Encoding successful: $verification_logical_state")
        println(io, "# Raw gate sequence of size $(sum(MCTS_gate_counts))")
        println(io, "# Executable (DQC) circuit operations of size $(MCTS_gate_counts) (excl. SWAPS)")
    end
    
    println("MCTS OUT GATES: $(MCTS_circuit.gates) with fidelity $(MCTS_circuit.fidelity)")
    return MCTS_circuit.gates
end


end