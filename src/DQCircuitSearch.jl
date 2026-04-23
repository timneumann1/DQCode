module DQCircuitSearch

include("types.jl")
include("trivariate_bicycle_code.jl")
include("helper.jl")
include("experiment_config.jl")
#include("circsim.jl")
include("baseline_encoding.jl")
#include("plots.jl")
#include("dtsimulation.jl")
include("genetic.jl")
include("mcts.jl")
include("qec_tools.jl")
include("dqc_state_prep_sim.jl")

using .Types
using .TrivariateBicycleCode
using .Helper: tableau_to_bitmatrix, data_qubit_partitioning, perm_to_transpositions, create_lookup_array, verify_success, execute_circuit
using .ExperimentConfig: experiment_configurations#distributed_qec_code, type_two_register_sizes, opt_params, genetic_params, mcts_params, gate_set#, noise_model, n_shots
using .BaselineEncoding: baseline_encoding
using .Genetic: genetic_search
using .MonteCarloTreeSearch: monte_carlo_tree_search

using Logging
using Serialization
using QECCore
using QECCore: distance
using QuantumClifford
using QuantumClifford: MixedDestabilizer, Stabilizer, Tableau, stabilizerview, logicalxview, logicalzview, canonicalize_rref!, tab, AbstractOperation
using QuantumClifford.ECC: DistanceMIPAlgorithm
using HiGHS
using JuMP

# -------------------------------------
# ------------ SETUP ------------------
# -------------------------------------

function create_code_network_data()
    
    for (exp_label, cfg) in experiment_configurations()
        @info "Writing experiment configuration for $exp_label configuration to $(cfg.folder)"
        network_specs = network_setup(cfg.code, cfg.qpu_sizes)
        code_params = code_setup(cfg.code)
        mkpath(cfg.folder)
        serialize( joinpath(cfg.folder, "network_specs.jls"), network_specs )
        serialize( joinpath(cfg.folder, "code_params.jls"), code_params )
        open(joinpath(cfg.folder, "network_specs.txt"), "w") do io
            for fn in fieldnames(typeof(network_specs))
                println(io, fn, " = ", repr(getfield(network_specs, fn)))
            end
        end
        open(joinpath(cfg.folder, "code_params.txt"), "w") do io
            for fn in fieldnames(typeof(code_params))
                println(io, fn, " = ", repr(getfield(code_params, fn)))
            end
        end
            end
end

# ---------- Quantum Network ----------

function network_setup(qec_code, register_sizes)

    @assert sum(register_sizes) == code_n(qec_code) "$(code_n(qec_code))"
    #println("Stabilizers of $qec_code: $(Stabilizer(qec_code))")
    permutation = data_qubit_partitioning(register_sizes, Stabilizer(qec_code))
    # Note: Whereas for partitioning, we use the original stabiliser formalism (used in QEC cycles), 
    # for optimisation, we use the canonical form again for commensurability -- the target state is identical since 
    # it is fixed by the entire stabiliser subgroup
    inv_perm = invperm(permutation)
    mapping = perm_to_transpositions(deepcopy(permutation))    
    #num_registers = length(networking_params.register_sizes)
    register_lookup_array, data_qubits, num_data_qubits, num_comm_qubits_per_register = create_lookup_array(register_sizes)
    @assert num_data_qubits == code_n(qec_code)
    num_qubits = num_data_qubits + num_comm_qubits_per_register * length(register_sizes)
    all_qubits = collect(1:num_qubits)
    comm_qubits = setdiff(all_qubits, data_qubits)
    comm_idx = [ ( index + num_comm_qubits_per_register * (register_lookup_array[index]-1) ) for index in 1:num_data_qubits]
    @assert comm_idx == data_qubits # both arrays are capturing the same mapping
    comm_inv_perm_idx = [ (inv_perm[index] + num_comm_qubits_per_register * (register_lookup_array[inv_perm[index]]-1) ) for index in 1:num_data_qubits]

    #@info "Number of qubits is $num_qubits, number of comm. qubits per register is $num_comm_qubits_per_register"
    #println("Lookup Array: $register_lookup_array")
    #println("Data qubits: $data_qubits")

    network_specs = NetworkSpecifications(
        register_sizes, #register sizes of Type-II architecture (here: only fewn memory qubits per core), CircuitSim automatically adds comm. qubits (ancillas are only added in DTS)
        length(register_sizes), # number of registers
        permutation, # permutation for initial mapping
        mapping,     # corresponding transpositions for initial mapping
        inv_perm,    # inverse permutation (to be used in circuit execution)
        register_lookup_array,  # lookup array for core membership 
        data_qubits, # array of data qubit indices
        comm_qubits, # array of communication qubut indices
        num_data_qubits, # number of data qubits
        num_comm_qubits_per_register, # number of communication qubits per register
        num_qubits,     # total number of qubits
        comm_idx,       # array containing indices under communication qubit reindexing
        comm_inv_perm_idx, # array containing indices under communication qubit reindexing after inverse permutation
    )
    return network_specs

end

# ---------- QEC Code ----------

function code_setup(qec_code)

    code = MixedDestabilizer(qec_code)

    target_state = vcat(stabilizerview(code), logicalzview(code))
    target_canon = canonicalize!(copy(target_state))
    num_X_checks = count(any( tableau_to_bitmatrix(tab(target_canon)) .== 1, dims = 2))# count how many rows contain X stabiliser (clean for CSS codes) #count(i -> tableau_to_bitmatrix(tab(target_canon))[i,i]==1, 1:size(target_canon,1))
    @info "Number of X checks is $num_X_checks"
    target_canon_rref = canonicalize_rref!(copy(target_state))
    target_tableau = tab(target_canon_rref[1])
    target_bit_matrix = tableau_to_bitmatrix(target_tableau)
    code_distance = distance(qec_code, DistanceMIPAlgorithm(solver=HiGHS))

    #print("CODE STABILISERS LOW WEIGHT: \n$(Stabilizer(qec_code))\n")
    #println("Logical Z operators are \n$(logicalzview(code))\n")
    #println("Logical X operators are \n$(logicalxview(code))\n")
    #println("\nTarget state:$target_state\n")
    @info "QEC code: $(qec_code)[$(code_n(qec_code)), $(code_k(qec_code)), $code_distance]]-code"


    stabilizer_generators = [copy(p) for p in stabilizerview(code)]
    logical_Zs  = logicalzview(code)
    stabilizer_group = [stabilizer_generators[1] * stabilizer_generators[1]] # start with the identity element of the stabilizer group

    for gen in stabilizer_generators
        new_elements = [s * gen for s in stabilizer_group]
        append!(stabilizer_group, new_elements)
    end

    @assert length(Set(stabilizer_group)) == 2^(code_n(qec_code)-code_k(qec_code)) "Stabilizer group does not have the correct size of 2^|stabilizer_generators|"
    
    code_params = CodeParameters(
        qec_code,
        Stabilizer(qec_code),
        stabilizer_group,
        num_X_checks,
        logical_Zs,
        target_state,
        target_bit_matrix,
        code_n(qec_code),
        code_k(qec_code),
        code_distance
    )
    return code_params

end

# -------------------------------------
# ------------ OPTIMISATION -----------
# -------------------------------------

function run_circuit_search()
    # this function orchestrates an experiment, gathering code and networking parameters, and then
    # initialising the baseline, GA, (multiple) MCTS, and baseline>GA, MCTS>GA runs 
    # It saves the results in {QEC_code}>{Network Architecture}>{Optimiser}
    
    for (exp_label, cfg) in experiment_configurations()

        @info "Loading experiment configuration for $exp_label configuration from $(cfg.folder)" 

        if !isfile(joinpath(cfg.folder, "network_specs.jls")) || !isfile(joinpath(cfg.folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run create_code_network_data() first.")
        end

        network_specs = deserialize( joinpath(cfg.folder, "network_specs.jls"))
        code_params = deserialize( joinpath(cfg.folder, "code_params.jls"))

        baseline_gates = baseline_encoding(code_params, network_specs, cfg.folder)

        best_MCTS_gates = Vector{Gate}() 
        best_MCTS_telegate_count = typemax(Int)
        for _ in 1:3
            MCTS_gates, verification_logical_state, telegate_count = monte_carlo_tree_search(code_params, network_specs, cfg.opt_params, cfg.mcts_params, cfg.folder)
            if verification_logical_state && telegate_count < best_MCTS_telegate_count
                best_MCTS_gates = MCTS_gates
                best_MCTS_telegate_count = telegate_count
            end
        end

        genetic_search(code_params, network_specs, cfg.opt_params, cfg.genetic_params, cfg.folder)

        genetic_search(code_params, network_specs, cfg.opt_params, cfg.genetic_params, cfg.folder, warm_start = true, warm_start_gates = baseline_gates, label = "Baseline")

        if !isempty(best_MCTS_gates)
            genetic_search(code_params, network_specs, cfg.opt_params, cfg.genetic_params, cfg.folder, warm_start = true, warm_start_gates = best_MCTS_gates, label = "MCTS")
        end

    end
end

# ---------- Optimiser Runs ----------

function run_baseline_encoding(exp_label::String)
    configs = experiment_configurations()
    if haskey(configs, exp_label)
        cfg = configs[exp_label]
        if !isfile(joinpath(cfg.folder, "network_specs.jls")) || !isfile(joinpath(cfg.folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run create_code_network_data() first.")
        end
        network_specs = deserialize( joinpath(cfg.folder, "network_specs.jls"))
        code_params = deserialize( joinpath(cfg.folder, "code_params.jls"))
        baseline_encoding(code_params, network_specs, cfg.folder)
        return 42
    else
        error("This configuration label $exp_label was not found. Please add the respective data to the configuration file")
    end
end


function run_MCTS(exp_label::String)
    configs = experiment_configurations()
    if haskey(configs, exp_label)
        cfg = configs[exp_label]
        if !isfile(joinpath(cfg.folder, "network_specs.jls")) || !isfile(joinpath(cfg.folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run create_code_network_data() first.")
        end
        network_specs = deserialize( joinpath(cfg.folder, "network_specs.jls"))
        code_params = deserialize( joinpath(cfg.folder, "code_params.jls"))
        monte_carlo_tree_search(code_params, network_specs, cfg.opt_params, cfg.mcts_params, cfg.folder)
        return 42
    else
        error("This configuration label $exp_label was not found. Please add the respective data to the configuration file")
    end 
end

function run_genetic_search(exp_label::String)
    configs = experiment_configurations()
    if haskey(configs, exp_label)
        cfg = configs[exp_label]
        if !isfile(joinpath(cfg.folder, "network_specs.jls")) || !isfile(joinpath(cfg.folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run create_code_network_data() first.")
        end
        network_specs = deserialize( joinpath(cfg.folder, "network_specs.jls"))
        code_params = deserialize( joinpath(cfg.folder, "code_params.jls"))
        genetic_search(code_params, network_specs, cfg.opt_params, cfg.genetic_params, cfg.folder)

        ########################

        # baseline_gates = baseline_encoding(code_params, network_specs, cfg.folder)
        # MCTS_gates, verification_logical_state, telegate_count = monte_carlo_tree_search(code_params, network_specs, cfg.opt_params, cfg.mcts_params, cfg.folder)
        # genetic_search(code_params, network_specs, cfg.opt_params, cfg.genetic_params, cfg.folder, warm_start = true, warm_start_gates = baseline_gates, label = "Baseline")

        ##############


        return 42
    else
        error("This configuration label $exp_label was not found. Please add the respective data to the configuration file")
    end
end

# using .Circuit_Plots: plot_gate_teleportation
# export plot_gate_teleportation


end