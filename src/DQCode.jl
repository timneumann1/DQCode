module DQCode

include("types.jl")
include("trivariate_bicycle_code.jl")
include("helper.jl")
include("experiment_config.jl")
#include("circsim.jl")
include("baseline_encoding.jl")
include("encoding_gott.jl")

#include("plots.jl")
#include("dtsimulation.jl")

include("genetic.jl")
include("mcts.jl")

#include("qec_tools.jl")
#include("dqc_state_prep_sim.jl")

using .Types
using .TrivariateBicycleCode
using .Helper: tableau_to_bitmatrix, data_qubit_partitioning, perm_to_transpositions, create_lookup_array, verify_success, execute_circuit, code_dirname, save_txt
using .ExperimentConfig: experiment_configurations#distributed_qec_code, type_two_register_sizes, opt_params, genetic_params, mcts_params, gate_set#, noise_model, n_shots
using .EncodingGott: encoding_gott
using .Genetic: genetic_search
using .BaselineEncoding: run_qiskit_baseline, run_mqt_baseline
using .MonteCarloTreeSearch: monte_carlo_tree_search

using Logging
using Serialization
using Random

# setting the seed yields reproducibility among runs of the entire code, but within one 10 iteration execution of the MCTS solver, e.g.
# we will obtain different results (but everytime the same different ones)
# thus, we set one 'global' seed for each experiment function we are implemting 
using QECCore
using QECCore: distance
using QuantumClifford
using QuantumClifford: MixedDestabilizer, Stabilizer, Tableau, stabilizerview, logicalxview, logicalzview, canonicalize_rref!, tab, AbstractOperation
using QuantumClifford.ECC: DistanceMIPAlgorithm
using HiGHS
using JuMP
using DataFrames, CSV


#export create_code_network_data, network_setup, code_setup

# -------------------------------------
# ------------ SETUP ------------------
# -------------------------------------

function create_code_network_data(exp_label::String)
     
    # Creates all code and network data for specified experiment from the config file
    configs = experiment_configurations()
    if haskey(configs, exp_label)
        # Retrieve configuration
        cfg = configs[exp_label]

        # Create networking specifications and code parameters
        network_specs = _network_setup(cfg.code, cfg.qpu_sizes)
        code_params = _code_setup(cfg.code)    

        # Save data to data/ folder
        folder = joinpath(@__DIR__,"..","data", string(code_dirname(cfg.code)), string(cfg.qpu_sizes))
        mkpath(folder)
        @info "Writing experiment configuration for $exp_label configuration to data/ folder"
        serialize( joinpath(folder, "network_specs.jls"), network_specs )
        serialize( joinpath(folder, "code_params.jls"), code_params )        
        save_txt(folder, "network_specs.txt", network_specs)
        save_txt(folder, "code_params.txt", code_params)

    else
        error("The configuration label $exp_label was not found. Please add the respective data to the configuration file.")
    end
end

# ---------- Quantum Network ----------

function _network_setup(qec_code, register_sizes)
    # No need to set a seed since the only prob. part is the hypergraph partitioning, which has its own seed (already defined in the config file) -> a good reason to execute all the different code-architecture pieces separately (all will be reprodcuble, indpendet of order of execution)
    @assert sum(register_sizes) == code_n(qec_code) "$(code_n(qec_code))"
    mapping = data_qubit_partitioning(register_sizes, Stabilizer(qec_code))
    # Note: Whereas for partitioning, we use the original stabiliser formalism (used in QEC cycles), 
    # for optimisation, we use the canonical form again for commensurability -- the target state is identical since 
    # it is fixed by the entire stabiliser subgroup
    # mapping contains the permutation that indicates for each position, which qubit will be there.
    mapping_transpositions = perm_to_transpositions(deepcopy(mapping))  
    inv_map = invperm(mapping)
    register_lookup_array, num_data_qubits, num_comm_qubits_per_register = create_lookup_array(register_sizes)
    @assert num_data_qubits == code_n(qec_code)
    num_comm_qubits = num_comm_qubits_per_register * length(register_sizes)
    num_data_and_comm_qubits = num_data_qubits + num_comm_qubits
    data_qubits = collect(1:num_data_qubits)
    data_and_comm_qubits = collect(1:num_data_and_comm_qubits)
    comm_qubits = setdiff(data_and_comm_qubits, data_qubits)

    network_specs = NetworkSpecifications(
        register_sizes, #register sizes of Type-II architecture (here: only fewn memory qubits per core), CircuitSim automatically adds comm. qubits (ancillas are only added in DTS)
        length(register_sizes), # number of registers
        mapping, # permutation for initial mapping
        mapping_transpositions,     # corresponding transpositions for initial mapping, to be consumed from right to left
        inv_map,    # inverse permutation (to be used in circuit execution)
        register_lookup_array,  # lookup array for core membership 
        data_qubits, # array of data qubit indices
        comm_qubits, # array of communication qubut indices
        num_data_qubits, # number of data qubits
        num_comm_qubits, # number of communication qubits
        num_comm_qubits_per_register, # number of communication qubits per register
        num_data_and_comm_qubits,     # total number of qubits
        #dqc_data_idx,       # array containing data qubit indices under communication qubit reindexing, i.e., in the DQC setting
        #comm_inv_perm_idx, # array containing indices under communication qubit reindexing after inverse permutation
    )
    return network_specs

end

# ---------- QEC Code ----------

function _code_setup(qec_code)

    code = MixedDestabilizer(qec_code)

    target_state = vcat(stabilizerview(code), logicalzview(code))
    target_canon = canonicalize!(copy(target_state))
    # can change to rref as well??
    num_X_checks = count(any( tableau_to_bitmatrix(tab(target_canon)) .== 1, dims = 2))# count how many rows contain X stabiliser (clean for CSS codes) #count(i -> tableau_to_bitmatrix(tab(target_canon))[i,i]==1, 1:size(target_canon,1))
    #@info "Number of X checks is $num_X_checks"
    target_canon_rref = canonicalize_rref!(copy(target_state))
    target_tableau = tab(target_canon_rref[1])
    target_bit_matrix = tableau_to_bitmatrix(target_tableau)

    #stabilizer_generators = [copy(p) for p in stabilizerview(code)]
    logical_Zs  = logicalzview(code)

    #if code_n(qec_code) < 20
    code_distance = 0
    try
        code_distance = distance(qec_code, DistanceMIPAlgorithm(solver=HiGHS))
        @info "Setup complete: $(qec_code)[$(code_n(qec_code)), $(code_k(qec_code)), $code_distance]]-code"
    catch err
        @info "Setup complete: $(qec_code):\n[$(code_n(qec_code)), $(code_k(qec_code)), $code_distance]]-code"
        @warn "Code distance computation failed: setting d = 0" err
    end
        
    # stabilizer_group = [stabilizer_generators[1] * stabilizer_generators[1]] # start with the identity element of the stabilizer group

    # for gen in stabilizer_generators
    #     new_elements = [s * gen for s in stabilizer_group]
    #     append!(stabilizer_group, new_elements)
    # end

    # @assert length(Set(stabilizer_group)) == 2^(code_n(qec_code)-code_k(qec_code)) "Stabilizer group does not have the correct size of 2^|stabilizer_generators|"
    

    #print("CODE STABILISERS LOW WEIGHT: \n$(Stabilizer(qec_code))\n")
    #println("Logical Z operators are \n$(logicalzview(code))\n")
    #println("Logical X operators are \n$(logicalxview(code))\n")
    #println("\nTarget state:$target_state\n")
    
    code_params = CodeParameters(
        qec_code,
        #stabilizerview(code),# Stabilizer(qec_code),
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
# ------------ BASELINE ---------------
# -------------------------------------


function baseline_encoding_qiskit(exp_label::String)
    # this function orchestrates a baseline run, gathering code and networking parameters, and then
    # It saves the results in {QEC_code}>{Network Architecture}>{qiskit}
    Random.seed!(42) 
    configs = experiment_configurations()

    if haskey(configs, exp_label)
        cfg = configs[exp_label]
        @info "Loading experiment configuration for $exp_label configuration from $(cfg.folder)" 

        if !isfile(joinpath(cfg.folder, "network_specs.jls")) || !isfile(joinpath(cfg.folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run create_code_network_data($exp_label).")
        end

        network_specs = deserialize( joinpath(cfg.folder, "network_specs.jls"))
        code_params = deserialize( joinpath(cfg.folder, "code_params.jls"))

        # gottesman, +reichardt, + genetic

        qiskit_encoding_circ, gcounts = run_qiskit_baseline(code_params, network_specs, cfg.folder)

        ## apply GENETIC SEARCH ON COMPILED VERSION

        #verification_ga, gate_counts_ga = genetic_search(code_params, network_specs, cfg.genetic_params, dqc_compiled_encoding_circuit, cfg.folder)#, label = "DQC_Compiled_Gottesman")
        
        #df = DataFrame(method = ["Gottesman", "dqc_compiled", "GA"], verfied = [verification_logical_state, verification_logical_state_compiled, verification_ga],gate_counts = [gate_counts, gate_counts_compiled, gate_counts_ga] )
       # CSV.write(joinpath(cfg.folder, "gottesman_stats.csv"), df)
    else
        error("The configuration label $exp_label was not found. Please add the respective data to the configuration file first.")
    end
    
    return 42

end

function baseline_encoding_mqt(exp_label::String)
    # this function orchestrates a baseline run, gathering code and networking parameters, and then
    # It saves the results in {QEC_code}>{Network Architecture}>{qiskit}
    Random.seed!(42) 
    configs = experiment_configurations()

    if haskey(configs, exp_label)
        cfg = configs[exp_label]
        @info "Loading experiment configuration for $exp_label configuration from $(cfg.folder)" 

        if !isfile(joinpath(cfg.folder, "network_specs.jls")) || !isfile(joinpath(cfg.folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run create_code_network_data($exp_label).")
        end

        network_specs = deserialize( joinpath(cfg.folder, "network_specs.jls"))
        code_params = deserialize( joinpath(cfg.folder, "code_params.jls"))

        # gottesman, +reichardt, + genetic

        mqt_encoding_circ, gcounts = run_mqt_baseline(code_params, network_specs, cfg.folder)

        ## apply GENETIC SEARCH ON COMPILED VERSION

        #verification_ga, gate_counts_ga = genetic_search(code_params, network_specs, cfg.genetic_params, dqc_compiled_encoding_circuit, cfg.folder)#, label = "DQC_Compiled_Gottesman")
        
        #df = DataFrame(method = ["Gottesman", "dqc_compiled", "GA"], verfied = [verification_logical_state, verification_logical_state_compiled, verification_ga],gate_counts = [gate_counts, gate_counts_compiled, gate_counts_ga] )
       # CSV.write(joinpath(cfg.folder, "gottesman_stats.csv"), df)
    else
        error("The configuration label $exp_label was not found. Please add the respective data to the configuration file first.")
    end
    
    return 42

end




# -------------------------------------
# ------------ OPTIMISATION -----------
# -------------------------------------
# function run_gottesman_encoding(exp_label::String)
#     configs = experiment_configurations()
#     if haskey(configs, exp_label)
#         cfg = configs[exp_label]
#         if !isfile(joinpath(cfg.folder, "network_specs.jls")) || !isfile(joinpath(cfg.folder, "code_params.jls"))
#             error("The serialized specification and parameter files for this experiment are missing. Please run create_code_network_data() first.")
#         end
#         network_specs = deserialize( joinpath(cfg.folder, "network_specs.jls"))
#         code_params = deserialize( joinpath(cfg.folder, "code_params.jls"))
#         encoding_gott(code_params, network_specs, cfg.folder)
#         return 42
#     else
#         error("This configuration label $exp_label was not found. Please add the respective data to the configuration file")
#     end
# end


function circuit_search_gott(exp_label::String)
    # this function orchestrates an experiment, gathering code and networking parameters, and then
    # initialising the baseline, GA, (multiple) MCTS, and baseline>GA, MCTS>GA runs 
    # It saves the results in {QEC_code}>{Network Architecture}>{Optimiser}
    Random.seed!(42) 
    configs = experiment_configurations()

    if haskey(configs, exp_label)
        cfg = configs[exp_label]
        @info "Loading experiment configuration for $exp_label configuration from $(cfg.folder)" 

        if !isfile(joinpath(cfg.folder, "network_specs.jls")) || !isfile(joinpath(cfg.folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run create_code_network_data($exp_label).")
        end

        network_specs = deserialize( joinpath(cfg.folder, "network_specs.jls"))
        code_params = deserialize( joinpath(cfg.folder, "code_params.jls"))

        # gottesman, +reichardt, + genetic

        gottesman_encoding_circuit, dqc_compiled_encoding_circuit, verification_logical_state, verification_logical_state_compiled, gate_counts, gate_counts_compiled = encoding_gott(code_params, network_specs, zeros(code_params.k), cfg.folder)

        # apply GENETIC SEARCH ON COMPILED VERSION

        verification_ga, gate_counts_ga = genetic_search(code_params, network_specs, cfg.genetic_params, dqc_compiled_encoding_circuit, cfg.folder)#, label = "DQC_Compiled_Gottesman")
        
        df = DataFrame(method = ["Gottesman", "dqc_compiled", "GA"], verfied = [verification_logical_state, verification_logical_state_compiled, verification_ga],gate_counts = [gate_counts, gate_counts_compiled, gate_counts_ga] )
        CSV.write(joinpath(cfg.folder, "gottesman_stats.csv"), df)
    else
        error("The configuration label $exp_label was not found. Please add the respective data to the configuration file first.")
    end
    
    return 42

end


function circuit_search_MCTS(exp_label::String)#, num_MCTS_runs::Int)

    Random.seed!(42) 
    configs = experiment_configurations()
    
    if haskey(configs, exp_label)
        cfg = configs[exp_label]
        
        if !isfile(joinpath(cfg.folder, "network_specs.jls")) || !isfile(joinpath(cfg.folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run create_code_network_data($exp_label).")
        end
        #println(cfg.folder)
        network_specs = deserialize( joinpath(cfg.folder, "network_specs.jls"))
        code_params = deserialize( joinpath(cfg.folder, "code_params.jls"))

        #sweep over depth, num iterations and exploration constant
        # for each combination (5^3), perform 10 runs, storing the results for saving
        #depths = [2, 3, 4, 5]# 8, 16]
        #num_iterations = [Int(1e4), Int(1e5), Int(1e6)]
        #exploration_constants = [3.5, 5.0, 10.0]

        #results_df = DataFrame(#depth = Vector{Int64}(), exploration_constants = Vector{Float64}(),
                               #run = Vector{Int64}(),
        #                       gate_counts = Vector{Vector{Int64}}(), verified = Vector{Bool}() )

        #for d in depths, c in exploration_constants
        folder_exp = joinpath(cfg.folder, "MCTS")#, string(d, "_", c)) 
        mkpath(folder_exp)
            #cfg.mcts_params.depth = d
            #cfg.mcts_params.num_iterations = niter
            #cfg.mcts_params.exploration_constant = c

            #MCTS_stats_per_comb = Vector{Tuple{Vector{Int64}, Bool}}()
            #for run in 1:num_MCTS_runs
                #Setting a random seed means we will get the same results when running this again, but we can (and will) still get differnt resuts for each of the num_MCTS_runs runs
        verification_MCTS_logical_state, MCTS_gate_counts = monte_carlo_tree_search(code_params, network_specs, cfg.mcts_params, folder_exp)
        #push!(results_df, (d,c, MCTS_gate_counts, verification_logical_state))
                #push!(MCTS_stats_per_comb, (MCTS_gate_counts, verification_logical_state) )
            #end

        df = DataFrame(method = ["MCTS"], verified = [verification_MCTS_logical_state], gate_counts = [MCTS_gate_counts])
        #df = DataFrame(run = 1:num_MCTS_runs, gate_counts = [stats[1] for stats in MCTS_stats_per_comb], verified = [stats[2] for stats in MCTS_stats_per_comb] )
        #CSV.write(joinpath(folder_exp, "MCTS_$(d)_$(niter)_$(c)_stats.csv"), df)
        #end
        #mkpath(joinpath(cfg.folder, "MCTS"))
        CSV.write(joinpath(cfg.folder, "MCTS_stats.csv"), df)
        return 42
    else
        error("The serialized specification and parameter files for this experiment are missing. Please run create_code_network_data() first.")
    end

end
        # best_MCTS_gates = Vector{Gate}() 
        # best_MCTS_telegate_count = typemax(Int)
        # best_MCTS_dir = ""

        # # LOG INTERMEDITATE GATES
        # # Sweep parameters and do statistics (probably better outsource) DONT NEED THE BEST ONE ANYMORE!
        # #three paramters: exploration constant, num_iterations, w
        # for _ in 1:num_MCTS_runs
        #     MCTS_gates, verification_logical_state, MCTS_gate_counts, MCTS_dir = monte_carlo_tree_search(code_params, network_specs, cfg.mcts_params, cfg.folder)
        #     push!(MCTS_gate_counts_per_run, MCTS_gate_counts )
        #     if isempty(best_MCTS_gates)
        #         best_MCTS_gates = MCTS_gates
        #         best_MCTS_telegate_count = MCTS_gate_counts[3]
        #         best_MCTS_dir = MCTS_dir
        #     elseif verification_logical_state && telegate_count < best_MCTS_telegate_count
        #         best_MCTS_gates = MCTS_gates
        #         best_MCTS_telegate_count = MCTS_gate_counts[3]
        #     end
        # end


          # don't do MCTS + GA
        # if !isempty(best_MCTS_gates)
        #     genetic_search(code_params, network_specs, cfg.genetic_params, cfg.folder, warm_start = true, warm_start_gates = best_MCTS_gates, MCTS_dir = best_MCTS_dir, label = "MCTS")
        # else
        #     @info "None of the $num_MCTS_runs MCTS runs was successful, please repeat the experiment."
        # end

# ---------- Optimiser Runs ----------


# function run_genetic_search(exp_label::String)
#     configs = experiment_configurations()
#     if haskey(configs, exp_label)
#         cfg = configs[exp_label]
#         if !isfile(joinpath(cfg.folder, "network_specs.jls")) || !isfile(joinpath(cfg.folder, "code_params.jls"))
#             error("The serialized specification and parameter files for this experiment are missing. Please run create_code_network_data() first.")
#         end
#         network_specs = deserialize( joinpath(cfg.folder, "network_specs.jls"))
#         code_params = deserialize( joinpath(cfg.folder, "code_params.jls"))
#         #genetic_search(code_params, network_specs, cfg.genetic_params, cfg.folder)

#         ########################

#         #baseline_gates = baseline_encoding(code_params, network_specs, cfg.folder)
#         #genetic_search(code_params, network_specs, cfg.genetic_params, cfg.folder, warm_start = true, warm_start_gates = baseline_gates, label = "Baseline")

#         MCTS_gates, verification_logical_state, telegate_count, MCTS_dir = monte_carlo_tree_search(code_params, network_specs, cfg.mcts_params, cfg.folder)
#         #MCTS_gates = deserialize( joinpath(folder, "MCTS_gates.jls"))
#         genetic_search(code_params, network_specs, cfg.genetic_params, cfg.folder, warm_start = true, warm_start_gates = MCTS_gates, MCTS_dir = MCTS_dir, label = "MCTS")

#         ##############


#         return 42
#     else
#         error("This configuration label $exp_label was not found. Please add the respective data to the configuration file")
#     end
# end


# using .Circuit_Plots: plot_gate_teleportation
# export plot_gate_teleportation


end