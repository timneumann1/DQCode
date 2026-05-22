module DQCode

include("types.jl")
include("trivariate_bicycle_code.jl")
include("symplectic_double_code.jl") 
include("helper.jl")
include("experiment_config.jl")
include("baseline_encoding.jl")
include("encoding_gott.jl")
include("genetic.jl")
include("mcts.jl")
include("dqc_simulator.jl")

using .Types
using .TrivariateBicycleCode
using .SymplecticDoubleCode
using .Helper: tableau_to_bitmatrix, data_qubit_partitioning, perm_to_transpositions, create_lookup_array, verify_success, execute_circuit, code_dirname, save_txt, save_circuit_diagram, qc_circuit_to_qasm
using .ExperimentConfig: experiment_configurations#distributed_qec_code, type_two_register_sizes, opt_params, genetic_params, mcts_params, gate_set#, noise_model, n_shots
using .EncodingGott: encoding_gott
using .Genetic: genetic_search
using .BaselineEncoding: run_qiskit_baseline, run_mqt_baseline
using .MonteCarloTreeSearch: monte_carlo_tree_search
using .DQCodeSimulator: dqc_ft_encoding_simulation

using QECCore
using QuantumClifford
using QuantumClifford: MixedDestabilizer, Stabilizer, Tableau, stabilizerview, logicalxview, logicalzview, canonicalize_rref!, tab, AbstractOperation
using QuantumClifford.ECC: DistanceMIPAlgorithm
using HiGHS
using JuMP
using DataFrames, CSV
using Logging
using Serialization
using Random

const MQT_PATH = "/Users/tim/Tim/projects/mqt/qecc/"

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
    num_X_checks = count(any( tableau_to_bitmatrix(tab(target_canon)) .== 1, dims = 2))# count how many rows contain X stabiliser (clean for CSS codes) #count(i -> tableau_to_bitmatrix(tab(target_canon))[i,i]==1, 1:size(target_canon,1))
    target_canon_rref = canonicalize_rref!(copy(target_state))
    target_tableau = tab(target_canon_rref[1])
    target_bit_matrix = tableau_to_bitmatrix(target_tableau)
    logical_Zs  = logicalzview(code)
    code_distance = _code_distance(qec_code)
    @info "Setup complete: $(qec_code)[$(code_n(qec_code)), $(code_k(qec_code)), $code_distance]]-code"
        
    code_params = CodeParameters(
        qec_code,
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

function _code_distance(qec_code)
    
    try
        return QECCore.distance(qec_code)
    catch err
        if err isa MethodError
            return QECCore.distance(qec_code, DistanceMIPAlgorithm(solver=HiGHS))
        else
            @warn "Code distance computation failed: setting d = 0" err
            rethrow()
        end
    end
end

# ----------------------------------------------
# ------------ BASELINE ENCODING ---------------
# ----------------------------------------------


function baseline_encoding_qiskit(exp_label::String)
    # this function orchestrates a baseline run, gathering code and networking parameters, and then
    # saves the results to data/{QEC_code}>{Network Architecture}>{qiskit}
    Random.seed!(42) 
    configs = experiment_configurations()

    if haskey(configs, exp_label)
        cfg = configs[exp_label]
        folder = joinpath(@__DIR__,"..","data", string(code_dirname(cfg.code)), string(cfg.qpu_sizes))

        if !isfile(joinpath(folder, "network_specs.jls")) || !isfile(joinpath(folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run create_code_network_data($exp_label).")
        end

        network_specs = deserialize( joinpath(folder, "network_specs.jls"))
        code_params = deserialize( joinpath(folder, "code_params.jls"))

        encoding_circ, gcounts, verification_logical_state = run_qiskit_baseline(code_params, network_specs)

        # Save data to data/ folder
        dir = joinpath(folder, "qiskit_encoding")
        mkpath(dir)
        @info "Saving results to $dir ..."
        serialize( joinpath(dir, "qiskit_encoding_circuit.jls"), encoding_circ )
        save_circuit_diagram(encoding_circ, dir, "qiskit_encoding_circuit.png")
        df = DataFrame(method = ["qiskit_encoding"], verified = [verification_logical_state], gate_counts = [gcounts])
        CSV.write(joinpath(dir, "qiskit_encoding_stats.csv"), df)
    else
        error("The configuration label $exp_label was not found. Please add the respective data to the configuration file first.")
    end
    
    return dir

end

function baseline_encoding_mqt(exp_label::String, mqt_path:: String, prep_method::String)
    # this function orchestrates a baseline run, gathering code and networking parameters, and then
    # saves the results to data/{QEC_code}>{Network Architecture}>{qiskit}
    Random.seed!(42) 
    configs = experiment_configurations()

    if haskey(configs, exp_label)
        cfg = configs[exp_label]
        folder = joinpath(@__DIR__,"..","data", string(code_dirname(cfg.code)), string(cfg.qpu_sizes))

        if !isfile(joinpath(folder, "network_specs.jls")) || !isfile(joinpath(folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run create_code_network_data($exp_label).")
        end

        network_specs = deserialize( joinpath(folder, "network_specs.jls"))
        code_params = deserialize( joinpath(folder, "code_params.jls"))

        encoding_circ, gcounts, verification_logical_state = run_mqt_baseline(code_params, network_specs, mqt_path, prep_method)

        # Save data to data/ folder
        dir = joinpath(folder, "mqt_encoding")
        mkpath(dir)
        @info "Saving results to $dir ..."
        serialize( joinpath(dir, "mqt_encoding_circuit.jls"), encoding_circ )
        save_circuit_diagram(encoding_circ, dir, "mqt_encoding_circuit.png")
        df = DataFrame(method = ["mqt_encoding"], verified = [verification_logical_state], gate_counts = [gcounts])
        CSV.write(joinpath(dir, "mqt_encoding_stats.csv"), df)

    else
        error("The configuration label $exp_label was not found. Please add the respective data to the configuration file first.")
    end
    
    return dir

end


# -----------------------------------------
# ------- ENCODING OPTIMISATION -----------
# -----------------------------------------

# ---------- Genetic Algorithm ------------

function circuit_search_gott_ga(exp_label::String)
    # this function orchestrates an experiment, gathering code and networking parameters, and then
    # initialising the Gottesman encoding warmstart Genetic Algorithm 
    # It saves the results in {QEC_code}>{Network Architecture}>{Optimiser}
    Random.seed!(42) 
    configs = experiment_configurations()

    if haskey(configs, exp_label)
        cfg = configs[exp_label]
        folder = joinpath(@__DIR__,"..","data", string(code_dirname(cfg.code)), string(cfg.qpu_sizes))

        if !isfile(joinpath(folder, "network_specs.jls")) || !isfile(joinpath(folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run create_code_network_data($exp_label).")
        end

        network_specs = deserialize( joinpath(folder, "network_specs.jls"))
        code_params = deserialize( joinpath(folder, "code_params.jls"))

        # Gottesman standard encoding + DQC compilation
        gottesman_encoding_circuit, dqc_compiled_encoding_circuit, verification_logical_state, verification_logical_state_compiled, gate_counts, gate_counts_compiled = encoding_gott(code_params, network_specs, zeros(code_params.k))

        # Warmstart GA with DQC-compiled circuit

        GA_encoding_circuit, verification_ga, gate_counts_ga, fitness_evolution, fidelity_evolution, gate_count_evolution = genetic_search(code_params, network_specs, cfg.genetic_params, dqc_compiled_encoding_circuit)
        
        # ----- Data Storage ----------
        dir = joinpath(folder, "warmstart_ga")
        mkpath(dir)

        serialize( joinpath(dir, "gott_encoding_circuit.jls"), gottesman_encoding_circuit )
        serialize( joinpath(dir, "gott_circuit_dqc_compiled.jls"), dqc_compiled_encoding_circuit )
        serialize(joinpath(dir, "GA_circuit.jls"), GA_encoding_circuit)
        save_txt(dir, "genetic_algorithm_parameters.txt", cfg.genetic_params)

        save_circuit_diagram(gottesman_encoding_circuit, dir, "gott_encoding_circuit.png")
        save_circuit_diagram(dqc_compiled_encoding_circuit, dir, "gott_circuit_dqc_compiled.png")
        save_circuit_diagram(GA_encoding_circuit, dir, "GA_circuit.png")


        df = DataFrame(method = ["Gottesman", "dqc_compiled", "GA"], verfied = [verification_logical_state, verification_logical_state_compiled, verification_ga],gate_counts = [gate_counts, gate_counts_compiled, gate_counts_ga] )
        CSV.write(joinpath(dir, "warm_start_ga_stats.csv"), df)

        df_ga = DataFrame(
            fitness_evolution   = fitness_evolution,
            fidelity_evolution  = fidelity_evolution,
            single_count = [gc[1] for gc in gate_count_evolution],
            two_qubit_count = [gc[2] for gc in gate_count_evolution],
            telegate_count = [gc[3] for gc in gate_count_evolution]
        )

        CSV.write(joinpath(dir, "genetic_evolution.csv"), df_ga)

        return dir
    else
        error("The configuration label $exp_label was not found. Please add the respective data to the configuration file first.")
    end
    

end

# ---------- Monte Carlo Tree Search ----------

function circuit_search_mcts(exp_label::String)

    Random.seed!(42) 
    configs = experiment_configurations()
    
    if haskey(configs, exp_label)
        cfg = configs[exp_label]
        folder = joinpath(@__DIR__,"..","data", string(code_dirname(cfg.code)), string(cfg.qpu_sizes))

        if !isfile(joinpath(folder, "network_specs.jls")) || !isfile(joinpath(folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run create_code_network_data($exp_label).")
        end
        network_specs = deserialize( joinpath(folder, "network_specs.jls"))
        code_params = deserialize( joinpath(folder, "code_params.jls"))

    
        MCTS_circuit, verification_MCTS_logical_state, MCTS_gate_counts, fidelity_evolution, gate_count_evolution, reward_evolution = monte_carlo_tree_search(code_params, network_specs, cfg.mcts_params)
       
        # ----- Data Storage ----------
        dir = joinpath(folder, "mcts")
        mkpath(dir)

        serialize( joinpath(dir, "MCTS_circuit.jls"), MCTS_circuit )
        save_circuit_diagram(MCTS_circuit, dir, "MCTS_circuit.png")
        save_txt(dir, "mcts_parameters.txt", cfg.mcts_params)

        df_evol = DataFrame(
            fidelity_evolution = fidelity_evolution, 
            gate_count_evolution = gate_count_evolution,
            reward_evolution = reward_evolution
        )

        CSV.write(joinpath(dir, "mcts_evolution.csv"), df_evol)

        df = DataFrame(method = ["MCTS"], verified = [verification_MCTS_logical_state], gate_counts = [MCTS_gate_counts])
        CSV.write(joinpath(dir, "mcts_stats.csv"), df)

        return dir
    else
        error("The serialized specification and parameter files for this experiment are missing. Please run create_code_network_data() first.")
    end

end


# -----------------------------------------
# ---------- DQC Execution ----------------
# -----------------------------------------

function dqc_simulation(exp_label::String, mqt_path::String, circuit_path::String, num_samples, ps, p_bells, telegate_idle_depth, method::String)
    Random.seed!(42) 
    configs = experiment_configurations()
    
    if haskey(configs, exp_label)
        cfg = configs[exp_label]
        folder = joinpath(@__DIR__,"..","data", string(code_dirname(cfg.code)), string(cfg.qpu_sizes))

        if !isfile(joinpath(folder, "network_specs.jls")) || !isfile(joinpath(folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run create_code_network_data($exp_label).")
        end
        network_specs = deserialize( joinpath(folder, "network_specs.jls"))
        code_params = deserialize( joinpath(folder, "code_params.jls"))

        circ_path = joinpath(folder,circuit_path) 
        encoding_circuit = deserialize(circ_path)

        if method == "none"
            
            data, data_circuit, DQC_circuit_noiseless = dqc_non_ft_encoding_simulation(code_params, network_specs, encoding_circuit)

            # ----- Data Storage ----------
            dir = joinpath(folder, "simulation_non_FT")
            mkpath(dir)
            df = DataFrame(data)
            CSV.write(joinpath(dir, "dqc_sim_data.csv"), df)
            serialize( joinpath(dir, "data_circuit.jls"), data_circuit )
            serialize( joinpath(dir, "DQC_circuit_zero_noise.jls"), DQC_circuit_noiseless )
            save_circuit_diagram(DQC_circuit_noiseless, dir, "DQC_circuit_zero_noise.png")

        else

            data, data_circuit, quantum_clifford_verification_circ, DQC_circuit_noiseless, num_ancillas, num_x_anc, num_z_anc, ancilla_map = dqc_ft_encoding_simulation(num_samples, ps, p_bells, telegate_idle_depth, code_params, network_specs, mqt_path, encoding_circuit, method)
            
            # ----- Data Storage ----------
            dir = joinpath(folder, "simulation_FT")
            mkpath(dir)
            df = DataFrame(data)
            CSV.write(joinpath(dir, "dqc_sim_data.csv"), df)
            serialize( joinpath(dir, "data_circuit.jls"), data_circuit )
            serialize( joinpath(dir, "verification_circuit.jls"), quantum_clifford_verification_circ )
            serialize( joinpath(dir, "DQC_circuit_zero_noise.jls"), DQC_circuit_noiseless )
            save_circuit_diagram(DQC_circuit_noiseless, dir, "DQC_circuit_zero_noise.png")
            ancilla_info = (; num_ancillas, num_x_anc, num_z_anc, ancilla_map)
            save_txt(dir, "ancilla_info.txt", ancilla_info)
            
        end

        return dir
    else
        error("The configuration label $exp_label was not found. Please add the respective data to the configuration file first.")
    end
end


end