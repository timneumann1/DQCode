# DQCode.jl — Optimising fault-tolerant zero-state encoding on distributed quantum architectures
# Author: Tim Neumann <tneumann@tudelft.nl / tim01.neumann@gmail.com>

"""
DQCode — utilities to set up CSS codes and distributed Type-II architectures, 
optimise logical zero state encoding, and simulate encoding circuits in a DQC setting.
"""
module DQCode

export create_code_network_data, baseline_encoding_qiskit, baseline_encoding_mqt, circuit_search_gott_ga, 
    circuit_search_mcts, dqc_simulation, resource_estimation

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
include("resource_estimate.jl")

using .Types
using .TrivariateBicycleCode
using .SymplecticDoubleCode
using .Helper: tableau_to_bitmatrix, data_qubit_partitioning, perm_to_transpositions, create_lookup_array,
                verify_success, code_dirname, save_txt, save_circuit_diagram, qc_circuit_to_qasm
using .ExperimentConfig: experiment_configurations
using .EncodingGott: encoding_gott
using .Genetic: genetic_search
using .BaselineEncoding: run_qiskit_baseline, run_mqt_baseline
using .MonteCarloTreeSearch: monte_carlo_tree_search
using .DQCodeSimulator: dqc_ft_encoding_simulation, dqc_non_ft_encoding_simulation 
using .ResourceEstimation: estimate_resources_encoding_circuit, estimate_resources_measurement_based_encoding

using QECCore
using QuantumClifford: MixedDestabilizer, Stabilizer, Tableau, stabilizerview, 
                        logicalxview, logicalzview, canonicalize!, canonicalize_rref!, tab
using QuantumClifford.ECC: DistanceMIPAlgorithm
using HiGHS, JuMP
using DataFrames, CSV
using Logging
using Serialization
using Random

# -------------------------------------
# ------------ SETUP ------------------
# -------------------------------------

"""
    create_code_network_data(exp_label::String)

Create all code and network data for the specified experiment configuration and
save it to `data/<code_dirname>/<qpu_sizes>`.

### Input

- `exp_label` -- label of the experiment defined in `experiment_configurations()`

### Output

The name of the folder to which the data is saved.
"""
function create_code_network_data(exp_label::String)::String
    code_architecture_setup,_,_ = experiment_configurations()
    if haskey(code_architecture_setup, exp_label)
        cfg = code_architecture_setup[exp_label]
        network_specs = _network_setup(cfg.code, cfg.qpu_sizes)
        code_params = _code_setup(cfg.code)    
        # ----- Data Storage ----------
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
    return folder
end


"""
    _network_setup(qec_code::AbstractCSSCode, register_sizes::Vector{Int})::NetworkSpecifications

Perform hypergraph partitioning to map the data qubits of the `qec_code` onto a 
distributed Type-II architecture consisting of multiple registers, and store the 
network specifications in the corresponding type structure.

### Input

- `qec_code` -- quantum error correction code 
- `register_sizes` -- array containing the sizes (number of memory qubits) of each QPU register

### Output 

A `NetworkSpecifications` struct containing the generated network layout, including the qubit mapping,
inverse permutation, lookup arrays, and subsets of data versus communication qubits.

### Notes 

The mapping extracted from hypergraph partitioning indicates for each index of the list,
which data qubit will be mapped to this position.

Whereas for partitioning, we use the original stabiliser formalism (that will later be used in QEC cycles),
for optimisation, we use the canonical form of the tableau for commensurability, by which the target state
of the `qec_code` is uniquely identified.
"""
function _network_setup(qec_code::AbstractCSSCode, register_sizes::Vector{Int})::NetworkSpecifications
    @assert sum(register_sizes) == code_n(qec_code) "$(code_n(qec_code))"
    mapping = data_qubit_partitioning(register_sizes, Stabilizer(qec_code)) # hypergraph partitioning comes with an in-built seed
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
        register_sizes,         # register sizes of Type-II architecture
        length(register_sizes), # number of registers
        mapping,                # initial placement of qubits
        mapping_transpositions, # corresponding transpositions for mapping, to be applied from right to left
        inv_map,                # inverse mapping (to be used in circuit execution)
        register_lookup_array,  # lookup array for core membership 
        data_qubits,            # array of data qubit indices
        comm_qubits,            # array of communication qubut indices
        num_data_qubits,        # number of data qubits
        num_comm_qubits,        # number of communication qubits
        num_comm_qubits_per_register, # number of communication qubits per register
        num_data_and_comm_qubits,     # total number of data and communication qubits
    )
    return network_specs
end


"""
    _code_setup(qec_code::AbstractCSSCode)::CodeParameters

Prepare the quantum error correction code parameters, including canonical representations
and distance computation, and store them in the corresponding type structure.

### Input

- `qec_code` -- quantum error correction code

### Output 

A `CodeParameters` struct containing the code definitions, number of X checks, 
logical Z operators, canonical target state tableau, bit matrix representation, 
and the code properties (n, k, distance).

### Notes

Since we are working with CSS codes, the number of X checks equals the number of 
Hadamards we need in the initial layer, assuming that our encoding circuit follows 
a H-CNOT template.
"""
function _code_setup(qec_code::AbstractCSSCode)::CodeParameters
    code = MixedDestabilizer(qec_code)
    target_state = vcat(stabilizerview(code), logicalzview(code))
    target_canon = canonicalize!(copy(target_state))
    num_X_checks = count(any( tableau_to_bitmatrix(tab(target_canon)) .== 1, dims = 2))
    target_canon_rref = canonicalize_rref!(copy(target_state))
    target_tableau = tab(target_canon_rref[1])
    target_bit_matrix = tableau_to_bitmatrix(target_tableau)
    logical_Zs  = logicalzview(code)
    code_distance = _code_distance(qec_code)
    @info "Setup complete: $(qec_code)[$(code_n(qec_code)), $(code_k(qec_code)), $code_distance]]-code"
    code_params = CodeParameters(
        qec_code,       # QEC CSS code for logical state preparation
        num_X_checks,   # number of X checks in canonical tableau
        logical_Zs,     # canonical logical Z operator
        target_state,   # target logical zero state of the code
        target_bit_matrix,  # bit matrix encoding the tableau of the target state
        code_n(qec_code),   # number of physical qubits the code is defined on
        code_k(qec_code),   # logical space dimension of the code
        code_distance       # distance of the QEC code
    )
    return code_params
end


"""
    _code_distance(qec_code::AbstractCSSCode)::Int

Compute the code distance for the given quantum error correction code

### Input

- `qec_code` -- quantum error correction code

### Output 

The computed distance of the code.

### Notes

This function first attempts to find the distance using `QECCore.distance()`. 
If this specific method is not defined in the (QuantumClifford or custom) code
construction file (raising a `MethodError`), it falls back to computing the
exact distance via mixed-integer linear programming using
`DistanceMIPAlgorithm(solver=HiGHS)`. Any other errors are logged and rethrown.
"""
function _code_distance(qec_code::AbstractCSSCode)::Int
    try
        return QECCore.distance(qec_code)
    catch err
        if err isa MethodError
            return QECCore.distance(qec_code, DistanceMIPAlgorithm(solver=HiGHS))
        else
            @warn "code distance computation failed: setting d = 0" err
            rethrow()
        end
    end
end


# ----------------------------------------------
# ------------ BASELINE ENCODING ---------------
# ----------------------------------------------

"""
    baseline_encoding_qiskit(exp_label::String)::String

Orchestrate a baseline encoding circuit generation using the Qiskit synthesis method.
This function loads pre-existing code and network parameters, executes the baseline run, and 
saves the generated circuit and compilation statistics to `data/<code_dirname>/<qpu_sizes>/<qiskit>`.

### Input

- `exp_label` -- label of the experiment defined in `experiment_configurations()`

### Output 

The name of the folder to which the data is saved.

### Notes

This function requires that the serialized specification and parameter files 
(`network_specs.jls` and `code_params.jls`) have already been created for this 
experiment. If they are missing, you must run `create_code_network_data(exp_label)` first.
"""
function baseline_encoding_qiskit(exp_label::String)::String
    Random.seed!(42) 
    code_architecture_setup,_,_ = experiment_configurations()
    if haskey(code_architecture_setup, exp_label)
        cfg = code_architecture_setup[exp_label]
        folder = joinpath(@__DIR__,"..","data", string(code_dirname(cfg.code)), string(cfg.qpu_sizes))
        if !isfile(joinpath(folder, "network_specs.jls")) || !isfile(joinpath(folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run `qdc_setup` with $exp_label.")
        end
        network_specs = deserialize( joinpath(folder, "network_specs.jls"))
        code_params = deserialize( joinpath(folder, "code_params.jls"))
        encoding_circ, gcounts, verification_logical_state = run_qiskit_baseline(code_params, network_specs)
        # ----- Data Storage ----------
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


"""
    baseline_encoding_mqt(exp_label::String, mqt_path:: String, prep_method::String)::String

Orchestrate a baseline encoding circuit generation using the MQT encoding synthesis method.
This function loads pre-existing code and network parameters, executes the baseline run, and 
saves the generated circuit and compilation statistics to `data/<code_dirname>/<qpu_sizes>/<mqt_encoding>`.

### Input

- `exp_label` -- label of the experiment defined in `experiment_configurations()`

### Output 

The name of the folder to which the data is saved.

### Notes

This function requires that the serialized specification and parameter files 
(`network_specs.jls` and `code_params.jls`) have already been created for this 
experiment. If they are missing, you must run `create_code_network_data(exp_label)` first.
"""
function baseline_encoding_mqt(exp_label::String, mqt_path:: String, prep_method::String)::String
    Random.seed!(42) 
    code_architecture_setup,_,_ = experiment_configurations()
    if haskey(code_architecture_setup, exp_label)
        cfg = code_architecture_setup[exp_label]
        folder = joinpath(@__DIR__,"..","data", string(code_dirname(cfg.code)), string(cfg.qpu_sizes))
        if !isfile(joinpath(folder, "network_specs.jls")) || !isfile(joinpath(folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run `qdc_setup` with $exp_label.")
        end
        network_specs = deserialize( joinpath(folder, "network_specs.jls"))
        code_params = deserialize( joinpath(folder, "code_params.jls"))
        encoding_circ, gcounts, verification_logical_state = run_mqt_baseline(code_params, network_specs, mqt_path, prep_method)
        # ----- Data Storage ----------
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

"""
    circuit_search_gott_ga(exp_label::String)::String

Orchestrate an encoding circuit search using a Genetic Algorithm (GA) warm-started with a DQC-compiled
Gottesman encoding circuit. Create the standard Gottesman encoding circuit (https://arxiv.org/abs/quant-ph/9705052),
compile it for the DQC setting using a DQC-adapted overlap method based on https://arxiv.org/abs/1106.2190,
and run a genetic algorithm to optimize the circuit. Save the resulting circuits, diagrams, and evolutionary
statistics to `data/<code_dirname>/<qpu_sizes>/<warmstart_ga>`.

### Input

- `exp_label` -- label of the experiment defined in `experiment_configurations()`

### Output 

The name of the folder to which the data is saved.

### Notes

This function requires that the serialized specification and parameter files 
(`network_specs.jls` and `code_params.jls`) have already been created for this 
experiment. If they are missing, run `create_code_network_data(exp_label)` first.
"""
function circuit_search_gott_ga(exp_label::String)::String
    Random.seed!(42) 
    code_architecture_setup, genetic_parameters,_ = experiment_configurations()
    if haskey(code_architecture_setup, exp_label) && haskey(genetic_parameters, exp_label)
        cfg = code_architecture_setup[exp_label]
        genetic_params = genetic_parameters[exp_label]
        folder = joinpath(@__DIR__,"..","data", string(code_dirname(cfg.code)), string(cfg.qpu_sizes))
        if !isfile(joinpath(folder, "network_specs.jls")) || !isfile(joinpath(folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run `qdc_setup` with $exp_label.")
        end
        network_specs = deserialize( joinpath(folder, "network_specs.jls"))
        code_params = deserialize( joinpath(folder, "code_params.jls"))
        (
            gottesman_encoding_circuit, 
            dqc_compiled_encoding_circuit, 
            verification_logical_state, 
            verification_logical_state_compiled, 
            gate_counts, 
            gate_counts_compiled
        ) = encoding_gott(code_params, network_specs, zeros(Int, code_params.k))
        (
            GA_encoding_circuit,
            verification_ga,
            gate_counts_ga, 
            fitness_evolution,
            fidelity_evolution, 
            gate_count_evolution
        ) = genetic_search(code_params, network_specs, genetic_params, dqc_compiled_encoding_circuit)
        # ----- Data Storage ----------
        dir = joinpath(folder, "warmstart_ga")
        mkpath(dir)
        serialize( joinpath(dir, "gott_encoding_circuit.jls"), gottesman_encoding_circuit )
        serialize( joinpath(dir, "gott_circuit_dqc_compiled.jls"), dqc_compiled_encoding_circuit )
        @info "Saving best-performing individual from genetic search, conditioned on `fidelity=1.0`."
        serialize(joinpath(dir, "GA_circuit.jls"), GA_encoding_circuit)
        save_txt(dir, "genetic_algorithm_parameters.txt", genetic_params)
        save_circuit_diagram(gottesman_encoding_circuit, dir, "gott_encoding_circuit.png")
        save_circuit_diagram(dqc_compiled_encoding_circuit, dir, "gott_circuit_dqc_compiled.png")
        save_circuit_diagram(GA_encoding_circuit, dir, "GA_circuit.png")
        df = DataFrame(method = ["Gottesman", "dqc_compiled", "GA"], verified = [verification_logical_state, verification_logical_state_compiled, verification_ga],
                        gate_counts = [gate_counts, gate_counts_compiled, gate_counts_ga] )
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


"""
    circuit_search_mcts(exp_label::String)::String

Orchestrate a Monte Carlo Tree Search for optimised encoding circuits. Save the resulting circuits, 
diagrams, and statistics to `data/<code_dirname>/<qpu_sizes>/<mcts>`.

### Input

- `exp_label` -- label of the experiment defined in `experiment_configurations()`

### Output 

The name of the folder to which the data is saved.

### Notes

This function requires that the serialized specification and parameter files 
(`network_specs.jls` and `code_params.jls`) have already been created for this 
experiment. If they are missing, run `create_code_network_data(exp_label)` first.
"""
function circuit_search_mcts(exp_label::String)::String
    Random.seed!(42) 
    code_architecture_setup, _, mcts_parameters = experiment_configurations()
    if haskey(code_architecture_setup, exp_label) && haskey(mcts_parameters, exp_label)
        cfg = code_architecture_setup[exp_label]
        mcts_params = mcts_parameters[exp_label]
        folder = joinpath(@__DIR__,"..","data", string(code_dirname(cfg.code)), string(cfg.qpu_sizes))
        if !isfile(joinpath(folder, "network_specs.jls")) || !isfile(joinpath(folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run `qdc_setup` with $exp_label.")
        end
        network_specs = deserialize( joinpath(folder, "network_specs.jls"))
        code_params = deserialize( joinpath(folder, "code_params.jls"))
        (
            MCTS_circuit, 
            verification_MCTS_logical_state, 
            MCTS_gate_counts, fidelity_evolution, 
            gate_count_evolution, 
            reward_evolution
        ) = monte_carlo_tree_search(code_params, network_specs, mcts_params)
        # ----- Data Storage ----------
        dir = joinpath(folder, "mcts")
        mkpath(dir)
        serialize( joinpath(dir, "MCTS_circuit.jls"), MCTS_circuit )
        save_circuit_diagram(MCTS_circuit, dir, "MCTS_circuit.png")
        save_txt(dir, "mcts_parameters.txt", mcts_params)
        df_mcts = DataFrame(
            fidelity_evolution = fidelity_evolution, 
            single_count = [gc[1] for gc in gate_count_evolution],
            two_qubit_count = [gc[2] for gc in gate_count_evolution],
            telegate_count = [gc[3] for gc in gate_count_evolution],
            reward_evolution = reward_evolution
        )
        CSV.write(joinpath(dir, "mcts_evolution.csv"), df_mcts)
        df = DataFrame(method = ["MCTS"], verified = [verification_MCTS_logical_state], gate_counts = [MCTS_gate_counts])
        CSV.write(joinpath(dir, "mcts_stats.csv"), df)
        return dir
    else
        error("The serialized specification and parameter files for this experiment are missing. Please run `qdc_setup` with $exp_label.")
    end

end


# -----------------------------------------
# ---------- DQC Execution ----------------
# -----------------------------------------

"""
    dqc_simulation(
        exp_label::String, 
        mqt_path::String, 
        circuit_path::String, 
        num_samples::Int, 
        ps::Vector{Float64}, 
        p_bells::Vector{Float64},
        telegate_idle_depth::Int,
        p_single_ratio::Float64,
        p_idle_ratio::Float64,
        method::String
    )::String

Orchestrate a DQC simulation for Fault Tolerant logical zero state encoding circuits, 
based on a given noise model. Collect simulation statistics, including the logical error rate 
of state preparation, and save to to `data/<code_dirname>/<qpu_sizes>/<simulation_(non_)FT>`.

### Input

- `exp_label` -- label of the experiment defined in `experiment_configurations()`
- `mqt_path` -- file path to the MQT encoding executable (used for FT gadget)
- `circuit_path` -- relative path (within the experiment data folder) to the `.jls` circuit file to be simulated
- `num_samples` -- total number of Monte Carlo sampling runs in the simulation
- `ps` -- array of physical initialisation, measurement and two-qubit gate noise probability p
- `p_bells` -- array of error probabilities representing the fidelity of Bell pair generation
- `telegate_idle_depth` -- the assumed depth that qubits idle while waiting for a telegate to be completed
- `p_single_ratio` -- ratio relating the physical error probability of single-qubit gates to p
- `p_idle_ratio` -- ratio relating the physical error probability of idle operations to p
- `method` -- the compilation and simulation method to be used; 
            `"none"` for non-FT encoding simulation, `"optimal"` or `"heuristic"` for FT encoding simulation

### Output 

The name of the folder to which the data is saved.

### Notes

This function requires that the serialized specification and parameter files 
(`network_specs.jls` and `code_params.jls`) have already been created for this 
experiment. If they are missing, run `create_code_network_data(exp_label)` first.
"""
function dqc_simulation(
    exp_label::String, 
    mqt_path::String, 
    circuit_path::String, 
    num_samples::Int, 
    ps::Vector{Float64}, 
    p_bells::Vector{Float64},
    telegate_idle_depth::Int,
    p_single_ratio::Float64,
    p_idle_ratio::Float64,
    method::String
)::String
    Random.seed!(42) 
    code_architecture_setup,_,_ = experiment_configurations()
    if haskey(code_architecture_setup, exp_label)
        cfg = code_architecture_setup[exp_label]
        folder = joinpath(@__DIR__,"..","data", string(code_dirname(cfg.code)), string(cfg.qpu_sizes))
        if !isfile(joinpath(folder, "network_specs.jls")) || !isfile(joinpath(folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run `qdc_setup` with $exp_label.")
        end
        network_specs = deserialize( joinpath(folder, "network_specs.jls"))
        code_params = deserialize( joinpath(folder, "code_params.jls"))
        circ_path = joinpath(folder,circuit_path) 
        encoding_circuit = deserialize(circ_path)
        if method == "none"
            data, data_circuit, DQC_circuit, full_circuit = dqc_non_ft_encoding_simulation(
                num_samples, ps, p_bells, telegate_idle_depth, p_single_ratio, 
                p_idle_ratio, code_params, network_specs, encoding_circuit)
            # ----- Data Storage ----------
            dir = joinpath(folder, "simulation_non_FT")
            mkpath(dir)
            df = DataFrame(data)
            CSV.write(joinpath(dir, "dqc_sim_data.csv"), df)
            serialize( joinpath(dir, "data_circuit.jls"), data_circuit )
            save_circuit_diagram(data_circuit, dir, "data_circuit.png")
            serialize( joinpath(dir, "DQC_circuit.jls"), DQC_circuit )
            save_circuit_diagram(DQC_circuit, dir, "DQC_circuit.png")
            serialize( joinpath(dir, "full_circuit.jls"), full_circuit )
            save_circuit_diagram(full_circuit, dir, "full_circuit.tex")
        else
            (
                data, 
                data_circuit, 
                quantum_clifford_verification_circ, 
                DQC_circuit, 
                full_circuit, 
                num_ancillas, 
                num_z_anc, 
                num_x_anc, 
                ancilla_map
            ) = dqc_ft_encoding_simulation(num_samples, ps, p_bells, telegate_idle_depth, p_single_ratio, p_idle_ratio, 
                                            code_params, network_specs, mqt_path, encoding_circuit, method)
            # ----- Data Storage ----------
            dir = joinpath(folder, "simulation_FT")
            mkpath(dir)
            df = DataFrame(data)
            CSV.write(joinpath(dir, "dqc_sim_data.csv"), df)
            serialize( joinpath(dir, "data_circuit.jls"), data_circuit )
            save_circuit_diagram(data_circuit, dir, "data_circuit.png")
            serialize( joinpath(dir, "verification_circuit.jls"), quantum_clifford_verification_circ )
            save_circuit_diagram(quantum_clifford_verification_circ, dir, "verification_circuit.png")
            serialize( joinpath(dir, "DQC_circuit.jls"), DQC_circuit )
            save_circuit_diagram(DQC_circuit, dir, "DQC_circuit.png")
            save_circuit_diagram(DQC_circuit, dir, "DQC_circuit.tex")
            serialize( joinpath(dir, "full_circuit.jls"), full_circuit )
            save_circuit_diagram(full_circuit, dir, "full_circuit.tex")
            ancilla_info = (; num_ancillas, num_z_anc, num_x_anc, ancilla_map)
            save_txt(dir, "ancilla_info.txt", ancilla_info)
            serialize(joinpath(dir, "ancilla_map.jls"), ancilla_map)
        end
        return dir
    else
        error("The configuration label $exp_label was not found. Please add the respective data to the configuration file first.")
    end
end

# -----------------------------------------
# --------- Resource Estimation -----------
# -----------------------------------------

"""
    resource_estimation(exp_label::String)::String

Estimate the resources required for (i) fault-tolerant encoding circuits and (ii) fault-tolerant distributed
stabiliser measurement-based initialisation of the logical zero state of a given experiment configuration. 
Save the statistics (ancillas, depth, gate counts, measurements) to `data/<code_dirname>/<qpu_sizes>/simulation_FT`.

### Input

- `exp_label` -- label of the experiment defined in `experiment_configurations()`

### Output 

The name of the folder to which the data is saved.

### Notes

This function requires that the serialized specification and parameter files (`network_specs.jls` and `code_params.jls`) 
as well as the DQC simulation ancilla information (`ancilla_map.jls`) have already been created for this. 
If they are missing, run `create_code_network_data(exp_label)`, `circuit_search_gott_ga(exp_label)` and `dqc_simulation(...)` first.
"""
function resource_estimation(exp_label::String)::String
    code_architecture_setup,_,_ = experiment_configurations()
    if haskey(code_architecture_setup, exp_label)
        cfg = code_architecture_setup[exp_label]
        folder = joinpath(@__DIR__,"..","data", string(code_dirname(cfg.code)), string(cfg.qpu_sizes))
        if !isfile(joinpath(folder, "network_specs.jls")) || !isfile(joinpath(folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run `qdc_setup` with $exp_label.")
        end
        network_specs = deserialize( joinpath(folder, "network_specs.jls"))
        code_params = deserialize( joinpath(folder, "code_params.jls"))
        dir = joinpath(folder, "simulation_FT")
        ancilla_map = deserialize(joinpath(dir, "ancilla_map.jls"))
        qpu_sizes = cfg.qpu_sizes .+ Int(network_specs.num_comm_qubits/network_specs.num_registers) # accounting for data + communication qubits
        (
            num_ancillas_circ, 
            circuit_depth, 
            total_gate_counts_circ, 
            total_number_measurements_circ, 
            qpu_core_sizes_circ,
            acceptance_ratio,
            p, p_bell
        ) = estimate_resources_encoding_circuit(dir, copy(qpu_sizes), ancilla_map) # FT encoding circuit resource estimate
        # ----- Data Storage ----------
        resources_info_circ = (; network_specs.num_comm_qubits, network_specs.num_registers, 
                                    num_ancillas_circ, circuit_depth, total_gate_counts_circ,
                                    total_number_measurements_circ, qpu_core_sizes_circ,
                                    acceptance_ratio, p, p_bell)
        save_txt(dir, "resources_info_circ.txt", resources_info_circ)
        resources_info_circ_df = DataFrame(
            method = "encoding_circ",
            num_comm_qubits = [network_specs.num_comm_qubits],
            num_registers = [network_specs.num_registers],
            num_ancillas = [num_ancillas_circ],
            depth_cx_layers = [circuit_depth[1]],
            depth_telegate_layers = [circuit_depth[2]],
            single_qubit_gates = [total_gate_counts_circ[1]],
            cx_gates = [total_gate_counts_circ[2]],
            telegates = [total_gate_counts_circ[3]],
            measurements = [total_number_measurements_circ],
            qpu_core_sizes = [join(qpu_core_sizes_circ, ";")],
            acceptance_ratio = acceptance_ratio,
            p=p,
            p_bell=p_bell
        )
        CSV.write(joinpath(dir, "resources_info_circ.csv"), resources_info_circ_df)
        (
            num_ancillas_meas, 
            depth_meas, 
            gate_counts_meas, 
            total_number_measurements_meas, 
            qpu_core_sizes_meas
        ) =  estimate_resources_measurement_based_encoding(network_specs, code_params, copy(qpu_sizes)) # Measurement-based encoding resource estimate
        # ----- Data Storage ----------
        resources_info_meas = (; network_specs.num_comm_qubits, network_specs.num_registers, num_ancillas_meas,
                                    depth_meas,  gate_counts_meas, total_number_measurements_meas, qpu_core_sizes_meas)
        save_txt(dir, "resources_info_meas.txt", resources_info_meas)
        resources_info_meas_df = DataFrame(
            method = "stabiliser_measurement",
            num_comm_qubits = [network_specs.num_comm_qubits],
            num_registers = [network_specs.num_registers],
            num_ancillas = [num_ancillas_meas],
            depth_cx_layers = [depth_meas[1]],
            depth_telegate_layers = [depth_meas[2]],
            single_qubit_gates = [gate_counts_meas[1]],
            cx_gates = [gate_counts_meas[2]],
            telegates = [gate_counts_meas[3]],
            measurements = [total_number_measurements_meas],
            qpu_core_sizes = [join(qpu_core_sizes_meas, ";")],
        )
        CSV.write(joinpath(dir, "resources_info_meas.csv"), resources_info_meas_df)
    else
        error("The configuration label $exp_label was not found. Please add the respective data to the configuration file first.")
    end
    return dir
end



end