# resource_estimate.jl

"""
Functions for comparing resources used in fault-tolerant encoding circuit and measurement-based 
initialisation via distributed stabiliser measurements
"""
module ResourceEstimation

export estimate_resources_encoding_circuit, estimate_resources_measurement_based_encoding

using ..Types
using ..Helper
using ..DQCodeSimulator: build_layers

using QuantumClifford: Stabilizer, AbstractOperation, sCNOT, AbstractSingleQubitOperator, AbstractTwoQubitOperator, AbstractMeasurement
using Serialization
using CSV, DataFrames


"""
    estimate_resources_encoding_circuit(folder::String, qpu_sizes::Vector{Int},
                                        ancilla_map::Vector{Int})::Tuple{Int, Vector{Int}, Vector{Int}, Int, Vector{Int}}

Load and summarise the resource requirements of the fault-tolerant encoding circuit from simulation data.

### Input

- `folder`      -- path to the directory containing `dqc_sim_data.csv`
- `qpu_sizes`   -- baseline data qubit counts per register, used to compute total qubit requirements
- `ancilla_map` -- vector mapping each ancilla qubit to its host register index

### Output

Returns a 5-element tuple containing the number of ancilla qubits, circuit depth as
`[depth_cx_layers, depth_telegate_layers]`, gate counts as `[num_single_qubit, num_cx, num_telegates]`,
number of measurements and total qubit count per register (code + communication + ancilla qubits).
"""
function estimate_resources_encoding_circuit(folder::String, qpu_sizes::Vector{Int}, 
                                                ancilla_map::Vector{Int})::Tuple{Int, Tuple{Int, Int}, Tuple{Int, Int, Int}, Int, Vector{Int}}
    df = CSV.read(joinpath(folder, "dqc_sim_data.csv"), DataFrame)
    depth_cx_layers = df.depth_cx_layers[1]
    depth_telegate_layers = df.depth_telegate_layers[1]
    total_single_qubit_count = df.total_single_qubit_count[1]
    total_two_qubit_count = df.total_two_qubit_count[1]
    total_telegate_count = df.total_telegate_count[1]
    num_meas = df.num_meas[1]
    ancilla_distr = zeros(Int, length(qpu_sizes))
    for core in ancilla_map # iterate over all entries in the ancilla map list
        ancilla_distr[core] += 1 # increase the respective core in ancilla_dist
    end
    qpu_core_sizes = copy(qpu_sizes) .+ ancilla_distr 
    @info "Gathering statistics for FT encoding circuit initialisation..."
    @info "The FT encoding circuit requires $(length(ancilla_map)) ancillas, $total_single_qubit_count single qubit gates,
            $total_two_qubit_count CX gates, $total_telegate_count telegates as well as $num_meas measurements."
    @info "The depth is ($depth_cx_layers, $depth_telegate_layers), and the total qubit counts per core is $qpu_core_sizes.\n"
    return length(ancilla_map), (depth_cx_layers, depth_telegate_layers), (total_single_qubit_count, total_two_qubit_count, total_telegate_count), num_meas, qpu_core_sizes
end
   
    
"""
    estimate_resources_measurement_based_encoding(network_specs::NetworkSpecifications, 
                                                    code_params::CodeParameters, qpu_sizes::Vector{Int})::Tuple{Int, Vector{Int}, Vector{Int}, Int, Vector{Int}}

Estimate the resource requirements for performing one round of distributed stabiliser measurements to initialise 
a QEC code (up to Pauli frame tracking based on the measurement outcome) on a Type-II distributed quantum architecture.

### Input

- `network_specs`  -- network specification object defining the Type-II DQC architecture
- `code_params`    -- code parameters defining the QEC code
- `qpu_sizes`      -- data qubit capacity of the cores 

### Outputs

Returns a 5-element tuple containing the number of ancilla qubits, circuit depth as
`[depth_cx_layers, depth_telegate_layers]`, gate counts as `[num_single_qubit, num_cx, num_telegates]`,
number of measurements and total qubit count per register (code + communication + ancilla qubits).

### Notes

We first construct the entire measurement-based initialisation circuit and count the gates in the process thereof, then we use
the circuit to build ASAP layers (+ verify the gate count) and extract the number of regular and telegate layers.

Since we are starting in the computational all-zero state, it suffices to measure the X-type stabilisers of a code. 

Furthermore, we draw upon the same assumptions as for the encoding circuit approach when it comes to ancillas:
For every measurement (i.e., every stabiliser), we utilise one fresh ancilla in order to allow for maximal parallelism, 
and we place the ancilla in the core that has the largest support in the respective stabiliser.

In this function, we are only counting the resource overhead of performing one round of distributed stabiliser measuremnts.
In order to make this encoding procedure fault-tolerant, we need to perform `code_distance` rounds of distributed stabiliser 
measurements, and then perform decoding on the spatio-temporal syndrome graph.

For simplicity, we neglect the effect of having negative phases in the `measurements` vector (which amount to single-qubit gate
corrections), since we can also re-interpret the measurement outcome for these stabiliser measurements.
"""
function estimate_resources_measurement_based_encoding(network_specs::NetworkSpecifications, 
                                                    code_params::CodeParameters, qpu_sizes::Vector{Int})::Tuple{Int, Vector{Int}, Vector{Int}, Int, Vector{Int}}
    gate_counts = [0,0,0]
    num_meas = 0
    qpu_core_sizes = copy(qpu_sizes)
    stabilisers = Stabilizer(code_params.qec_code)
    lZs = code_params.logical_Zs
    measurements = collect(stabilisers)
    append!(measurements, lZs)
    @info "Gathering statistics for measurement-based initialisation..."
    measurement_circ = Vector{AbstractOperation}()
    num_ancillas = 0
    ancilla_map = Vector{Int}()
    for (stab_meas_idx, check) in enumerate(measurements)
        msm_cores =  Dict{Int, Int}() # track core with greatest support in stabiliser (each check consumes one ancilla)
        pauli = copy(check)
        pauli_support = Vector{Int}()
        pauli_type = "X"
        for idx in collect(1:length(pauli)) 
            qubit = pauli[idx]
            if qubit != (false, false)
                core = network_specs.register_lookup_array[network_specs.inv_map[idx]]
                msm_cores[core] = get!(msm_cores, core, 0) + 1
                push!(pauli_support, idx)
                if qubit == (false, true) 
                    pauli_type = "Z"
                end
                # the `pauli_type` is consistent within each `pauli` by the virtue of the canonical form of stabiliser generators and logical Z operators of CSS codes
            end
        end
        if pauli_type=="X"
            gate_counts[1] += 2 # to account for the Hadamard gates on the ancilla
            num_meas += 1 
            num_ancillas += 1
        else 
            continue # we only apply the X type checks
        end
        mapped_core = first(keys(msm_cores))
        mapped_count = msm_cores[mapped_core] 
        for (c, cnt) in msm_cores
            if cnt > mapped_count # ties are broken in order of appearance in `msm_cores`
                mapped_core = c
                mapped_count = cnt
            end
        end
        qpu_core_sizes[mapped_core] +=1 # ancilla qubit for this `pauli` stabiliser is mapped to `mapped_core`
        push!(ancilla_map, mapped_core)
        for idx in pauli_support # traverse through supported qubit in `pauli` and determine the number of required telegates
            mapped_q = network_specs.inv_map[idx]
            if network_specs.register_lookup_array[mapped_q] == mapped_core
                gate_counts[2] += 1
            else 
                gate_counts[3] += 1
            end
            push!(measurement_circ, sCNOT(mapped_q, network_specs.num_data_qubits+stab_meas_idx))
        end
    end
    # the resulting `measurement_circ` captures only the CNOT interactions; the indices of the gates are the slots of the 
    # mapped qubits, so we determine the depth of this circuit object directly
    layers_measurement_initialisation = build_layers(measurement_circ, network_specs.num_data_qubits+num_ancillas)
    depth = [0,0] # [depth_cx_layers, depth_telegate_layers]
    gate_counts_meas = [0,0,0]
    for layer in layers_measurement_initialisation # traversing layers as in `dqc_simulator.jl` to check for depth, but without applying noise channels
        telegate_pairs = Dict{Tuple{Int,Int}, Int}()
        for gate in layer
            T = typeof(gate)
            if T <: AbstractTwoQubitOperator
                control = gate.q1 
                target =  gate.q2
                # for `measurement_circ`, the control qubit is always a data qubit, and the target is always an ancilla
                control_register = network_specs.register_lookup_array[control] 
                target_register = ancilla_map[target-network_specs.num_data_qubits] # we appended the ancillas directly after the data qubits 
                if control_register == target_register 
                    gate_counts_meas[2] += 1
                else
                    pair_key = minmax(control_register, target_register)
                    telegate_pairs[pair_key] = get(telegate_pairs, pair_key, 0) + 1
                    gate_counts_meas[3] += 1
                end # we ignore ancilla qubit measurements in this treatment since they don't contribute depth
            else
                throw("the reduced `measurement_circ` should only contain two-qubit gates")
            end
        end
        num_telegates_layer = isempty(telegate_pairs) ? 0 : maximum(values(telegate_pairs))
        if num_telegates_layer == 0
            depth[1] +=1
        else
            depth[2] += num_telegates_layer 
        end
    end
    @assert gate_counts_meas[2] == gate_counts[2] && gate_counts_meas[3] == gate_counts[3]
    @info "One round of Measurement-based encoding requires $(length(ancilla_map)) ancillas, $(gate_counts_meas[1]) single qubit gates,
            $(gate_counts_meas[2]) CX gates, $(gate_counts_meas[2]) telegates as well as $num_meas measurements."
    @info "The depth is ($(depth[1]), $(depth[2])), and the total qubit counts per core is $qpu_core_sizes.\n"
    return num_ancillas, depth, gate_counts, num_meas, qpu_core_sizes
end



end

