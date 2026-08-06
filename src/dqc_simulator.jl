# dqc_simulator.jl

"""
Functions that construct DQC-executable circuits from an optimised encoding circuit 
for a logical zero state, simulate the DQC execution including telegates under circuit-level
and Bell pair initialisation noise, and collect corresponding statistics, e.g. on logical error rates.

Credit: 
- The pipeline for the evaluation of the logical rate (`dqc_logical_evaluation`) is largely based on the implementation
provided as part of the QuantumClifford library (https://github.com/QuantumSavory/QuantumClifford.jl/blob/2664eba07a461441ea051a17cad7c725f9576176/src/ecc/decoder_pipeline.jl),
with functions for the noiseless encoding circuits being based on https://github.com/QuantumSavory/QuantumClifford.jl/blob/master/src/ecc/circuits.jl.
[The QuantumClifford.jl library is licensed under a MIT license.]
- The `build_layers` function is largely based on the `collect_circuit_layers` function exposed in the Munich Quantum Toolkit -- QECC library
(https://github.com/munich-quantum-toolkit/qecc/blob/b74a88d2e521270ff00089a77927a28290d4c59a/src/mqt/qecc/circuit_synthesis/circuit_utils.py).
[The MQT QECC library is licensed under a MIT license.]
"""
module DQCodeSimulator

export dqc_ft_encoding_simulation, dqc_non_ft_encoding_simulation
export construct_DQC_executable_circuit
export dqc_logical_evaluation
export build_layers

using ..Types
using ..Helper

using QuantumClifford
using QuantumClifford: AbstractOperation, AbstractStabilizer, AbstractMeasurement,
                        apply_single_x!, apply_single_y!, apply_single_z!
using QECCore
using QuantumClifford.ECC: DecoderCorrectionGate, CSSTableDecoder, decode,  AbstractSyndromeDecoder, 
                            faults_matrix, ClassicalTableDecoder, parity_checks
using Combinatorics: combinations
using ProgressMeter
using PyCall
using StatsBase
using LinearAlgebra

import QuantumClifford: affectedqubits, applynoise! 

qiskit = pyimport("qiskit")
qasm2 = pyimport("qiskit.qasm2")


"""
    dqc_non_ft_encoding_simulation(num_samples::Int, ps::Vector{Float64}, p_bells::Vector{Float64}, 
                                        telegate_idle_depth::Int, p_single_ratio::Float64, p_idle_ratio::Float64,
                                        code_params::CodeParameters, network_specs::NetworkSpecifications, 
                                        circuit::Vector{AbstractOperation})::Tuple{Any, Vector{AbstractOperation}, Vector{AbstractOperation}, Vector{AbstractOperation}}

Simulate the raw encoding circuit retrieved from the optimisation step in a DQC environment 
under circuit-level Bell pair initialisation noise, gathering statistics about logical error rate, fidelity etc.

### Input

- `num_samples` -- number of Monte Carlo trajectories per `(p, p_bell)` point
- `ps` -- vector of two-qubit gate error rates to sweep over
- `p_bells` -- vector of Bell pair initialisation error rates to sweep over
- `telegate_idle_depth` -- number of idle noise layers applied during a telegate round
- `p_single_ratio` -- ratio of single-qubit gate error rate to `p`, i.e. `p_single = p * p_single_ratio`
- `p_idle_ratio` -- ratio of idle error rate to `p`, i.e. `p_idle = p * p_idle_ratio`
- `code_params` -- code parameters including the target logical state and QEC structure
- `network_specs` -- network topology including qubit layout, register assignments, and communication qubits
- `circuit` -- optimised encoding circuit

### Output

Returns a 4-element tuple `(data, data_circuit, DQC_circuit, full_circuit)`, where `data` gathers the relevant
evaluation data, including `logical_error_rate`, `acceptance_ratio` and gate counts, `data_circuit` is the input
encoding circuit, `DQC_circuit` is the DQC-executable circuit constructed from `data_circuit`, and `full_circuit`
is the full circuit including noise channels and noise-free stabiliser measurements.

### Notes

We first check for validity of the noiseless circuit, also with resepct to the DQC evaluation, and then perform a DQC evaluation.

The effective idle error probability during a telegate is approximated as `p_idle_telegate_layer =1-(1-p_idle)^telegate_idle_depth`.
We consider the measurement, intitialisation and two-qubit gate error probability to be equal to `p`,
and disregard the effect of classical communication and crosstalk.
"""
function dqc_non_ft_encoding_simulation(num_samples::Int, ps::Vector{Float64}, p_bells::Vector{Float64}, 
                                        telegate_idle_depth::Int, p_single_ratio::Float64, p_idle_ratio::Float64,
                                        code_params::CodeParameters, network_specs::NetworkSpecifications, 
                                        circuit::Vector{AbstractOperation})::Tuple{Any, Vector{AbstractOperation}, Vector{AbstractOperation}, Vector{AbstractOperation}}
    data_circuit = deepcopy(circuit)
    #------------ Noiseless Testing -------------
    @info "Noiseless testing of raw encoding circuit ..."
    noise_verif = NoiseSpecs(1e2,0,0,0,0,0) # performing a small number of runs without noise
    quantum_clifford_verification_circ = Vector{AbstractOperation}() # no verification circuit → non-FT
    num_ancillas = 0
    ancilla_map = Vector{Int}()
    DQC_circuit,_,_,_,_ = construct_DQC_executable_circuit(data_circuit, quantum_clifford_verification_circ, num_ancillas, ancilla_map, network_specs, noise_verif)
    initial_state = Register(one(MixedDestabilizer, network_specs.num_data_and_comm_qubits+num_ancillas),
                                network_specs.num_comm_qubits + num_ancillas)
    state, stat = mctrajectory!(initial_state, vcat(DQC_circuit, VerifyOp(code_params.target_state, network_specs.data_qubits))) 
    @assert stat == true_success_stat "Adding the verifcation circuit compromised the data circuit"
    (
        logical_failures,
        logical_error_rate, 
        avg_fidelity, 
        z_error_pre_decoding_rate, 
        x_error_pre_decoding_rate, 
        acceptance_ratio, 
        discarded_runs,
        _,_,_,_,_,_
    ) = dqc_logical_evaluation(data_circuit, quantum_clifford_verification_circ, num_ancillas, ancilla_map,
                                code_params, network_specs, noise_verif)
    @assert z_error_pre_decoding_rate == 0.0
    @assert x_error_pre_decoding_rate == 0.0
    @assert avg_fidelity == 1.0
    @assert logical_error_rate == 0.0 "Logical error rate is non-zero in noiseless setting; please check your circuit"  
    @assert acceptance_ratio == 1.0 "Not all runs accepted in noiseless setting; please check your verification circuit"
    @info "Noiseless testing of FT encoding circuit successful."
    #------------- Noisy DQC simulation --------------
    data = NamedTuple[]
    progress = Progress(length(ps) * length(p_bells); desc="(p, p_bell) sweep", dt=1)
    full_circuit = nothing
    @assert isempty(quantum_clifford_verification_circ) # we append an empty verification circuit in the non-FT setting
    for (p, p_bell) in Iterators.product(ps, p_bells)
        p_idle = p*p_idle_ratio
        p_single = p*p_single_ratio
        p_idle_telegate_layer = 1-(1-p_idle)^telegate_idle_depth # compute the telegate_layer idle-error probability
        noise_model = NoiseSpecs(num_samples,p,p_idle,p_idle_telegate_layer,p_single,p_bell) # for each (p, p_bell), initialise NoiseSpecs(...)
        (
            logical_failures, 
            logical_error_rate, 
            avg_fidelity, 
            z_error_pre_decoding_rate, 
            x_error_pre_decoding_rate, 
            acceptance_ratio, 
            discarded_runs, 
            n_samples, 
            depth, 
            encoding_circ_gate_counts, 
            gate_counts, 
            num_meas, 
            full_circuit
        ) = dqc_logical_evaluation(data_circuit, quantum_clifford_verification_circ, num_ancillas, 
                                    ancilla_map, code_params, network_specs, noise_model) 
        push!(data, (p=p, p_bell=p_bell, p_single_ratio=p_single_ratio, p_idle_ratio=p_idle_ratio, telegate_idle_depth=telegate_idle_depth, logical_error_rate=logical_error_rate, 
                    acceptance_ratio=acceptance_ratio, discarded_runs=discarded_runs, n_samples=n_samples, logical_failures=logical_failures,
                    avg_fidelity=avg_fidelity, z_error_pre_decoding_rate=z_error_pre_decoding_rate, x_error_pre_decoding_rate=x_error_pre_decoding_rate,
                    depth_cx_layers=depth[1], depth_telegate_layers=depth[2], encoding_circ_gate_counts = encoding_circ_gate_counts,
                    total_single_qubit_count=gate_counts[1], total_two_qubit_count=gate_counts[2], total_telegate_count=gate_counts[3], num_meas=num_meas))
        next!(progress; showvalues=[(:p, p), (:p_bell, p_bell), (:logical_error_rate, logical_error_rate), (:acceptance_ratio,acceptance_ratio)])
    end
    return data, data_circuit, DQC_circuit, full_circuit
end



"""
    dqc_ft_encoding_simulation(num_samples::Int, ps::Vector{Float64}, p_bells::Vector{Float64}, 
                                    telegate_idle_depth::Int, p_single_ratio::Float64, p_idle_ratio::Float64,
                                    code_params::CodeParameters, network_specs::NetworkSpecifications, mqt_path::String,
                                    circuit::Vector{AbstractOperation}, method::String)::Tuple{Any, Vector{AbstractOperation}, Vector{AbstractOperation}, Vector{AbstractOperation}, 
                                                                                                    Vector{AbstractOperation}, Int, Int, Int, Vector{Int}}

Append the optimised encoding circuit with a FT gadget and simulate the resulting fault-tolerant encoding circuit
in a DQC environment under circuit-level Bell pair initialisation noise, gathering statistics about
logical error rate, fidelity etc.

### Input

- `num_samples` -- number of Monte Carlo trajectories per `(p, p_bell)` point
- `ps` -- vector of two-qubit gate error rates to sweep over
- `p_bells` -- vector of Bell pair initialisation error rates to sweep over
- `telegate_idle_depth` -- number of idle noise layers applied during a telegate round
- `p_single_ratio` -- ratio of single-qubit gate error rate to `p`, i.e. `p_single = p * p_single_ratio`
- `p_idle_ratio` -- ratio of idle error rate to `p`, i.e. `p_idle = p * p_idle_ratio`
- `code_params` -- code parameters including the target logical state and QEC structure
- `network_specs` -- network topology including qubit layout, register assignments, and communication qubits
- `circuit` -- optimised encoding circuit
- `method` -- Method for verification circuit synthesis, to be passed to MQT QECC library

### Output

Returns a 8-element tuple, where `data` gathers the relevant evaluation data, including `logical_error_rate`, 
`acceptance_ratio` and gate counts, `data_circuit` is the input encoding circuit, `DQC_circuit` is the
DQC-executable circuit constructed from `data_circuit`, `full_circuit` is the full circuit including noise 
channels and noise-free stabiliser measurements, and the remaining elements collect data about the verification 
ancilla register.

### Notes

# Procedure: We pass a qasm string with the optimised encoding circuit, and get back a qasm string with the verification 
    #(in MQT QECC, the result was a qiskit circuit, which was then converted to qasm in order to make the output stream usable, and then converted to qiskit here again)
    # The verification measurements are composed of some stabilisers of the zero state (X checks, Z checks and logical Zs, which are all stabilising), correcting for hook errors
    # our initial definition of logical still applies and can be used for the noiseless extraction of logical Z observable error rate
We then parse through this, doing a lot of indexing acrobatics to finally have the encoding circuit in our desired format, where
data qubits come first, then communication qubits, then ancilla qubits (and lastly ancillas for the noisefree syndrome extraction)

We first check for validity of the noiseless circuit, also with resepct to the DQC evaluation, and then perform a DQC evaluation.

The effective idle error probability during a telegate is approximated as `1-(1-p_idle)^telegate_idle_depth`.
We consider the measurement, intitialisation and two-qubit gate error probability to be equal to `p`,
and disregard the effect of crosstalk.
"""
function dqc_ft_encoding_simulation(num_samples::Int, ps::Vector{Float64}, p_bells::Vector{Float64}, 
                                    telegate_idle_depth::Int, p_single_ratio::Float64, p_idle_ratio::Float64,
                                    code_params::CodeParameters, network_specs::NetworkSpecifications, mqt_path::String,
                                    circuit::Vector{AbstractOperation}, method::String)::Tuple{Any, Vector{AbstractOperation}, Vector{AbstractOperation}, Vector{AbstractOperation}, 
                                                                                                    Vector{AbstractOperation}, Int, Int, Int, Vector{Int}}
    data_circuit = deepcopy(circuit)
    qasm = qc_circuit_to_qasm(circuit)
    python_bin = joinpath(mqt_path, ".venv/bin/python3")
    script_path =  joinpath(mqt_path, "scripts/verification_circuit.py")
    # --------- Verification Circuit --------------
    @info "Retrieving verification circuit..."
    verification_circ_qasm = readchomp(`$(python_bin) $(script_path) $qasm $(code_params.distance) $(method)`)
    verification_circ = qasm2.loads(verification_circ_qasm) 
    @info "Retrieved verification circuit"
    quantum_clifford_verification_circ = Vector{AbstractOperation}()
    ancilla_data_interactions = Dict{Tuple{String, Int}, Vector{Int}}() # will contain ancilla qubits, e.g., (Z_anc, 1) as keys, and interacting data qubit indices as values   
    num_q = 0
    num_z_anc = 0
    num_x_anc = 0
    num_flags = 0
    for reg in verification_circ.qregs
        if reg.name == "q"
            num_q = reg.size
        elseif reg.name == "z_anc"
            num_z_anc = reg.size
        elseif reg.name == "x_anc"
            num_x_anc = reg.size
        elseif reg.name == "flag"
            num_flags = reg.size # this is an upper bound, as set in MQT QECC (will be truncated later)
        end
    end
    @assert num_q == network_specs.num_data_qubits
    for instruction in verification_circ.data # parse through qiskit circuit as returned from MQT QECC library
        gate = instruction.operation.name
        if gate == "h"
            qubits = instruction.qubits
            bit_info = verification_circ.find_bit(qubits[1])
            reg_name = String(bit_info.registers[1][1].name)
            if reg_name == "q"
                continue
            end
            index = Int(bit_info.index) # refers to (Python) indexing qubit taking into account data and ancilla qubits (from MQT QECC)
            push!(quantum_clifford_verification_circ, sHadamard(index+network_specs.num_comm_qubits+1)) # need to add number of communication qubits
        elseif gate == "cx"
            qubits = instruction.qubits
            control_info = verification_circ.find_bit(qubits[1])
            control_reg_name = String(control_info.registers[1][1].name)
            control = Int(control_info.index)
            target_info = verification_circ.find_bit(qubits[2])
            target_reg_name = String(target_info.registers[1][1].name) # holds the name of the register (`target_index_reg` later holds the index within this reigster)
            target = Int(target_info.index)
            if control_reg_name == "q" && target_reg_name == "q" # skip encoding circuit instructions (we are only interested in the verification circuit part)
                continue
            elseif control_reg_name == "q" # control qubit is data qubit, target qubit is ancilla
                if target_reg_name == "z_anc"
                    target_index_reg = target-network_specs.num_data_qubits # index within `z_anc` register
                elseif target_reg_name == "x_anc"
                    target_index_reg = target-network_specs.num_data_qubits-num_z_anc # index within `x_anc` register
                elseif target_reg_name == "flag"
                    target_index_reg = target-network_specs.num_data_qubits-num_z_anc-num_x_anc # index within `flag` register
                end
                # record information about interaction between data and ancilla qubits
                key = (target_reg_name,target_index_reg+1)
                push!(get!(ancilla_data_interactions, key, Int[]), network_specs.register_lookup_array[network_specs.inv_map[control+1]])
                push!(quantum_clifford_verification_circ, sCNOT(control+1,target+network_specs.num_comm_qubits+1))
            elseif target_reg_name == "q" # control qubit is ancilla, target qubit is data qubit           
                if control_reg_name == "z_anc"
                    control_index_reg = control-network_specs.num_data_qubits # index within `z_anc` register
                elseif control_reg_name == "x_anc"
                    control_index_reg = control-network_specs.num_data_qubits-num_z_anc # index within `x_anc` register
                elseif control_reg_name == "flag"
                    control_index_reg =  control-network_specs.num_data_qubits-num_z_anc-num_x_anc # index within `flag` register
                end
                # record information about interaction between data and ancilla qubits
                key = (control_reg_name,control_index_reg+1)
                push!(get!(ancilla_data_interactions, key, Int[]), network_specs.register_lookup_array[network_specs.inv_map[target+1]])
                push!(quantum_clifford_verification_circ, sCNOT(control+network_specs.num_comm_qubits+1,target+1))
            else # both qubits are ancilla qubits
                if control_reg_name == "z_anc"
                    control_index_reg = control-network_specs.num_data_qubits # index within `z_anc` register
                elseif control_reg_name == "x_anc"
                    control_index_reg = control-network_specs.num_data_qubits-num_z_anc # index within `x_anc` register
                elseif control_reg_name == "flag"
                    control_index_reg = control-network_specs.num_data_qubits-num_z_anc-num_x_anc # index within `flag` register
                end
                if target_reg_name == "z_anc"
                    target_index_reg = target-network_specs.num_data_qubits # index within `z_anc` register
                elseif target_reg_name == "x_anc"
                    target_index_reg = target-network_specs.num_data_qubits-num_z_anc # index within `x_anc` register
                elseif target_reg_name == "flag"
                    target_index_reg = target-network_specs.num_data_qubits-num_z_anc-num_x_anc # index within `flag` register
                end
                control_key = (control_reg_name,control_index_reg+1)
                target_key = (target_reg_name,target_index_reg+1)
                # if both qubits are ancilla qubits, we still record information about their mapping by appending
                # the current mode of each ancilla to the dictionary of the other ancilla
                control_interactions = get!(ancilla_data_interactions, control_key, Int[])
                control_mode = isempty(control_interactions) ? 1 : mode(control_interactions)
                target_interactions = get!(ancilla_data_interactions, target_key, Int[])
                target_mode = isempty(target_interactions) ? 1 : mode(target_interactions)
                push!(get!(ancilla_data_interactions, control_key, Int[]), target_mode )
                push!(get!(ancilla_data_interactions, target_key, Int[]), control_mode )
                push!(quantum_clifford_verification_circ, sCNOT(control+network_specs.num_comm_qubits+1,target+network_specs.num_comm_qubits+1))
            end
        elseif gate=="measure" 
            qubits = instruction.qubits
            bit_info = verification_circ.find_bit(qubits[1])
            reg_name = String(bit_info.registers[1][1].name)
            @assert reg_name != "q" # only ancillas will be measured
            index = Int(bit_info.index)
            push!(quantum_clifford_verification_circ, sMZ(index+network_specs.num_comm_qubits+1, index-network_specs.num_data_qubits+network_specs.num_comm_qubits+1)) 
        else
            throw("Verification circuit contains gates that are currently not supported.")
        end
    end
    # determine the optimal placement for the ancilla qubits based on the number of interactions with the QPUs, thereby minimising telegate overhead
    ancilla_map = Vector{Int}() 
    anc_order = Dict("z_anc" => 1, "x_anc" => 2, "flag" => 3)
    sorted_dict = sort(collect(ancilla_data_interactions); by = x -> (anc_order[x[1][1]], x[1][2]))
    for (index, (ancilla, interactions)) in enumerate(sorted_dict)
        push!(ancilla_map, mode(interactions))
    end
    num_ancillas = length(ancilla_map) # collects all ancillas that have partaken in some interaction
    @info "Verification circuit uses $num_ancillas ancillas, with $num_z_anc z-ancillas, $num_x_anc x-ancillas, and $(num_ancillas-num_z_anc-num_x_anc) flag qubits"
    #------------ Noiseless Testing -------------
    @info "Noiseless testing of FT encoding circuit ..."
    noise_verif = NoiseSpecs(1e2,0,0,0,0,0) # perform a small number of runs without noise
    DQC_circuit,_,_,_,_ = construct_DQC_executable_circuit(data_circuit, quantum_clifford_verification_circ, num_ancillas, 
                                                            ancilla_map, network_specs, noise_verif)
    initial_state = Register(one(MixedDestabilizer, network_specs.num_data_and_comm_qubits+num_ancillas),
                                network_specs.num_comm_qubits + num_ancillas)
    state, stat = mctrajectory!(initial_state, vcat(DQC_circuit, VerifyOp(code_params.target_state, network_specs.data_qubits))) 
    @assert stat == true_success_stat "Adding the verifcation circuit compromised the data circuit"
    ( 
        logical_failures, 
        logical_error_rate, 
        avg_fidelity, 
        z_error_pre_decoding_rate, 
        x_error_pre_decoding_rate, 
        acceptance_ratio, discarded_runs,
        _,_,_,_,_,_
    ) = dqc_logical_evaluation(data_circuit, quantum_clifford_verification_circ, num_ancillas, 
                                ancilla_map, code_params, network_specs, noise_verif)
    @assert z_error_pre_decoding_rate == 0.0
    @assert x_error_pre_decoding_rate == 0.0
    @assert avg_fidelity == 1.0
    @assert logical_error_rate == 0.0 "Logical error rate is non-zero in noiseless setting; please check your circuit"  
    @assert acceptance_ratio == 1.0 "Not all runs accepted in noiseless setting; please check your verification circuit"
    @info "Noiseless testing of FT encoding circuit successful."
    #------------- Noisy DQC simulation --------------
    data = NamedTuple[]
    progress = Progress(length(ps) * length(p_bells); desc="(p, p_bell) sweep", dt=1)
    full_circuit = nothing
    for (p, p_bell) in Iterators.product(ps, p_bells)
        p_idle = p*p_idle_ratio
        p_single = p*p_single_ratio
        p_idle_telegate_layer = 1-(1-p_idle)^telegate_idle_depth 
        noise_model = NoiseSpecs(num_samples,p,p_idle,p_idle_telegate_layer,p_single,p_bell) # for each (p, p_bell), initialise NoiseSpecs(...)
        (
            logical_failures, 
            logical_error_rate, 
            avg_fidelity, 
            z_error_pre_decoding_rate, 
            x_error_pre_decoding_rate, 
            acceptance_ratio, 
            discarded_runs, 
            n_samples, 
            depth, 
            encoding_circ_gate_counts, 
            gate_counts, 
            num_meas, 
            full_circuit
        ) = dqc_logical_evaluation(data_circuit, quantum_clifford_verification_circ, num_ancillas, 
                                    ancilla_map, code_params, network_specs, noise_model)
        push!(data, (p=p, p_bell=p_bell, p_single_ratio=p_single_ratio, p_idle_ratio=p_idle_ratio, telegate_idle_depth=telegate_idle_depth,
                    logical_error_rate=logical_error_rate, acceptance_ratio=acceptance_ratio, discarded_runs=discarded_runs, n_samples=n_samples, logical_failures=logical_failures,
                    avg_fidelity=avg_fidelity, z_error_pre_decoding_rate=z_error_pre_decoding_rate, x_error_pre_decoding_rate=x_error_pre_decoding_rate,
                    depth_cx_layers=depth[1], depth_telegate_layers=depth[2], encoding_circ_gate_counts = encoding_circ_gate_counts,
                    total_single_qubit_count=gate_counts[1], total_two_qubit_count=gate_counts[2], total_telegate_count=gate_counts[3], num_meas=num_meas))
        next!(progress; showvalues=[(:p, p), (:p_bell, p_bell), (:logical_error_rate, logical_error_rate), (:acceptance_ratio,acceptance_ratio)])
    end
    return data, data_circuit, quantum_clifford_verification_circ, DQC_circuit, full_circuit, num_ancillas, num_z_anc, num_x_anc, ancilla_map
end



"""
    dqc_logical_evaluation(data_circuit::Vector{AbstractOperation}, verification_circuit::Vector{AbstractOperation}, 
                                num_ancillas::Int, ancilla_map::Vector{Int}, code_params::CodeParameters, network_specs::NetworkSpecifications,
                                noise::NoiseSpecs)::Tuple{Vector{Int}, Float64, Float64, Float64, Float64, Float64, Int, Int, Vector{Int}, Vector{Int}, Vector{Int}, Int, Vector{AbstractOperation}}

Evaluate the performance of the FT encoding circuit under circuit-level plus Bell state initialisation error,
based on KPIs of logical error rate, fidelity, etc., extracted with noise-free decoding from noisy Monte Carlo simulations.                              

### Input

- `data_circuit` -- optimised raw encoding circuit 
- `verification_circuit` -- FT gadget / verification circuit appended after the encoding circuit
- `num_ancillas` -- number of ancilla qubits allocated for verification circuit measurements
- `ancilla_map` -- mapping of ancilla qubit indices to their assigned cores
- `code_params` -- code parameters
- `network_specs` -- network specificatoins
- `noise` -- noise model specifying the number of samples and circuit-level noise error rates as well as Bell pair initialisation error probability

### Output

Returns a 13-element tuple collecting data about the DQC simulation, including `logical_failures`, `logical_error_rate`, 
`avg_fidelity`, `acceptance_ratio`, etc. The returned `circuit` contains the DQC state preparation and verification circuit as well as 
the syndrome circuit (but not the logical Z-measurement circuit).

### Notes

In this function, we first retrieve the DQC-executable circuit `DQC_circuit` from `construct_DQC_executable_circuit`. Then,
we setup a noisefree QEC cycle / syndrome decoding routine.  For the noise-free decoding, we use a CSS Lookup Table Decoder. 

Simulating a number of Monte Carlo trajectories, we sample from the noise distribution and gather statistics about the quality of FT encoding.
Here, we first discard runs according to `verification_bits` (these runs do not further impact the fidelity or logical error count),
then identify pre-decoding X- and Z- errors, and then determine logical Z observable errors, i.e., errors caused by logical X or Y operators (note
that for the logical Z state, logical Z errors cannot occur).
Therefore, we collect the `error_guess` by providing the decoder with the measured `syndrome`. The array `error_guess` collects `n` guesses for X-errors
and `n` guesses for Z errors -- we are primarily interested in the first `n` guesses, i.e., whether the decoder predicts that a certain physical X error happened
(the second half of the array contains guesses for Z errors, which cannot induce a logical error on the logical Z observable).
If the decoder cannot infer an error guess, we register a logical error on all logical qubits, and a fidelity of `0.0` for this run (since the
resulting state lives outside of the codespace and cannot be corrected). 

If the decoder can make an error guess, we expose two different methods, the first (Z-observable inconsistency) to approximate the true logical initialisation error
rate, and the second (fidelity) as further metric for comparison.

    (i) We noiselessly measure the logical Z observables and compute the logical error rate based on the `measured_logical_Z_bits` as well as 
        `faults_matrix_z` (a `k \times 2n` matrix), which collects information about which of the physical errors (first `n` columns regarding physical `X`, 
        last `n` columns regarding physical `Z` errors, matching the ordering in `error_guess`) flip which of the `k` logical
        Z operators / logical Z observables via anti-commutation. The logic here is as follows: 
        If the j-th `sum_mod = sum_{q=1:2n) ( faults_matrix_z[j, q] * error_guess[q] ) mod 2` is 1, then there is an odd number of indices that are jointly supported by both `faults_matrix_z`
        and `error_guess` for the given `j`th logical Z operator. In this case, correcting based on `error_guess` will flip the `k`th logical Z operator. This
        is only desirable if there actually was a logical Z observable error on the `k`th logical Z observable, i.e., `sum_mod = measured_logical_Z_bits[j]`. Otherwise,
        a logical failure is recorded.
        If the j-th `sum_mod = sum_q ( faults_matrix_z[j, q] * error_guess[q] ) mod 2` is 0, then there is an even number of indices that are jointly supported by both `faults_matrix_z`
        and `error_guess` for the given `j`th logical Z operator. In this case, correcting based on `error_guess` will not flip the `k`th logical Z operator by 
        error degeneracy of the code (applying an even number of operations that anti-commute with the logical operator overall cancels out). This is desirable as 
        long as there was actually no logical Z observable error on the `k`th logical Z observable, i.e., `sum_mod = measured_logical_Z_bits[j] again. Otherwise,
        a logical failure is recorded.
        In the end, the logical error rate is derived as the mean of all `k` logical Z observable error rates, conditioned on the accepted runs.
            
    (ii) For the fidelity measurement, we apply the correction gate inferred by the LUT decoder, and then measure the fidelity of the resulting state with our target state.
        Therefore, we apply the correction gate (https://github.com/QuantumSavory/QuantumClifford.jl/blob/444f341a50d2926b16b63d98586b8b06a7b6ac10/src/ecc/decoder_correction_gate.jl)
        that maps the quantum state back to the codespace (since there is a correctable syndrome, we know that the inferred correction maps back to the codespace) and 
        call `dot` between the corrected `corrected_state_data_qubits` state and the target state `code_params.target_state`. We update fidelities (the squared dot, or inner, product)
        per sample in order to avoid having to save the fidelity for each of the samples. This is done using the update equation

            `avg_fidelity = (avg_fidelity*(samples-1)+x)/samples = avg_fidelity+(x-avg_fidelity)/samples → avg_fidelity += (x-avg_fidelity)/samples`,
        
        where `samples` is replaced with `sample-discarded_runs`. This method computes fidelity of the entire resulting state with our target state (not per logical qubit), hence we
        expect it to lower bound the mean LIER.

In the entire procedure, there are two main sources of uncertainty: the probabilistic noise sampling of errors in the circuit, and the measurement of the logical Z operators. 
In a real-life experiment, both sources of errors exist (where probabilistically measuring the correct state implies that the system indeed is in the correct state afterwards).
However, by the virtue of our noise channel and correction gates being Pauli gates, the measurement in method (i) is actually deterministic, such that only the first source of uncertainty 
(noise sampling) remains. (For the fidelity computation in method (ii), such measurement stochasticity does not exist anyways.)

Capping a simulation when a specified number of logical errors has been registered helps to decrease overall runtime (this reduces execution time at
the cost of granularity in the high-noise regime, which is reasonable).

Overall, we can then use the information extracted from this function to compare logical initialisation error rates with the physical initialisation error rate `p`.
"""
function dqc_logical_evaluation(data_circuit::Vector{AbstractOperation}, verification_circuit::Vector{AbstractOperation}, 
                                num_ancillas::Int, ancilla_map::Vector{Int}, code_params::CodeParameters, network_specs::NetworkSpecifications,
                                noise::NoiseSpecs)::Tuple{Vector{Int}, Float64, Float64, Float64, Float64, Float64, Int, Int, Vector{Int}, Vector{Int}, Vector{Int}, Int, Vector{AbstractOperation}}
    (
        DQC_circuit, 
        depth, 
        encoding_circ_gate_counts, 
        gate_counts, 
        num_meas
    ) = construct_DQC_executable_circuit(data_circuit, verification_circuit, num_ancillas, 
                                            ancilla_map, network_specs, noise)
    # Lookup Table Decoder for noiseless syndrome decoding
    css_lut_decoder = CSSTableDecoder(code_params.qec_code, error_weight =  Int(floor((code_params.distance-1)/2)) )
    # determine x-type stabilisers for downstream analysis (by the virtue of CSS codes stabiliser types in canonical form are cleanly separated)
    H = parity_checks(css_lut_decoder) # stacks `Hx` and `Hz` horizontally (X-type stabiliser checks come first)
    nrows = length(H)
    ncols = nqubits(H)
    Z_type_indices = Set{Int}()
    X_type_indices = Set{Int}()
    for i in 1:nrows, j in 1:ncols
        if H[i][j] == (false, true) 
            push!(Z_type_indices,i)  
        elseif H[i][j] == (true, false) 
            push!(X_type_indices,i)
        end
    end 
    @assert maximum(X_type_indices) < minimum(Z_type_indices) "the X-type parity checks should be listed first to guarantee correct functioning of our syndrome extraction pipeline"
    @assert setdiff(Z_type_indices, X_type_indices) == Z_type_indices "there were mixed Pauli strings amongst the stabiliser generators - please check the CSS code structure"
    # determine `faults_matrix_z` -- the fault matrix is a (2k)x(2n) dimensional matrix → last k rows specify logical Z part
    faults_matrix_z = css_lut_decoder.faults_matrix[end÷2+1:end,:] 
    k = size(faults_matrix_z, 1) 
    two_n = size(faults_matrix_z, 2) # two_n=2*code_params.n
    @assert k == code_params.k
    @assert two_n == 2*code_params.n
    noisefree_syndrome_circ, num_noisefree_syndrome_ancillas, noisefree_syndrome_bits = syndrome_circuit(H, network_specs.num_data_and_comm_qubits + num_ancillas + 1, 
                                                                                                            network_specs.num_comm_qubits + num_ancillas + 1, network_specs)
    # by ordering of `H`, X-type syndrome measurements are performed first
    noisefree_logical_Z_circ, num_noisefree_logical_Z_ancillas, noisefree_logical_Z_bits = syndrome_circuit(code_params.logical_Zs,
                                                                                                            network_specs.num_data_and_comm_qubits + num_ancillas + num_noisefree_syndrome_ancillas + 1, 
                                                                                                            last(noisefree_syndrome_bits)+1, network_specs )
    total_number_qubits = network_specs.num_data_and_comm_qubits + num_ancillas +  num_noisefree_syndrome_ancillas + num_noisefree_logical_Z_ancillas
    total_number_classical_regs = network_specs.num_comm_qubits + num_ancillas + num_noisefree_syndrome_ancillas + num_noisefree_logical_Z_ancillas
    circuit = vcat(DQC_circuit, noisefree_syndrome_circ)
    discarded_runs = 0
    z_errors_pre_decoding = 0
    x_errors_pre_decoding = 0
    logical_failures = zeros(code_params.k) # track logical Z-observable failures per qubit
    avg_fidelity = 0.0
    n_samples = noise.n_samples
    initial_state = Register(one(MixedDestabilizer,total_number_qubits), total_number_classical_regs)
    for sample in 1:n_samples 
        state_post_syndrome, stats = mctrajectory!(copy(initial_state), circuit)
        verification_bits = @view state_post_syndrome.bits[network_specs.num_comm_qubits+1:network_specs.num_comm_qubits+num_ancillas] # extract verification bits
        syndrome = @view state_post_syndrome.bits[noisefree_syndrome_bits] # extract syndrome bits
        if any(verification_bits)  # determine whether or not to discard the run (any bit is `true` → verification stabiliser measurement was triggered, indicating a harmful error)
            discarded_runs +=1
            continue
        end
        state_post_logical, stats = mctrajectory!(copy(state_post_syndrome), noisefree_logical_Z_circ)
        measured_logical_Z_bits = @view state_post_logical.bits[noisefree_logical_Z_bits] # extract logical Z bits
        if any(syndrome[collect(Z_type_indices)]) || any(measured_logical_Z_bits)
            x_errors_pre_decoding +=1 # track pre-encoding X error if any Z-type stabiliser among the parity checks or the logical Z operators is triggered
        end
        if any(syndrome[collect(X_type_indices)]) 
            z_errors_pre_decoding +=1 # track pre-encoding Z error if any X-type stabiliser among the parity checks is triggered
        end
        error_guess = decode(css_lut_decoder, syndrome) # retrieve the error guess from the LUT decoder, given the syndrome
        # `decode` function expected syndrome of form [`syndrome X-type checks` | `syndrome Z-type checks`], which is satisfied by the ordering of `H`
        # `error_guess` is of the form (guess_x, guess_z), where `guess_x` is the error guess of X errors, and `guess_z` for Z errors
        if isnothing(error_guess) # if no error guess can be extracted from the lookup table, we collect a logical error on all logical qubits and fidelity 0.0
            logical_failures .+= 1
            avg_fidelity -= avg_fidelity / (sample-discarded_runs) # if execution reaches this line, we know that sample \neq discarded_runs holds
            continue
        end
        #-------- (i) Analytical logical error rate computation --------
        for j in 1:k # iterate over the k logical Z operators
            sum_mod = 0
            # iterate over all the `two_n=2*code_params.n` physical qubits error locations, the physical Z-errors stored at indices `code_params.n+1:two_n` should be zero, since they commute with the logical Z operators
            @inbounds @simd for q in 1:two_n
                sum_mod += faults_matrix_z[j, q] * error_guess[q] 
            end
            sum_mod %= 2 
            if sum_mod != measured_logical_Z_bits[j] # `measured_logical_Z_bits` captures the measured Z-faults
                logical_failures[j] += 1
            end
        end
        #-------- (ii) Simulated correction gate --------
        correction_gate = DecoderCorrectionGate(css_lut_decoder, network_specs.data_qubits, noisefree_syndrome_bits ) 
        corrected_state,_ = mctrajectory!(copy(state_post_syndrome),[correction_gate]) 
        corrected_state_data_qubits = stabilizerview( ptrace(copy(corrected_state.stab), collect(network_specs.num_data_qubits+1:total_number_qubits)) )
        post_decoding_fidelity = dot(corrected_state_data_qubits, code_params.target_state)^2
        avg_fidelity += (post_decoding_fidelity-avg_fidelity)/(sample-discarded_runs)
        if mean(logical_failures) >= 500
            @info "For physical error rate $(noise.p) and Bell state error rate $(noise.p_bell), we collected $(maximum(logical_failures)) logical failures after $sample iterations."
            n_samples = sample
            break
        end
    end
    acceptance_ratio = 1 - discarded_runs/n_samples
    logical_error_rate = mean(logical_failures)/(n_samples-discarded_runs) 
    z_error_pre_decoding_rate = z_errors_pre_decoding/(n_samples-discarded_runs)
    x_error_pre_decoding_rate = x_errors_pre_decoding/(n_samples-discarded_runs)
    return logical_failures, logical_error_rate, avg_fidelity, z_error_pre_decoding_rate, x_error_pre_decoding_rate, acceptance_ratio, discarded_runs, n_samples, depth, encoding_circ_gate_counts, gate_counts, num_meas, circuit
end


"""
    syndrome_circuit(parity_check_tableau::Stabilizer, ancillary_index::Int, bit_index::Int,
                            network_specs::NetworkSpecifications)::Tuple{Vector{AbstractOperation}, Int, UnitRange{Int64} }

Construct a syndrome measurement circuit for a given parity check tableau.

### Input

- `parity_check_tableau` -- tableua of stabiliser generators/parity checks to be measured
- `ancillary_index` -- starting index for ancilla qubits in the global qubit register
- `bit_index` -- starting index for classical bits in the classical bit register
- `network_specs` -- network specfications

### Output

Returns a 3-element tuple, where `syndrome_circ` is the resulting syndrome measurement circuit,
`ancillaries` counts the number of ancilla qubits consumed, and `bit_range` encodes the indices of 
the syndrome bits in the classical bit register.

### Notes

We iterate over each check in `parity_check_tableau` and append the corresponding ancilla-based
Pauli measurement, offsetting ancilla and classical bit indices accordingly.
"""
function syndrome_circuit(parity_check_tableau::Stabilizer, ancillary_index::Int, bit_index::Int,
                            network_specs::NetworkSpecifications)::Tuple{Vector{AbstractOperation}, Int, UnitRange{Int64} }
    syndrome_circ = AbstractOperation[]
    ancillaries = 0
    bits = 0
    for check in parity_check_tableau
        append!(syndrome_circ, perfect_ancillary_pauli_measurement(check, ancillary_index+ancillaries, bit_index+bits, network_specs))
        ancillaries +=1
        bits +=1
    end
    return syndrome_circ, ancillaries, bit_index:bit_index+bits-1
end


"""
    perfect_ancillary_pauli_measurement(p::PauliOperator, ancillary_index::Int, bit_index::Int, network_specs::NetworkSpecifications)::Vector{AbstractOperation}

Construct a noiseless ancilla-based circuit to measure a single Pauli operator `p` via
an auxiliary qubit, recording the outcome in a classical bit.

### Input

- `p` -- Pauli operator to be measured
- `ancillary_index` -- index of the ancilla qubit in the qubit register
- `bit_index` -- index of the classical bit to which the measurement outcome is written
- `network_specs` -- network specifications

### Output

Returns the circuit performing the perfect ancillary Pauli measurement.

### Notes

This circuit acts on the original qubit indexing. We use an arbitrary gate set (`sXCX`, `sCNOT` and `sYCX`)
as well as a deterministic phase correction gate since only validity (not implementability) is of concern 
here (in `resource_estimation.jl`, the syndrome circuit overhead is analysed).
"""
function perfect_ancillary_pauli_measurement(p::PauliOperator, ancillary_index::Int, bit_index::Int, network_specs::NetworkSpecifications)::Vector{AbstractOperation}
    circuit = AbstractOperation[]
    num_data_qubits = nqubits(p)
    @assert num_data_qubits == network_specs.num_data_qubits
    for qubit in 1:num_data_qubits
        if p[qubit] == (1,0)
            push!(circuit, sXCX(qubit, ancillary_index)) # X-controlled X     
        elseif p[qubit] == (0,1)
            push!(circuit, sCNOT(qubit, ancillary_index)) # Z-controlled X
        elseif p[qubit] == (1,1)
            push!(circuit, sYCX(qubit, ancillary_index)) # Y-controlled X
        end
    end
    p.phase[] == 0 || push!(circuit, sX(ancillary_index))
    mz = sMZ(ancillary_index, bit_index) # don't need to reset the ancilla, since every stabiliser measurement is measured separately
    push!(circuit, mz)
    return circuit
end



"""
    construct_DQC_executable_circuit(data_circuit::Vector{AbstractOperation}, verification_circuit::Vector{AbstractOperation}, 
                                            num_ancillas::Int, ancilla_map::Vector{Int}, 
                                            n::NetworkSpecifications, noise::NoiseSpecs)::Tuple{Vector{AbstractOperation}, Vector{Int}, Vector{Int}, Vector{Int}, Int}

Construct a DQC-executable circuit from the optimised encoding circuit `data_circuit` and the FT gadget `verification_circuit`,
applying circuit-level and Bell state initialisation noise according to the noise specifications in `noise`.

### Input

- `data_circuit` -- optimised raw encoding circuit 
- `verification_circuit` -- FT gadget / verification circuit appended after the encoding circuit
- `num_ancillas` -- number of ancilla qubits allocated for verification circuit measurements
- `ancilla_map` -- mapping of ancilla qubit indices to their assigned cores
- `n` -- network specificatoins
- `noise` -- noise model specifying the number of samples and circuit-level noise error rates as well as Bell pair initialisation error probability

### Output

Returns a 5-tuple containing the executable `circuit`` and statistics about the circuit (`depth`, `encoding_circ_gate_counts`, `gate_counts`, `num_meas`).

### Notes

We then treat `data_circuit` and `verification_circuit` separately, for each building an ASAP scheduling of circuit operations based on 
the data qubits defining the code, and then performing the operations per `layer`. Here, we need to take into account that even though 
per layer, the affected qubits of gate operations are disjoint by design, there can be conflicts regarding the communication qubits in use.
Since we do not allocate one communication qubit per data qubit, but merely assign one communication qubit for each other register, per register,
this requires us to track a dictionary of `telegate_pairs`.

When building the overall circuit, we append the telegate gadget via `add_telegate` and add the respective noise channels (consisting of 
initialisation, measurement, single- and two-qubit gate, idling and telegate-layer idling errors).
Here, for every layer, we count from how many consecutive telegate creation wait-time layers (`num_telegates_in_layer`) all non-telegate qubits suffer, and 
add the corresponding idling noise with probability `1-(1-noise.p_idle_telegate_layer)^num_telegates_in_layer)`. The same idling noise on the 
telegate qubits that need to wait due to communication qubit conflicts is applied dynamically.

We eventually seek to compare the logical initialisation error rate with the physical. Since we are exclusively initialising in the 
zero state, and measuring in the Z basis, when applying a depolarising channel with error probability `p`, out of the three Pauli errors {`X`,`Y`,`Z`},
only `X` and `Y` lead to a physical initialisation or measurement error. Therefore, our noise model captures an `effective' initialisation and 
measurement error of `2/3*p`.

In `n.mapping`, value `j` at index `i` indicates that qubit `j` is mapped to physical slot `i`. As discussed in `helper.jl`, 
the list of transpositions `n.mapping_transpositions` (created by the function `perm_to_transpositions`), can be applied from
right to left to create this DQC mapping. This amount to applying `sSWAP` operations in the reverse order of `n.mapping_transpositions` 
(where traversing the transpositions from right to left implements the permutation by construction), before adding the circuit operations
in `data_circuit` (the non-FT encoding circuit) and `verification_circuit` (the FT gadget), and then undoing these SWAPS by applying the  mirror-image afterwards. 
To identify the DQC index of an operation, we then use `n.inv_map[op.q]`, where `op.q` is the unmapped index of the affected qubit in quantum operation `q` (to 
apply `sHadamard(1)`, for example, we need to identify the DQC slot to which qubit `1` has been mapped; this is achieved by retrieving index `1`  of the 
inverse permutation of `n.mapping`). In general, data qubits sit at indices `1:n.num_data_qubits` and are mapped onto the DQC architecture for visualisation
purposes, while communication qubits and ancilla qubits are appended to the indices following index `n.num_data_qubits`, while tracking their core assignment with `ancilla_map`.
"""

function construct_DQC_executable_circuit(data_circuit::Vector{AbstractOperation}, verification_circuit::Vector{AbstractOperation}, 
                                            num_ancillas::Int, ancilla_map::Vector{Int}, 
                                            n::NetworkSpecifications, noise::NoiseSpecs)::Tuple{Vector{AbstractOperation}, Vector{Int}, Vector{Int}, Vector{Int}, Int}
    circuit = Vector{AbstractOperation}() # initialise the DQC-executable circuit
    for (i,j) in reverse(n.mapping_transpositions) # DQC mapping
        push!(circuit, sSWAP(i,j))  
    end     
    #----------- Make encoding circuit DQC-executable ------------   
    layers_enc_circ = build_layers(data_circuit, n.num_data_qubits)    
    add_noise(circuit, [n.inv_map[data_q] for data_q in collect(1:n.num_data_qubits)], noise.p) # initialisation error probability `p`
    depth = [0,0] # (# regular layers, # telegates layers)
    gate_counts = [0,0,0]
    num_meas = 0
    for layer in layers_enc_circ    
        affected_qubits = Set( Iterators.flatten( [affectedqubits(gate) for gate in layer] ) ) 
        idle_qubits = setdiff(1:n.num_data_qubits, affected_qubits)
        idle_qubits_DQC = [n.inv_map[idle_q] for idle_q in idle_qubits] # DQC indices of idling qubits in `layer`    
        telegate_pairs = Dict{Tuple{Int,Int}, Int}()
        telegate_qubits = Set{Int64}()
        for gate in layer
            T = typeof(gate)
            if T <: AbstractSingleQubitOperator
                qubit = affectedqubits(gate)[1]
                DQC_qubit = n.inv_map[qubit]
                push!(circuit, sHadamard(DQC_qubit))
                gate_counts[1] += 1
                add_noise(circuit, [DQC_qubit], noise.p_single) # single-qubit gate error probability `p_single`
            elseif T <: AbstractTwoQubitOperator
                control = affectedqubits(gate)[1]
                target = affectedqubits(gate)[2]
                DQC_control = n.inv_map[control]
                DQC_target = n.inv_map[target]
                control_register = n.register_lookup_array[DQC_control] 
                target_register = n.register_lookup_array[DQC_target] 
                if control_register == target_register 
                    push!(circuit, sCNOT(DQC_control, DQC_target))
                    gate_counts[2] += 1
                    add_noise(circuit, [DQC_control, DQC_target], noise.p; two_qubits = true) # two-qubit gate error probability `p`
                else # telegate
                    pair_key = minmax(control_register, target_register)
                    prev_consecutive_telegates_pair_key = get(telegate_pairs, pair_key, 0)
                    prev_max_telegate_layers = isempty(telegate_pairs) ? 0 : maximum(values(telegate_pairs))
                    # noise before the telegate application, due to other telegates using the same communication qubit pair
                    add_noise(circuit, [DQC_control, DQC_target], 1-(1-noise.p_idle_telegate_layer)^prev_consecutive_telegates_pair_key ) 
                    telegate_pairs[pair_key] = prev_consecutive_telegates_pair_key + 1
                    circuit = add_telegate(circuit, DQC_control, DQC_target, control_register, target_register, n, noise) # perform telegate between qubits
                    gate_counts[3] += 1                    
                    if prev_consecutive_telegates_pair_key + 1 > prev_max_telegate_layers # the updated count in `telegate_pairs[pair_key]` is the new maximum
                        add_noise(circuit, collect(telegate_qubits), noise.p_idle_telegate_layer) # post-apply one layer of `p_idle_telegate_layer` noise on all other `telegate_qubits` encountered so far
                    else
                        # post-apply `maximum(values(telegate_pairs))-(prev_consecutive_telegates_pair_key+1)` layers of `p_idle_telegate_layer` noise on the current `DQC_control` and `DQC_target`
                        @assert prev_max_telegate_layers == maximum(values(telegate_pairs))
                        add_noise(circuit, [DQC_control, DQC_target], 1-(1-noise.p_idle_telegate_layer)^(prev_max_telegate_layers-prev_consecutive_telegates_pair_key-1)) 
                    end
                    push!(telegate_qubits, DQC_control)
                    push!(telegate_qubits, DQC_target)
                end
            else
                throw("Circuit contains gates that have not been classified as Single- or Two-Qubit gate so far.")
            end
        end
        num_telegates_in_layer = isempty(telegate_pairs) ? 0 : maximum(values(telegate_pairs))
        if num_telegates_in_layer == 0
            add_noise(circuit, idle_qubits_DQC, noise.p_idle) 
            depth[1] +=1
        else
            # `telegate_qubits` are excluded from this idling noise application, since idling noise has been applied to them dynamically already
            add_noise(circuit, setdiff(1:n.num_data_qubits, Set(telegate_qubits)), 1-(1-noise.p_idle_telegate_layer)^num_telegates_in_layer)  
            depth[2] += num_telegates_in_layer 
        end
    end
    #----------- Make verification gadget DQC-executable ------------
    encoding_circ_gate_counts = copy(gate_counts)
    layers_ver_circ = build_layers(verification_circuit, n.num_data_and_comm_qubits + num_ancillas)
    all_qubits = 1:(n.num_data_and_comm_qubits + num_ancillas)
    ancilla_qubits = setdiff(all_qubits, 1:n.num_data_and_comm_qubits)
    add_noise(circuit, ancilla_qubits, noise.p) #  init noise on ancilla qubits used for verification measurements
    for layer in layers_ver_circ # analogous to traversal of `layers_enc_circ`
        affected_qubits = Set( Iterators.flatten( [affectedqubits(gate) for gate in layer] ) ) 
        idle_qubits = setdiff(all_qubits, union(Set(n.comm_qubits), affected_qubits))
        idle_qubits_DQC = Vector{Int}()
        for idle_q in idle_qubits
            # comm qubits are not treated as idle, since they get reinitialised upon every use
            if idle_q in n.data_qubits
                push!(idle_qubits_DQC, n.inv_map[idle_q] )
            elseif idle_q in ancilla_qubits
                push!(idle_qubits_DQC, idle_q)
            else 
                throw("Communication qubits are not idling since they will be re-initialised before usage")
            end
        end
        telegate_pairs = Dict{Tuple{Int,Int}, Int}()
        telegate_qubits = Set{Int64}()
        for gate in layer
            T = typeof(gate)
            if T <: AbstractSingleQubitOperator
                qubit = affectedqubits(gate)[1]
                # In the verifciation circuit, Hadamard gates are exclusively applied to ancillas
                push!(circuit, sHadamard(qubit))
                gate_counts[1] += 1
                add_noise(circuit, [qubit], noise.p_single) # ancilla qubits experience the same single-qubit gate error probability `p_single`
            elseif T <: AbstractTwoQubitOperator
                control = affectedqubits(gate)[1]
                target = affectedqubits(gate)[2]
                DQC_control = -1
                DQC_target = -1
                control_register = -1
                target_register = -1
                if control in ancilla_qubits 
                    DQC_control = control
                    control_register = ancilla_map[control-n.num_data_and_comm_qubits] 
                else
                    DQC_control = n.inv_map[control]
                    control_register = n.register_lookup_array[DQC_control] 
                end
                if target in ancilla_qubits
                    DQC_target = target
                    target_register = ancilla_map[target-n.num_data_and_comm_qubits] 
                else
                    DQC_target = n.inv_map[target]
                    target_register = n.register_lookup_array[DQC_target] 
                end
                @assert DQC_control > 0 
                @assert DQC_target > 0 
                @assert control_register > 0 
                @assert target_register > 0 
                if control_register == target_register 
                    push!(circuit, sCNOT(DQC_control, DQC_target))
                    gate_counts[2] += 1
                    add_noise(circuit, [DQC_control, DQC_target], noise.p; two_qubits = true) # two-qubit gate error probability `p`
                else # telegate
                    pair_key = minmax(control_register, target_register)
                    prev_consecutive_telegates_pair_key = get(telegate_pairs, pair_key, 0)
                    prev_max_telegate_layers = isempty(telegate_pairs) ? 0 : maximum(values(telegate_pairs))
                    # noise before the telegate application, due to other telegates using the same communication qubit pair
                    add_noise(circuit, [DQC_control, DQC_target], 1-(1-noise.p_idle_telegate_layer)^prev_consecutive_telegates_pair_key )
                    telegate_pairs[pair_key] = prev_consecutive_telegates_pair_key + 1
                    circuit = add_telegate(circuit, DQC_control, DQC_target, control_register, target_register, n, noise)
                    gate_counts[3] += 1
                    if prev_consecutive_telegates_pair_key + 1 > prev_max_telegate_layers # the updated count in `telegate_pairs[pair_key]` is the new maximum
                        add_noise(circuit, collect(telegate_qubits), noise.p_idle_telegate_layer) # post-apply one layer of `p_idle_telegate_layer` noise on all other `telegate_qubits` encountered so far
                    else
                        # post-apply `maximum(values(telegate_pairs))-(prev_consecutive_telegates_pair_key+1)` layers of `p_idle_telegate_layer` noise on the current `DQC_control` and `DQC_target`
                        @assert prev_max_telegate_layers == maximum(values(telegate_pairs))
                        add_noise(circuit, [DQC_control, DQC_target], 1-(1-noise.p_idle_telegate_layer)^(prev_max_telegate_layers-prev_consecutive_telegates_pair_key-1)) 
                    end
                    push!(telegate_qubits, DQC_control)
                    push!(telegate_qubits, DQC_target)
                end
            elseif T <: AbstractMeasurement
                add_noise(circuit, affectedqubits(gate), noise.p) # measurement noise error probability `p` 
                push!(circuit, gate) # only ancillas are ever measured, so no DQC mapping required
                num_meas += 1
            else
                throw("Circuit contains gates that have not been classified as Single- or Two-Qubit gate so far.")
            end
        end
        num_telegates_in_layer = isempty(telegate_pairs) ? 0 : maximum(values(telegate_pairs))
        if num_telegates_in_layer == 0
            add_noise(circuit, idle_qubits_DQC, noise.p_idle) 
            depth[1] +=1
        else
            # `telegate_qubits` are excluded from this idling noise application, since idling noise has been applied to them dynamically already
            add_noise(circuit, setdiff(vcat(1:n.num_data_qubits, n.num_data_and_comm_qubits+1:length(all_qubits)), Set(telegate_qubits)), 1-(1-noise.p_idle_telegate_layer)^num_telegates_in_layer) 
            depth[2] += num_telegates_in_layer 
        end
    end
    for (i,j) in n.mapping_transpositions # revert back to original permutation
        push!(circuit, sSWAP(i, j)) 
    end
    return circuit, depth, encoding_circ_gate_counts, gate_counts, num_meas
end


"""
    build_layers(gates::Vector{AbstractOperation}, num_qubits::Int)::Vector{Vector{AbstractOperation}}

Build layers for as-soon-as-possible (ASAP) scheduling of quantum operations.

### Input

`gates` -- vector of quantum operations to be partitioned/scheduled in layers
`num_qubits` -- number of qubits on which the `gates` act

### Output

Returns a vector of layers that capture the scheduling of quantum operations in ASAP layers.
"""
function build_layers(gates::Vector{AbstractOperation}, num_qubits::Int)::Vector{Vector{AbstractOperation}}
    gate_list = copy(gates)
    layers = Vector{Vector{AbstractOperation}}()
    while !isempty(gate_list)
        idx = 1
        layer = Vector{AbstractOperation}()
        qubit_used_in_layer = falses(num_qubits) 
        del_gates = Vector{Int}()
        while idx <= length(gate_list)
            gate = gate_list[idx]
            if !( any(qubit_used_in_layer[affectedqubits(gate)]) ) 
                push!(layer, gate)
                push!(del_gates, idx)
            end
            qubit_used_in_layer[affectedqubits(gate)] .= true    
            idx += 1
        end
        push!(layers, layer)
        deleteat!(gate_list, del_gates)
    end
    return layers
end


"""
    add_telegate(circuit::Vector{AbstractOperation}, DQC_control::Int, DQC_target::Int,
                        control_register::Int, target_register::Int, n::NetworkSpecifications, noise::NoiseSpecs)::Vector{AbstractOperation}

Add a telegate based on the EJPP protocol (https://doi.org/10.1103/PhysRevA.62.052317) to the quantum circuit.

### Input

- `circuit` -- vector encoding the quantum circuit to which the telegate is to be added
- `DQC_control` -- DQC-mapped index of the control qubit for the telegate
- `DQC_target` -- DQC-mapped index of the target qubit for the telegate
- `n` -- network specificatoins
- `noise` -- noise model specifying the number of samples and circuit-level noise error rates as well as Bell pair initialisation error probability

### Output

Returns the DQC `circuit` with the appended telegate.

### Notes

When performing the EJPP protocol, we make use of communication qubits enabling the communication between `control_register` and `target_register`.
The index of a communication qubit in register `r` is determined as 
    
            `num_data_quits + num_comm_qubits_per_register*"number of registers preceding register `r`" + "offset within `r`" `.

To identify the indices in the classical register at which to store measurement results, we compute subtract the number of data qubits from
the index of the communication qubit.

When simulating probabilistic Bell pair creation, we abstract away the physical processes that lead to entanglement generation (such as photons 
being emitted and interacting with a Beam splitter), and instead prepare a perfect Bell state with `[sHadamard(DQC_control), sCNOT(DQC_control, DQC_target)]`, 
followed by a two-qubit correlated depolarising noise channel. Similar to the discussion in `construct_DQC_executable_circuit`, out of the 15 non-trivial
2-qubit Pauli strings, only 12 are harmful to the state (`XX` and `ZZ` stabilise the state), such that each two-qubit depolarising noise channel yields 
an efftive error probability of `12/15*p_bell`. 

After the completion of a EJPP protocol, we assume that the qubits are reset to the physical |0> state again. We deem this perfect initialisation a
reasonable simplification since any pre-entanglement initialisation error can be assumed to be absorbed by the Bell pair initialisation error `p_bell`.

For two-qubit gates between data and communication qubits (which potentially is an inter-species operation on a trapped-ion device, for examply), we
assume the standard two-qubit gate noise error probability `p`.
"""
function add_telegate(circuit::Vector{AbstractOperation}, DQC_control::Int, DQC_target::Int,
                        control_register::Int, target_register::Int, n::NetworkSpecifications, noise::NoiseSpecs)::Vector{AbstractOperation}
    if control_register == target_register
        throw("oops, you called `add_telegate` on a qubits in the same register")
    end
    control_comm_index = n.num_data_qubits+ (n.num_comm_qubits_per_register * (control_register-1)) + (control_register < target_register ? target_register-1 : target_register ) 
    target_comm_index = n.num_data_qubits+  (n.num_comm_qubits_per_register * (target_register-1)) + (target_register < control_register ? control_register-1 : control_register )
    # -------- EJPP Protocol --------
    # ---- I. Bell pair creation ----
    push!(circuit, sHadamard(control_comm_index))
    push!(circuit, sCNOT(control_comm_index, target_comm_index))
    add_noise(circuit, [control_comm_index, target_comm_index], noise.p_bell; two_qubits = true) # Bell state initialisation noise `p_bell`
    add_noise(circuit, [DQC_control], noise.p_idle_telegate_layer) # idling noise with probability `p_idle_telegate_layer` on data qubits partaking in telegate
    add_noise(circuit, [DQC_target], noise.p_idle_telegate_layer)  
    # ---- II. CNOT(control, comm_c) + CNOT(comm_t, target) ----
    push!(circuit, sCNOT(DQC_control, control_comm_index))
    push!(circuit, sCNOT(target_comm_index, DQC_target))
    add_noise(circuit, [DQC_control, control_comm_index], noise.p; two_qubits = true) # two-qubit gate noise with error probability `p`
    add_noise(circuit, [target_comm_index, DQC_target], noise.p; two_qubits = true) 
    # ---- III. comm_c + comm_t Measurement ---- 
    classical_register_index_control = control_comm_index - n.num_data_qubits 
    meas_control = sMRZ(control_comm_index, classical_register_index_control ) 
    add_noise(circuit, [control_comm_index], noise.p) # communication qubit measurement noise with error probability `p` 
    push!(circuit, meas_control)
    classical_register_index_target = target_comm_index - n.num_data_qubits
    meas_target = sMRZ(target_comm_index, classical_register_index_target)
    push!(circuit, sHadamard(target_comm_index)) # measuring in the X-basis amounts to applying `sHadamard` and measuring in the Z-basis
    add_noise(circuit, [target_comm_index], noise.p_single) 
    add_noise(circuit, [target_comm_index], noise.p) # communication qubit measurement noise with error probability `p`
    push!(circuit, meas_target)
    add_noise(circuit, [DQC_control], noise.p_idle) # idling depolarising_noise with error probability `p_idle`
    add_noise(circuit, [DQC_target], noise.p_idle) 
    # ---- IV. Conditional Operations on data qubits ----
    push!(circuit, ConditionalGate(sX(DQC_target),sId1(DQC_target), meas_control.bit)) # perform operation conditioned on the measurement bit
    push!(circuit, ConditionalGate(sZ(DQC_control),sId1(DQC_control), meas_target.bit))
    add_noise(circuit, [DQC_control], noise.p_single) 
    add_noise(circuit, [DQC_target], noise.p_single) 
    return circuit
end



# ------------ Noise functions for executable circuit ---------------

# Depolarising channel: https://github.com/QuantumSavory/QuantumClifford.jl/blob/74ee758e87f5d7b1255d6747b346cff15ee10cea/src/noise.jl#L63-73
"""
    add_noise(circuit::Vector{AbstractOperation}, qubits::Vector{Int}, prob::Float64; two_qubits = false)::Vector{AbstractOperation}

Add (correlated or uncorrelated) depolarising noise with a specified error probability `prob` on a set of `qubits`.

### Input

- `circuit` -- vector encoding the quantum circuit to which the noise is to be added
- `qubits` -- qubits onto which the noise channel is applied
- `prob` -- depolarising noise channel error proability
- `two_qubits` -- (optional, default: `false`) controls whether or not to apply correlated noise on two qubits

### Output

Returns the DQC `circuit` with the appended noise operation.

### Notes

The error channel `UnbiasedUncorrelatedNoise` applies a single-qubit Pauli gate (`X`, `Y` or `Z`) on each qubit in `qubits` with probability `prob/3`.

The error channel `TwoQubitDepolarisingNoise` applies a non-trivial two-qubit Pauli gate (`IX`, `IY`, ..., `ZY`, `ZZ`) on the two specificied 
qubits in `qubits` with probability `prob/15`.
"""
function add_noise(circuit::Vector{AbstractOperation}, qubits::Vector{Int}, prob::Float64; two_qubits = false)::Vector{AbstractOperation}
    if prob<0 || prob > 1
        throw("Please provide a valid noise probabilty in [0,1]")
    end
    if two_qubits
        @assert length(qubits) == 2 "you are trying to apply a two-qubit correlated noise channel to a system of size $(length(qubits))"
        noise = NoiseOp(TwoQubitDepolarisingNoise(prob), qubits) # apply correlated noise
    else
        if !isempty(qubits)
            noise = NoiseOp(UnbiasedUncorrelatedNoise(prob), qubits) # apply uncorrelated noise
        else
            return circuit # no qubits suffer from noise 
        end
    end
    push!(circuit, noise)
    return circuit
end


"""
    TwoQubitDepolarisingNoise

Type that defines a two-qubit correlated depolarising noise channel.

### Fields

- `p` -- error probability for depolarising noise channel
"""
struct TwoQubitDepolarisingNoise{T} <: QuantumClifford.AbstractNoise
    p::T
end

TwoQubitDepolarisingNoise(p::Integer) = TwoQubitDepolarisingNoise(float(p))


"""
    applynoise!(s::AbstractStabilizer, noise::TwoQubitDepolarisingNoise, indices::Tuple{Int, Int})::AbstractStabilizer

Implement the circuit execution structure of two-qubit correlated depolarising noise of error
probability `p`.

### Input

- `s` -- stabiliser object on which the noise channel is performed during circuit execution
- `noise` -- noise object of type `TwoQubitDepolarisingNoise` capturing the error probability
- `indices` -- qubit indices on which the correlated noise is applied

### Output

Returns the stabiliser object after the application of the noise channel.
"""
function applynoise!(s::AbstractStabilizer, noise::TwoQubitDepolarisingNoise, indices::Tuple{Int, Int})::AbstractStabilizer
    infid = noise.p/15
    i,j = indices[1], indices[2]
    r = rand()
    if r<infid
        apply_single_x!(s,i)
    elseif r<2infid
        apply_single_z!(s,i)
    elseif r<3infid
        apply_single_y!(s,i)
    elseif r<4infid
        apply_single_x!(s,j)
    elseif r<5infid
        apply_single_z!(s,j)
    elseif r<6infid
        apply_single_y!(s,j)
    elseif r<7infid
        apply_single_x!(s,i)
        apply_single_x!(s,j)
    elseif r<8infid
        apply_single_x!(s,i)
        apply_single_z!(s,j)
    elseif r<9infid
        apply_single_x!(s,i)
        apply_single_y!(s,j)
    elseif r<10infid
        apply_single_z!(s,i)
        apply_single_x!(s,j)
    elseif r<11infid
        apply_single_z!(s,i)
        apply_single_y!(s,j)
    elseif r<12infid
        apply_single_z!(s,i)
        apply_single_z!(s,j)
    elseif r<13infid
        apply_single_y!(s,i)
        apply_single_x!(s,j)
    elseif r<14infid
        apply_single_y!(s,i)
        apply_single_y!(s,j)
    elseif r<15infid
        apply_single_y!(s,i)
        apply_single_z!(s,j)
    end
    return s
end


function affectedqubits(op::AbstractSingleQubitOperator)::Vector{Int}
    # returns the affected qubits of a single-qubit gate as vector
    qs = Int[]
    append!(qs, op.q)
    return qs
end

function affectedqubits(op::AbstractTwoQubitOperator)::Vector{Int}
    # returns the affected qubits of a two-qubit gate as vector
    qs = Int[]
    append!(qs, op.q1)
    append!(qs, op.q2)
    return qs
end

function affectedqubits(op::ConditionalGate)::Vector{Int}
    # returns the affected qubits of a conditional gate as vector
    qs = Int[]
    append!(qs, collect(affectedqubits(op.truegate)))
    if op.falsegate !== nothing
       append!(qs, collect(affectedqubits(op.falsegate)))
    end
    return unique(qs)
end

function affectedqubits(op::sMZ)::Vector{Int}
    # returns the affected qubits of a measurement operation as vector
    qs = Int[]
    append!(qs, op.qubit)
    return qs
end




end

