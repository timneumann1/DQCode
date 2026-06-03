# dqc_simulator.jl

"""
Functions that construct DQC-executable circuits from an optimised encoding circuit 
for a logical zero state, simulate the DQC execution including telegates under circuit-level
and Bell pair initialisation noise, and collect corresponding statistics, e.g. on logical error rates.

Credit: 
- The pipeline for logical evaluation (`dqc_logical_evaluation`) is largely based on the implementation
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

import QuantumClifford: apply!, affectedqubits, applynoise! 

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
We first check for validity of the noiseless circuit, also with resepct to the DQC evaluation, 
and then perform a DQC evaluation.

The effective idle error probability during a telegate is computed as `1-(1-p_idle)^telegate_idle_depth`.
We consider the measurement, intitialisation and two-qubit gate error probability to be equal to `p`,
and disregard the effect of crosstalk.
"""
function dqc_non_ft_encoding_simulation(num_samples::Int, ps::Vector{Float64}, p_bells::Vector{Float64}, 
                                        telegate_idle_depth::Int, p_single_ratio::Float64, p_idle_ratio::Float64,
                                        code_params::CodeParameters, network_specs::NetworkSpecifications, 
                                        circuit::Vector{AbstractOperation})::Tuple{Any, Vector{AbstractOperation}, Vector{AbstractOperation}, Vector{AbstractOperation}}
    data_circuit = deepcopy(circuit)
    #------------ Noiseless verification -------------
    @info "Noiseless testing of raw encoding circuit ..."
    noise_verif = NoiseSpecs(1e2,0,0,0,0,0) # performing a small number of runs without noise
    quantum_clifford_verification_circ = Vector{AbstractOperation}() # no verification circuit → non-FT
    DQC_circuit,_,_,_,_ = construct_DQC_executable_circuit(data_circuit, quantum_clifford_verification_circ, 0, [], network_specs, noise_verif)
    num_ancillas = 0
    ancilla_map = Vector{Int}()
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
    ) = dqc_logical_evaluation(data_circuit, quantum_clifford_verification_circ, 0, [],
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
        push!(data, (p=p, p_bell=p_bell, logical_error_rate=logical_error_rate, acceptance_ratio=acceptance_ratio, 
                    discarded_runs=discarded_runs, n_samples=n_samples, logical_failures=logical_failures,
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
                                    circuit::Vector{AbstractOperation}, method::String)::Tuple{Any, Vector{AbstractOperation}, Vector{AbstractOperation}, 
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
    # our initial definition of logical still applies and can be used for the noiseless extraction of logical X error rate
We then parse through this, doing a lot of indexing acrobatics to finally have the encoding circuit in our desired format, where
data qubits come first, then communication qubits, then ancilla qubits (and lastly ancillas for the noisefree syndrome extraction)

We first check for validity of the noiseless circuit, also with resepct to the DQC evaluation, 
and then perform a DQC evaluation

The effective idle error probability during a telegate is computed as `1-(1-p_idle)^telegate_idle_depth`.
We consider the measurement, intitialisation and two-qubit gate error probability to be equal to `p`,
and disregard the effect of crosstalk.
"""
function dqc_ft_encoding_simulation(num_samples::Int, ps::Vector{Float64}, p_bells::Vector{Float64}, 
                                    telegate_idle_depth::Int, p_single_ratio::Float64, p_idle_ratio::Float64,
                                    code_params::CodeParameters, network_specs::NetworkSpecifications, mqt_path::String,
                                    circuit::Vector{AbstractOperation}, method::String)::Tuple{Any, Vector{AbstractOperation}, Vector{AbstractOperation}, 
                                                                                                Vector{AbstractOperation}, Int, Int, Int, Vector{Int}}
    data_circuit = deepcopy(circuit)
    qasm = qc_circuit_to_qasm(circuit)
    python_bin = joinpath(mqt_path, ".venv/bin/python3")
    script_path =  joinpath(mqt_path, "scripts/verification_circuit.py")
    # --------- Verification Circuit --------------
    @info "Retrieving verification circuit"
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
    #------------ Noiseless verification -------------
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
        push!(data, (p=p, p_bell=p_bell, logical_error_rate=logical_error_rate, acceptance_ratio=acceptance_ratio, 
                    discarded_runs=discarded_runs, n_samples=n_samples, logical_failures=logical_failures,
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
`avg_fidelity`, `acceptance_ratio`, etc.

### Notes

In this function, we first retrieve the DQC-executable circuit `DQC_circuit` from `construct_DQC_executable_circuit`. Then,
we setup a noisefree QEC cycle / syndrome decoding routine.  For the noise-free decoding, we use a CSS Lookup Table Decoder. 

Simulating a number of Monte Carlo trajectories, we sample from the noise distribution and gather statistics about the quality of FT encoding.
Here, we first discard runs according to `verification_bits`, then identify pre-decoding X- and Z- errors, and then determine logical X errors, 
i.e., errors in the logical Z observables (note that for the logical Z state, logical Z errors cannot occur).
Therefore, we collect the `error_guess` by providing the decoder with the measured `syndrome` (`error_guess` collects `n`` guesses for X-errors
and `n` guesses for Z errors -- we are interested in the first `n` guesses, i.e., whether the decoder predicts that a certain physical X error happened).
If the decoder cannot infer an error guess, we register a logical error on all logical qubits, and a fidelity of `0.0` for this run (since the
resulting state lives outside of the codespace and cannot be corrected). 

If the decoder can make an error guess, we expose two different methods as proxies for the logical initialisation error rate:

    (i) We noiselessly measure the logical Z observables and compute the logical error rate based on the `measured_logical_Z_bits` as well as 
        `faults_matrix_z` (a `k \times n` matrix), which collects information about which (of the `n`) physical errors flip which (of the `k`) logical
        Z operators / logical Z observables via anti-commutation. The logic here is as follows: 
        If `sum_mod += faults_matrix_z[j, q] * error_guess[q]` is 1, then there is an odd number of indices that are jointly supported by both `faults_matrix_z`
        and `error_guess` for the given `j`th logical Z operator. In this case, correcting based on `error_guess` will flip the `k`th logical Z operator. This
        is only desirable if there actually was a logical X error on the `k`th logical Z observable, i.e., `sum_mod = measured_logical_Z_bits[j]`. Otherwise,
        a logical failure is recorded.
        If `sum_mod += faults_matrix_z[j, q] * error_guess[q]` is 0, then there is an even number of indices that are jointly supported by both `faults_matrix_z`
        and `error_guess` for the given `j`th logical Z operator. In this case, correcting based on `error_guess` will not flip the `k`th logical Z operator by 
        error degeneracy of the code (applying an even number of operations that anti-commute with the logical operator overall cancels out). This is desirable as 
        long as there was actually no logical X error on the `k`th logical Z observable, i.e., `sum_mod = measured_logical_Z_bits[j] again. Otherwise,
        a logical failure is recorded.
        In the end, the logical error rate is derived as the mean of all `k` logical Z observable error rates, conditioned on the accepted runs.
            
    (ii) Alternatively, we can also apply the correction gate inferred by the LUT decoder, and then measure the fidelity of the resulting state with our target state.
        Therefore, we apply the correction gate (https://github.com/QuantumSavory/QuantumClifford.jl/blob/444f341a50d2926b16b63d98586b8b06a7b6ac10/src/ecc/decoder_correction_gate.jl)
        that maps the quantum state back to the codespace (since there is a correctable syndrome, we know that the inferred correction maps back to the codespace) and 
        call `dot` between the corrected `corrected_state_data_qubits` state and the target state `code_params.target_state`. We update fidelities per sample in order to
        avoid having to save the fidelity for each of the samples. This is done using the update equation

            `avg_fidelity = (avg_fidelity*(samples-1)+x)/samples = avg_fidelity+(x-avg_fidelity)/samples → avg_fidelity += (x-avg_fidelity)/samples`,
        
        where `samples` is replaced with `sample-discarded_runs`.

In the entire procedure, there are two main sources of uncertainty: the probabilistic noise sampling of errors in the circuit, and the measurement of the logical Z operators. 
In a real-life experiment, both sources of errors exist (where probabilistically measuring the correct state implies that the system indeed is in the correct state afterwards),
and method (i) mimics this closely. With method (ii), we eliminate the second source of uncertainty by determining a precise fidelity. 

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
    H = parity_checks(css_lut_decoder)
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
    @assert setdiff(Z_type_indices, X_type_indices) == Z_type_indices # there were mixed Pauli strings amongst the stabiliser generators
    # determine `faults_matrix_z` -- the fault matrix is a (2k)x(2n) dimensional matrix → last k rows specify logical Z part
    faults_matrix_z = css_lut_decoder.faults_matrix[end÷2+1:end,:] 
    k = size(faults_matrix_z, 1) 
    n = size(faults_matrix_z, 2)
    @assert k == code_params.k
    @assert n == 2*code_params.n
    noisefree_syndrome_circ, num_noisefree_syndrome_ancillas, noisefree_syndrome_bits = syndrome_circuit(H, network_specs.num_data_and_comm_qubits + num_ancillas + 1, 
                                                                                                            network_specs.num_comm_qubits + num_ancillas + 1, network_specs)
    noisefree_logical_Z_circ, num_noisefree_logical_Z_ancillas, noisefree_logical_Z_bits = syndrome_circuit(code_params.logical_Zs,
                                                                                                            network_specs.num_data_and_comm_qubits + num_ancillas + num_noisefree_syndrome_ancillas + 1, 
                                                                                                            last(noisefree_syndrome_bits)+1, network_specs )
    total_number_qubits = network_specs.num_data_and_comm_qubits + num_ancillas +  num_noisefree_syndrome_ancillas + num_noisefree_logical_Z_ancillas
    total_number_classical_regs = network_specs.num_comm_qubits + num_ancillas + num_noisefree_syndrome_ancillas + num_noisefree_logical_Z_ancillas
    circuit = vcat(DQC_circuit, noisefree_syndrome_circ)
    discarded_runs = 0
    z_errors_pre_decoding = 0
    x_errors_pre_decoding = 0
    logical_failures = zeros(code_params.k) # track logical Z failures per qubit
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
        error_guess = decode(css_lut_decoder, syndrome) # retrieve the error guess from the LUT decoder 
        if isnothing(error_guess) # if no error guess can be extracted from the Lookup table, we collect a logical error on all logical qubits and fidelity 0.0
            logical_failures .+= 1
            avg_fidelity -= avg_fidelity / (sample-discarded_runs) # if execution reaches this line, we know that sample \neq discarded_runs holds
            continue
        end
        #-------- (i) Analytical logical error rate computation --------
        for j in 1:k # iterate over the k logical Z operators
            sum_mod = 0
            @inbounds @simd for q in 1:n # iterate over all the physical qubits/error locations
                sum_mod += faults_matrix_z[j, q] * error_guess[q] 
            end
            sum_mod %= 2 
            if sum_mod != measured_logical_Z_bits[j] # `measured_logical_Z_bits` captures the measured Z-faults
                logical_failures[j] += 1
            end
        end
        #-------- (ii) Simulated correction gate --------
        @info noisefree_syndrome_bits
        correction_gate = DecoderCorrectionGate(css_lut_decoder, network_specs.data_qubits, noisefree_syndrome_bits ) 
        corrected_state,_ = mctrajectory!(copy(state_post_syndrome),[correction_gate]) 
        corrected_state_data_qubits = stabilizerview( ptrace(copy(corrected_state.stab), collect(network_specs.num_data_qubits+1:total_number_qubits)) )
        post_decoding_fidelity = dot(corrected_state_data_qubits, code_params.target_state)
        avg_fidelity += (post_decoding_fidelity-avg_fidelity)/(sample-discarded_runs)
        if maximum(logical_failures) >= 1000
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



# ---------------------------------------------------------------
# ------------ Construction of executable circuit ---------------
# ---------------------------------------------------------------

function construct_DQC_executable_circuit(data_circuit, verification_circuit, num_ancillas, ancilla_map, n::NetworkSpecifications, noise::NoiseSpecs)
    
    # Map, then add verification on mapped, then unmap

    circuit = Vector{AbstractOperation}()
   
    # in the permutation, [1,9,...] indicates that the 9th element gets permuted to second position, "9 is mapped to 2"
    # Apply the inverse permutation of the mapping by applying transpoistions of inverse perm in left action <-> transpoistions of perm in right action via reverse(mapping) [the mapping contains transposition derived from the permutation, implementing it in left action]
    for (i,j) in reverse(n.mapping_transpositions)
        push!(circuit, sSWAP(i,j))  # We could also use comm_perm_idx or comm_inv_perm_idx, since the relabeling based on the permutation conjugtes and thus fixes the permutation induces by the transposition SWAPS
    end
        
    #go through circuit and accumulate ASAP layers (we may assume that no further grouping for telegates can be done)
    layers_enc_circ = build_layers(data_circuit, n.num_data_qubits)

    # Add initialisatin noise to all qubits
    # Depolarising channel: https://github.com/QuantumSavory/QuantumClifford.jl/blob/74ee758e87f5d7b1255d6747b346cff15ee10cea/src/noise.jl#L63-73
    add_noise(circuit, [n.inv_map[data_q] for data_q in collect(1:n.num_data_qubits)], noise.p) # initialisation noise p

    depth = [0,0] # first index encodes the number of ordinary layers, second index encodes the number of telegate layers (i.e., non-parallel telegates)
    gate_counts = [0,0,0]
    num_meas = 0
    # number of gates will be counted separately (see resource_estimate > helper)
    for layer in layers_enc_circ
    # for every layer, add gate noise for every normal gate and idling noise for any idle qubits (if there is telegates in the layer, we increase the probability of the noise channel)
    # also, insert the telegate gadgets in place (this includes the folloing noise: comm init noise, two qubit depolarising for each gate , measeuremtn noise, classical noise, one more gate noise for single quits )
    
        #find idling gates and apply idle_depolarising_noise
        affected_qubits = Set( Iterators.flatten( [affectedqubits(gate) for gate in layer] ) ) 
        idle_qubits = setdiff(1:n.num_data_qubits, affected_qubits)
        idle_qubits_DQC = [n.inv_map[idle_q] for idle_q in idle_qubits]
    
        #telegates_layer = false
        num_telegates_layer = 0
        telegate_pairs = Set{Tuple{Int,Int}}()
        telegate_qubits = Set{Int64}()

        for gate in layer
            T = typeof(gate)
            if T <: AbstractSingleQubitOperator
                qubit = affectedqubits(gate)[1]
                DQC_qubit = n.inv_map[qubit]
                push!(circuit, sHadamard(DQC_qubit))
                gate_counts[1] += 1
                add_noise(circuit, [DQC_qubit], noise.p_single) # Single-qubit noise
                
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
                    add_noise(circuit, [DQC_control, DQC_target], noise.p; two_qubits = true) # two-qubit noise
                else
                    # Perform telegate between control and target qubit in different registers
                    if num_telegates_layer == 0
                        num_telegates_layer += 1 # if this is the first telegate in the layer, we need to increase the telegate count; for any of the following telegates, we increase the number of layers only when the comm pair has already been used (since then we must wait)
                    end
                    #telegates_layer = true
                    if (control_register, target_register) ∈ telegate_pairs || (target_register, control_register) ∈ telegate_pairs # if the communication qubits for this register pair have been used, we need to increase the number of telegates for this layer, and apply as much noise as there have been previous telegates in this layer
                        #for _ in 1:num_telegates_layer
                        add_noise(circuit, [DQC_control, DQC_target], 1-(1-noise.p_idle_telegate_layer)^num_telegates_layer ) 
                        add_noise(circuit, collect(telegate_qubits), noise.p_idle_telegate_layer)
                        #end
                        num_telegates_layer +=1
                    end 
                    # if the telegate was not the first one in the layer, or the comm qubits had not been used, then we don't need to increase the num_telegates_layer count
                    circuit = add_telegate(circuit, DQC_control, DQC_target, control_register, target_register, n, noise)
                    gate_counts[3] += 1
                    # due to this telegate, previous telegate qubits incur noise (all non-telegate qubits are handled later)
                    # the current telegate qubits are not in the telegate_qubits set yet by construction, which is why we add the noise above separately.
                    
                    push!(telegate_pairs, (control_register, target_register))
                    push!(telegate_qubits, DQC_control)
                    push!(telegate_qubits, DQC_target)
                end
            else
                throw("Circuit contains gates that have not been classified as Single- or Two-Qubit gate so far.")
            end
        end
        if num_telegates_layer == 0
            add_noise(circuit, idle_qubits_DQC, noise.p_idle) # idling noise (if telegate: to all non-telegate qubits, no matter if idle or not (since telegate is much longer), otherwise: idling noise p on idle qubits)
            depth[1] +=1
        else
            #for _ in 1:num_telegates_layer 
            add_noise(circuit, setdiff(1:n.num_data_qubits, Set(telegate_qubits)), 1-(1-noise.p_idle_telegate_layer)^num_telegates_layer)  # the telegate qubits are excluded from this noise application on the genuine idle qubits
            #end
            depth[2] += num_telegates_layer 
            #@info "Number of telegates in this layer: $layer is $num_telegates_layer"
        end
    end

    #@info "Gate counts encoding circuit: $gate_counts"
    # Add the verification circuit before reversing the virtual mapping
    # telegates between ancillas and data qubits can use the comm qubits of the register, likewise for telegates between ancillas and flags

    encoding_circ_gate_counts = copy(gate_counts)

    layers_ver_circ = build_layers(verification_circuit, n.num_data_and_comm_qubits + num_ancillas)

    #@info "Mapping: $ancilla_map"

    all_qubits = 1:(n.num_data_and_comm_qubits + num_ancillas)
    ancilla_qubits = setdiff(all_qubits, 1:n.num_data_and_comm_qubits)

    #add_noise(circuit, [data_q for data_q in collect(1:n.num_data_qubits)], noise.idle_depolarising_noise) # since noise is applied to all qubits, we don't need to worry about mapping
    add_noise(circuit, ancilla_qubits, noise.p) #  init noise on ancilla qubits for verification measurements

    for layer in layers_ver_circ
        # analogous to raw encoding circuit
        affected_qubits = Set( Iterators.flatten( [affectedqubits(gate) for gate in layer] ) ) 
        idle_qubits = setdiff(all_qubits, union(Set(n.comm_qubits), affected_qubits))
        #println("all qubits: $all_qubits, idle: $idle_qubits, data: $(n.data_qubits)")
        idle_qubits_DQC = Vector{Int}()
        for idle_q in idle_qubits
            # comm qubits will never be treated as idle, since they get reinitialised every time
            if idle_q in n.data_qubits
                push!(idle_qubits_DQC, n.inv_map[idle_q] )
            elseif idle_q in ancilla_qubits
                push!(idle_qubits_DQC, idle_q)
            else 
                throw("Communication qubits are not idling since they will be re-initialised before usage")
            end
        end

        #telegates_layer = false
        num_telegates_layer = 0
        #tele_qubits = Vector{Int64}()
        telegate_pairs = Set{Tuple{Int,Int}}()
        telegate_qubits = Set{Int64}()
        for gate in layer
            T = typeof(gate)
            if T <: AbstractSingleQubitOperator
                qubit = affectedqubits(gate)[1]
                # In the verifciation circuit, Hadamard gates are ONLY applied to ancillas, which sit at their regular index
                push!(circuit, sHadamard(qubit))
                gate_counts[1] += 1
                add_noise(circuit, [qubit], noise.p_single) # single-qubit noise, ancilla qubits experience the same sort of noise, since they are of the same physical type
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
                    add_noise(circuit, [DQC_control, DQC_target], noise.p; two_qubits = true) # two-qubit noise
                else
                    if num_telegates_layer == 0
                        num_telegates_layer += 1 # if this is the first telegate in the layer, we need to increase the telegate count; for any of the following telegates, we increase the number of layers only when the comm pair has already been used (since then we must wait)
                    end
                    #telegates_layer = true
                    if (control_register, target_register) ∈ telegate_pairs || (target_register, control_register) ∈ telegate_pairs # if the communication qubits for this register pair have been used, we need to increase the number of telegates for this layer, and apply as much noise as there have been previous telegates in this layer
                        #for _ in 1:num_telegates_layer
                        add_noise(circuit, [DQC_control, DQC_target], 1-(1-noise.p_idle_telegate_layer)^num_telegates_layer ) # current telegate qubit obtain noise from previous telegate_layers (where they were not yet part of the telegate qubits by construction)
                        add_noise(circuit, collect(telegate_qubits), noise.p_idle_telegate_layer) # other telegate qubits (used so far) incur noise
                        #end
                        num_telegates_layer +=1
                    end 

                    circuit = add_telegate(circuit, DQC_control, DQC_target, control_register, target_register, n, noise)
                    #telegates_layer = true
                    #num_telegates_layer +=1
                    gate_counts[3] += 1
                    push!(telegate_pairs, (control_register, target_register))
                    push!(telegate_qubits, DQC_control)
                    push!(telegate_qubits, DQC_target)
                    
                end
            elseif T <: AbstractMeasurement
                add_noise(circuit, affectedqubits(gate), noise.p) # measurement noise the ancilla is of the same type, thus we have the same measurement noise
                push!(circuit, gate) # only ancillas are ever measured, so we don't need a remapping
                num_meas += 1
            else
                throw("Circuit contains gates that have not been classified as Single- or Two-Qubit gate so far.")
            end
        end
       # println("idle qubits: $idle_qubits_DQC")
        if num_telegates_layer == 0
            add_noise(circuit, idle_qubits_DQC, noise.p_idle) # idling noise (if telegate: to all non-telegate qubits, no matter if idle or not (since telegate is much longer), otherwise: idling noise p on idle qubits)
            depth[1] +=1
        else
            # the telegate qubits are excluded from this noise application on the genuine idle qubits (here, noise is also applied on the ancillas)
            add_noise(circuit, setdiff(vcat(1:n.num_data_qubits, n.num_data_and_comm_qubits+1:length(all_qubits)), Set(telegate_qubits)), 1-(1-noise.p_idle_telegate_layer)^num_telegates_layer) 
            # the set diff targets all non-telegate qubits (since noise is uncorrelated here, DQC vs. usually indices does not matter)
            depth[2] += num_telegates_layer 
            #@info "Number of telegates in this layer: $layer is $num_telegates_layer"
        end
        #telegates_layer ? add_noise(circuit, idle_qubits_DQC, noise.p_idle_telegate_layer) : add_noise(circuit, idle_qubits_DQC, noise.p) # idling noise
#        telegates_layer ? add_noise(circuit, setdiff(vcat(1:n.num_data_qubits, n.num_data_and_comm_qubits+1:length(all_qubits)), Set(tele_qubits)), noise.p_idle_telegate_layer) : add_noise(circuit, idle_qubits_DQC, noise.p_idle) # idling noise (if telegate: to all non-telegate qubits (data + ancilla), no matter if idle or not (since telegate is much longer), otherwise: idling noise p on idle qubits)

    end

    # Revert swapping for measurement of target state
    for (i,j) in n.mapping_transpositions # reverse-reverse -> revert back to original permutation
        push!(circuit, sSWAP(i, j)) 
    end

    return circuit, depth, encoding_circ_gate_counts, gate_counts, num_meas
end


# ------------ ASAP Layers for construction of executable circuit ---------------

function build_layers(gates,num_qubits)
    # adapted from MQT circuit_utils.py file
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
    #@info "Layers: $layers"
    return layers
end




# ------------ Telegates for construction of executable circuit ---------------


function add_telegate(circuit, DQC_control, DQC_target, control_register, target_register, n, noise)
    if control_register == target_register
        throw("Ooops, this is not a telegate.")
    end
    
    # control_comm_idx is num_data_qubits + control_q_per_core*number of control_registers that came before the active register + target register offset
    control_comm_index = n.num_data_qubits+ (n.num_comm_qubits_per_register * (control_register-1)) + (control_register < target_register ? target_register-1 : target_register ) 
    target_comm_index = n.num_data_qubits+  (n.num_comm_qubits_per_register * (target_register-1)) + (target_register < control_register ? control_register-1 : control_register )
    
    # -------- EJPP Protocol --------
    
    # ---- I. Bell pair creation ----

    # This H-CNOT only mimics the way Bell state entanglement is created; in reality, this is achieved via beam splitters or such.
    # In order to account for the differing physical circumstance, we apply a two-qubit depolarising channel with a specific noise probability afterwards
    push!(circuit, sHadamard(control_comm_index))
    push!(circuit, sCNOT(control_comm_index, target_comm_index))
    add_noise(circuit, [control_comm_index, target_comm_index], noise.p_bell; two_qubits = true) # Bell state initialisation noise
    add_noise(circuit, [DQC_control], noise.p_idle_telegate_layer) # idling noise on data qubit (was excluded in the telegate layer idling noise channel)
    add_noise(circuit, [DQC_target], noise.p_idle_telegate_layer) # "

    # ---- II. CNOT(control, comm_c) + CNOT(comm_t, target) ----

    push!(circuit, sCNOT(DQC_control, control_comm_index))
    push!(circuit, sCNOT(target_comm_index, DQC_target))
    add_noise(circuit, [DQC_control, control_comm_index], noise.p; two_qubits = true) # mixed-species CNOT noise
    add_noise(circuit, [target_comm_index, DQC_target], noise.p; two_qubits = true) # "


    # ---- III. comm1 + comm2 Measurement ----
 
    #pauli_string_control = build_pauli_string_measurement(n.num_data_and_comm_qubits, [control_comm_index])
    classical_register_index_control = control_comm_index - n.num_data_qubits #sum(n.register_sizes[1:control_register])
    #meas_control = PauliMeasurement(pauli_string_control, classical_register_index_control)
    meas_control = sMRZ(control_comm_index, classical_register_index_control ) 
    # The restoration can probably be neglected since in reality, the Bell pair will be created anew via photonic beam splitters. 
    # For the sake of simulation however, we assume lossless restoration. This is handled by using sMRZ rather than sMZ.
    add_noise(circuit, [control_comm_index], noise.p) # Comm qubit measurement noise
    push!(circuit, meas_control)

    classical_register_index_target = target_comm_index - n.num_data_qubits#sum(n.register_sizes[1:target_register])
    meas_target = sMRZ(target_comm_index, classical_register_index_target)#PauliMeasurement(pauli_string_target, classical_register_index_target)
    push!(circuit, sHadamard(target_comm_index)) # measuring in the X-basis requires us to Hadamard-transform and then measure in the Z basis
    add_noise(circuit, [target_comm_index], noise.p_single) # single-gate noise
    add_noise(circuit, [target_comm_index], noise.p) # Comm qubit measurement noise
    push!(circuit, meas_target)

    add_noise(circuit, [DQC_control], noise.p_idle) # usual idle depolarising_noise
    add_noise(circuit, [DQC_target], noise.p_idle) # "

    # ---- IV. Conditional Operations on data qubits ----
    # perform conditional operations, conditioned on the measurement bit: if state.bits[op.controlbit] is true, the measurment yielded eigenvalue -1, if it is false, it yielded +1
    push!(circuit, ConditionalGate(sX(DQC_target),sId1(DQC_target), meas_control.bit))  
    push!(circuit, ConditionalGate(sZ(DQC_control),sId1(DQC_control), meas_target.bit))

    add_noise(circuit, [DQC_control], noise.p_single) # single-qubit gate noise
    add_noise(circuit, [DQC_target], noise.p_single)  # "

    # Ideally, we would add noise conditional on whether or not we apply a gate. However, it is acceptable to simply assume that we apply an identity gate
    #add_noise(circuit, [DQC_control], noise.classical_comm_noise) # no matter if X or I applied, we assume some classical communication noise
    #add_noise(circuit, [DQC_target], noise.classical_comm_noise) # 

    
    #push!(circuit, Types.ConditionalGate(sX(control_comm_index),sId1(control_comm_index), meas_control.bit)) # restore the |0> state in the control comm qubit
    #push!(circuit, Types.ConditionalGate(sX(target_comm_index),sId1(target_comm_index), meas_target.bit))  # restore the |0> state in the target comm qubit

    return circuit

end

# function build_pauli_string_measurement(num_qubits::Int, qubits::Vector{Int})
#     pauli = I # we can safely assume that the first qubit is a data qubit, since this is only false whenever there are zero qubits
#     @inbounds for i in 2:(num_qubits) # traverses all data and comm qubits 
#         pauli = (i in qubits) ? pauli⊗Z : pauli⊗I
#     end
#     return pauli
# end


# ------------ Noise functions for executable circuit ---------------

function add_noise(circuit, qubits::Vector{Int}, prob::Float64; two_qubits = false) 
    """Depolarising noise on a set of qubits"""
    if prob<0 || prob > 1
        throw("Please provide a valid noise probabilty in [0,1]")
    end
    if two_qubits
        @assert length(qubits) == 2 "Trying to apply a two-qubit channel to a system of size $(length(qubits))"
        noise = NoiseOp(TwoQubitDepolarisingNoise(prob),qubits);
    else
        noise = NoiseOp(UnbiasedUncorrelatedNoise(prob),qubits);
    end

    push!(circuit, noise)
    return circuit
end


struct TwoQubitDepolarisingNoise{T} <: QuantumClifford.AbstractNoise
    p::T
end
TwoQubitDepolarisingNoise(p::Integer) = TwoQubitDepolarisingNoise(float(p))

function applynoise!(s::AbstractStabilizer, noise::TwoQubitDepolarisingNoise, indices::Tuple{Int, Int})
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
    s
end


#We can also use NoisyGate: https://github.com/QuantumSavory/QuantumClifford.jl/blob/74ee758e87f5d7b1255d6747b346cff15ee10cea/docs/src/noisycircuits_ops.md

# function apply!(state::Register, op::ConditionalGate)
#     if state.bits[op.controlbit]
#         apply!(state, op.truegate)
#     else
#         apply!(state, op.falsegate)
#     end
#     return state
# end

function affectedqubits(op::ConditionalGate)
    qs = Int[]
    append!(qs, collect(affectedqubits(op.truegate)))
    if op.falsegate !== nothing
       append!(qs, collect(affectedqubits(op.falsegate)))
    end
    return unique(qs)
end

function affectedqubits(op::AbstractSingleQubitOperator)
    qs = Int[]
    append!(qs, op.q)
    return qs
end

function affectedqubits(op::sMZ)
    qs = Int[]
    append!(qs, op.qubit)
    return qs
end

function affectedqubits(op::AbstractTwoQubitOperator)
    qs = Int[]
    append!(qs, op.q1)
    append!(qs, op.q2)
    return qs
end



end

