# baseline_encoding.jl

"""
Encoding circuit synthesis baselines as defined in
- Qiskit (https://github.com/Qiskit/qiskit/blob/main/qiskit/synthesis/clifford/clifford_decompose_full.py) and
- Munich Quantum Toolkit QECC Version 1.0.9. (https://arxiv.org/abs/2408.11894) as implemented in 
    https://github.com/munich-quantum-toolkit/qecc/tree/main/src/mqt/qecc/
"""
module BaselineEncoding

export run_qiskit_baseline, run_mqt_baseline

using ..Types
using ..Helper

using QuantumClifford
using Serialization
using PyCall
using JSON
using CSV, DataFrames

np = pyimport("numpy")
qiskit = pyimport("qiskit")
qi = pyimport("qiskit.quantum_info")
qasm2 = pyimport("qiskit.qasm2")
synth = pyimport("qiskit.synthesis")
plt = pyimport("matplotlib.pyplot")
Clifford = qi.Clifford
synth_clifford_full = synth.synth_clifford_full


"""
    run_qiskit_baseline(code_params::CodeParameters, network_specs::NetworkSpecifications)::Tuple{Vector{AbstractOperation}, Vector{Int}, Bool}

Synthesises an encoding circuit for the provided CSS code using Qiskit's `synth_clifford_full` decomposition.

### Input

- `code_params` -- parameters defining the target QEC code
- `network_specs` -- hardware networking specifications for the underlying DQC architecture

### Output

Returns a 3-element tuple containing the retrieved encoding circuit, its corresponding gate counts and verification.

### Notes

We invoke the circuit synthesis via PyCall, which we provide with the stabiliser of our `qec_code` in the 
required string format (see `stab_to_str`). The synthesis itself is called with the `greedy` method, and 
the resulting QASM string is converted back to a `QuantumClifford.AbstractOperation` object.
"""
function run_qiskit_baseline(code_params::CodeParameters, network_specs::NetworkSpecifications)::Tuple{Vector{QuantumClifford.AbstractOperation}, Vector{Int}, Bool}
    destabiliser = destabilizerview(MixedDestabilizer(code_params.qec_code))
    logical_x = logicalxview(MixedDestabilizer(code_params.qec_code))
    stabiliser = stabilizerview(MixedDestabilizer(code_params.qec_code))
    logical_z = logicalzview(MixedDestabilizer(code_params.qec_code))
    stab_str = stab_to_str(destabiliser, stabiliser, logical_x, logical_z)
    cliff = Clifford(stab_str) 
    qiskit_encoding_circ = synth_clifford_full(cliff, "greedy") 
    @info qiskit_encoding_circ.draw(output="text")
    quantum_clifford_encoding_circ = Vector{QuantumClifford.AbstractOperation}()
    num_q = 0
    for reg in qiskit_encoding_circ.qregs
        if reg.name == "q"
            num_q = reg.size
        else
            @error("Unknown register type encountered in Qiskit circuit")
        end
    end
    @assert num_q == network_specs.num_data_qubits
    for instruction in qiskit_encoding_circ.data
        gate = instruction.operation.name
        if gate == "h"
            qubits = instruction.qubits
            bit_info = qiskit_encoding_circ.find_bit(qubits[1])
            index = Int(bit_info.index)
            push!(quantum_clifford_encoding_circ, sHadamard(index+1))
        elseif gate == "x"
            qubits = instruction.qubits
            bit_info = qiskit_encoding_circ.find_bit(qubits[1])
            index = Int(bit_info.index)
            push!(quantum_clifford_encoding_circ, sX(index+1))
        elseif gate == "y"
            qubits = instruction.qubits
            bit_info = qiskit_encoding_circ.find_bit(qubits[1])
            index = Int(bit_info.index)
            push!(quantum_clifford_encoding_circ, sY(index+1))
        elseif gate == "z"
            qubits = instruction.qubits
            bit_info = qiskit_encoding_circ.find_bit(qubits[1])
            index = Int(bit_info.index)
            push!(quantum_clifford_encoding_circ, sZ(index+1))
        elseif gate == "s"
            qubits = instruction.qubits
            bit_info = qiskit_encoding_circ.find_bit(qubits[1])
            index = Int(bit_info.index)
            push!(quantum_clifford_encoding_circ, sPhase(index+1))
        elseif gate == "cx"
            qubits = instruction.qubits
            control_info = qiskit_encoding_circ.find_bit(qubits[1])
            control = Int(control_info.index)
            target_info = qiskit_encoding_circ.find_bit(qubits[2])
            target = Int(target_info.index)
            push!(quantum_clifford_encoding_circ, sCNOT(control+1,target+1))
        elseif gate == "swap"
            qubits = instruction.qubits
            control_info = qiskit_encoding_circ.find_bit(qubits[1])
            control = Int(control_info.index)
            target_info = qiskit_encoding_circ.find_bit(qubits[2])
            target = Int(target_info.index)
            push!(quantum_clifford_encoding_circ, sSWAP(control+1,target+1))
        else
            @error("Unknown gate type encountered in Qiskit circuit")
        end
    end
    verification_logical_state = Helper.verify_success(quantum_clifford_encoding_circ, code_params.target_state, network_specs)
    @info "Verification of Qiskit Greedy Clifford compilation encoding circuit successful: $verification_logical_state"
    gcounts = Helper.gate_counts(quantum_clifford_encoding_circ, network_specs)
    @info "DQC gate counts for Qiskit encoding is $gcounts"
    return quantum_clifford_encoding_circ, gcounts, verification_logical_state
end 


"""
    stab_to_str(destab::Stabilizer, stab::Stabilizer, logical_x::Stabilizer, logical_z::Stabilizer)::Vector{String}

Convert the components of a `MixedDestabilizer` tableau into a vector of Pauli string representations that
can be used as input to Qiskit's cirucit synthesis.

### Input

- `destab` -- the destabilizer generators of the target CSS code
- `stab` -- the stabilizer generators of the target CSS code
- `logical_x` -- the logical X operators
- `logical_z` -- the logical Z operators

### Output

Returns a vector of strings, where each string translates a stabilizer row into its complete
Pauli components (e.g., `["+XIXIXIX", "+IXXIIXX", ...]`).

### Notes

Because Qiskit uses reverse qubit ordering for its tableau representation, the Pauli sequence is transcribed right-to-left (`ncols:-1:1`).
The phase bits are translated to a prepended `+` or `-` symbol depending on the 0x00 / 0x02 `tab(x).phases` values.
Since we focus on CSS codes, stabilizers will only contain pure `X` strings, pure `Z` strings, and identities (`I`).

### Examples

The stabiliser 
    
        + X_X_X_X
        + _XX__XX
        + ___XXXX
        + Z_ZZ_Z_
        + ZZ__ZZ_
        + ZZ_Z__Z
        + ZZZZZZZ

is translated to the vector of strings
    
    ["XIXIXIX", "IXXIIXX", "IIIXXXX", "ZIZZIZI", "ZZIIZZI", "ZZIZIIZ", "ZZZZZZ"].

"""
function stab_to_str(destab::Stabilizer, stab::Stabilizer, logical_x::Stabilizer, logical_z::Stabilizer)::Vector{String}
    ncols = nqubits(stab)
    clifford = String[]
    function row_string(row, phase)
        if phase == 0x00
            row_str = "+"
        elseif phase == 0x02
            row_str = "-"
        end
        for j in ncols:-1:1 # qiskit uses reverse qubit ordering for its tableau representation
            pauli = row[j]
            if pauli == (true, false) 
                row_str *= "X"
            elseif pauli == (false, true)
                row_str *= "Z"
            elseif pauli == (true, true)
                row_str *= "Y"
            else
                row_str *= "I"
            end
        end
        return row_str
    end
    for i in 1:length(destab)
        row = row_string(destab[i], tab(destab).phases[i])
        push!(clifford, row)
    end
     for i in 1:length(logical_x)
        row = row_string(logical_x[i], tab(logical_x).phases[i])
        push!(clifford, row)
    end
    for i in 1:length(stab)
        row = row_string(stab[i], tab(stab).phases[i])
        push!(clifford, row)
    end
    for i in 1:length(logical_z)
        row = row_string(logical_z[i], tab(logical_z).phases[i])
        push!(clifford, row)
    end
    return clifford
end


"""
    run_mqt_baseline(code_params::CodeParameters, network_specs::NetworkSpecifications, 
                        mqt_path::String, prep_method::String)::Tuple{Vector{AbstractOperation}, Vector{Int}, Bool}

Synthesises an encoding circuit as generated by the Munich Quantum Toolkit QECC repository.

### Input

- `code_params` -- parameters defining the target QEC code
- `network_specs` -- hardware networking specifications 
- `mqt_path` -- absolute system path mapping to the local Python virtual environment containing the MQT QECC setup
- `prep_method` -- algorithm specification passed directly to MQT ("optimal" or "heuristic")

### Output

Returns a 3-element tuple containing the retrieved encoding circuit, its corresponding gate counts and verification.

### Notes

The MQT encoding is orchestrated via a Python subprocess calling `scripts/state_encoding.py`, mainly interacting with 
the circuit synthesis tools defined in https://github.com/munich-quantum-toolkit/qecc/tree/main/src/mqt/qecc/circuit_synthesis.
The resulting Python QASM string is converted to a `QuantumClifford.AbstractOperation` object and evaluated.
We restrict to all-positive stabiliser generator phases, and by the design of the MQT algorithm, we expect only 
Hadamard and CNOT gates to be returned.
""" 
function run_mqt_baseline(code_params::CodeParameters, network_specs::NetworkSpecifications,
                            mqt_path::String, prep_method::String)::Tuple{Vector{QuantumClifford.AbstractOperation}, Vector{Int}, Bool}
    stabiliser = stabilizerview(MixedDestabilizer(code_params.qec_code))
    if any(phases(stabiliser) .== 0x02) 
        error("Baseline encoding currently only implemented for all-positive phases")
    end
    hx = stab_to_hx(stabiliser)
    hx_json = JSON.json(hx) # for passing the check matrix to the Python library
    python_bin = joinpath(mqt_path, ".venv/bin/python3")
    script_path =  joinpath(mqt_path, "scripts/state_encoding.py")
    mqt_encoding_circ_qasm = readchomp(`$(python_bin) $(script_path) $hx_json $(code_params.distance) $(prep_method)`)    
    mqt_encoding_circ = qasm2.loads(mqt_encoding_circ_qasm)     
    quantum_clifford_encoding_circ = Vector{QuantumClifford.AbstractOperation}()
    num_q = 0
    for reg in mqt_encoding_circ.qregs
        if reg.name == "q"
            num_q = reg.size
        else
            @error("Unknown register type encountered in Qiskit circuit")
        end
    end
    @assert num_q == network_specs.num_data_qubits
    for instruction in mqt_encoding_circ.data
        gate = instruction.operation.name        
        if gate == "h"
            qubits = instruction.qubits
            bit_info = mqt_encoding_circ.find_bit(qubits[1])
            index = Int(bit_info.index)
            push!(quantum_clifford_encoding_circ, sHadamard(index+1))
        elseif gate == "cx"
            qubits = instruction.qubits
            control_info = mqt_encoding_circ.find_bit(qubits[1])
            control = Int(control_info.index)
            target_info = mqt_encoding_circ.find_bit(qubits[2])
            target = Int(target_info.index)
            push!(quantum_clifford_encoding_circ, sCNOT(control+1,target+1))
        else
            error("Unsupported gate type in MQT Encoding circuit: $gate")
        end
    end
    verification_logical_state = Helper.verify_success(quantum_clifford_encoding_circ, code_params.target_state, network_specs)
    @info "Verification of MQT encoding circuit successful: $verification_logical_state"
    gcounts = Helper.gate_counts(quantum_clifford_encoding_circ, network_specs)
    return quantum_clifford_encoding_circ, gcounts, verification_logical_state
end


"""
    stab_to_hx(stab::Stabilizer)::Matrix{Int8}

Convert the components of a `stabilizerview' object into a a bit matrix 
representing the X-type part of the stabiliser generators.

### Input

- `stab` -- the stabilizer generators of the target CSS code

### Output

Returns the bit matrix corresponding to the X-type part of the tableau `stab`

### Notes

By the virtue of encoding the logical zero state of a CSS code, it suffices to
implement the X-type stabilisers of the tableau. Because the MQT library uses bit 
matrices to specify tableaus, we convert a given stabiliser to the corresponding matrix. 
Phases are not treated here, so we require that the input `stab` has all-positive phases.

### Examples

The tableau 
    
        + X_X_X_X
        + _XX__XX
        + ___XXXX
        + Z_ZZ_Z_
        + ZZ__ZZ_
        + ZZ_Z__Z

returns the X-type check matrix 
    
    Int8[1 0 1 0 1 0 1; 
         0 1 1 0 0 1 1; 
         0 0 0 1 1 1 1]
"""
function stab_to_hx(stab::Stabilizer)::Matrix{Int8}
    nrows = length(stab)
    ncols = nqubits(stab)
    hx = zeros(Int8, 0, ncols)
    for i in 1:nrows
        row = zeros(Int8, ncols)
        x_row = false
        for j in 1:ncols
            pauli = stab[i][j]
            if pauli == (true, false) 
                row[j] = 1
                x_row = true
            else
                @assert pauli == (false, false) || (pauli == (false, true) && !x_row)
            end
        end
        if x_row
            hx = vcat(hx, reshape(row, 1, ncols))
        end
    end 
    return np.array(hx, dtype=np.int8)
end



end
