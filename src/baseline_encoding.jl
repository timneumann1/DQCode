module BaselineEncoding

include("types.jl")
include("trivariate_bicycle_code.jl")
include("helper.jl")

using .Types
using .Helper

using QuantumClifford
using Serialization
using PyCall
using JSON
using CSV, DataFrames

# Python packages, to be retrieved via PyCall
np = pyimport("numpy")
qiskit = pyimport("qiskit")
qi = pyimport("qiskit.quantum_info")
qasm2 = pyimport("qiskit.qasm2")
synth = pyimport("qiskit.synthesis")
plt = pyimport("matplotlib.pyplot")
Clifford = qi.Clifford
synth_clifford_full = synth.synth_clifford_full

export run_qiskit_baseline, run_mqt_baseline

function stab_to_str(destab::Stabilizer, stab::Stabilizer, logical_x::Stabilizer, logical_z::Stabilizer)

    """
      Since we are working with CSS codes, we can assume that any stabiliser that contains an X ONLY contains X, and no other stabiliser contins X
    accepts stabilizer like

        + X_X_X_X
        + _XX__XX
        + ___XXXX
        + Z_ZZ_Z_
        + ZZ__ZZ_
        + ZZ_Z__Z
        + ZZZZZZZ

    and (depending on string Bool) returns hx check matrix in string form like
    
    ["XIXIXIX", "IXXIIXX", "IIIXXXX", "ZIZZIZI", "ZZIIZZI", "ZZIZIIZ", "ZZZZZZ],
    """

    ncols = nqubits(stab)

    # Return string form like ["+XI", "+XX"]
    clifford = String[]

    function row_string(row, phase)
        if phase == 0x00
            row_str = "+"
        elseif phase == 0x02
            row_str = "-"
        end

        for j in ncols:-1:1 # qiskit ordering is reverse
            pauli = row[j]
            if pauli == (true, false)  # Check if X component is set at qubit j
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
        #print(destab[i][1])
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

function stab_to_hx(stab::Stabilizer)
    """
    Since we are working with CSS codes, we can assume that any stabiliser that contains an X ONLY contains X, and no other stabiliser contins X
    accepts stabilizer like

        + X_X_X_X
        + _XX__XX
        + ___XXXX
        + Z_ZZ_Z_
        + ZZ__ZZ_
        + ZZ_Z__Z

    and returns hx check matrix 
    
    in bit form (dtype np.int8) like

        hx = np.array([
        [0,0,1,1,0,0,1,1,0,0,0,0],
        [1,0,0,0,1,0,0,1,1,0,0,0],
        [0,1,0,0,0,1,1,0,1,0,0,0],
        [1,0,0,0,0,1,0,0,0,1,1,0],
        [0,1,0,1,0,0,0,0,0,0,1,1],
        [0,0,1,0,1,0,0,0,0,1,0,1],
    ], dtype=np.int8)

    """

    nrows = length(stab)
    ncols = nqubits(stab)

    # Return bit matrix as numpy int8 array
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
    
    # Convert to numpy array for compatibility with downstream Python code
    return np.array(hx, dtype=np.int8)
    

end

function run_mqt_baseline(code_params, network_specs, mqt_path, prep_method)
   
    stabiliser = stabilizerview(MixedDestabilizer(code_params.qec_code))
    if any(phases(stabiliser) .== 0x02) 
        error("Baseline encoding currently only implemented for all-positive phases")
    end
    hx = stab_to_hx(stabiliser)

    hx_json = JSON.json(hx) # for passing the check matrix to Python library

    python_bin = joinpath(mqt_path, ".venv/bin/python3")
    script_path =  joinpath(mqt_path, "scripts/state_encoding.py")
    mqt_encoding_circ_qasm = readchomp(`$(python_bin) $(script_path) $hx_json $(code_params.distance) $(prep_method)`)    
    mqt_encoding_circ = qasm2.loads(mqt_encoding_circ_qasm) 
    
    # Build QuanutmClifford circuit from QASM
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

    # The MQT library returns Hadamard and CNOT gates only
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


function run_qiskit_baseline(code_params, network_specs)
   
    destabiliser = destabilizerview(MixedDestabilizer(code_params.qec_code))
    logical_x = logicalxview(MixedDestabilizer(code_params.qec_code))
    stabiliser = stabilizerview(MixedDestabilizer(code_params.qec_code))
    logical_z = logicalzview(MixedDestabilizer(code_params.qec_code))
    stab_str = stab_to_str(destabiliser, stabiliser, logical_x, logical_z)
    
    cliff = Clifford(stab_str) 
    qiskit_encoding_circ = synth_clifford_full(cliff, "greedy") # "ag", "greedy"
    @info qiskit_encoding_circ.draw(output="text")
    
    # Build QuanutmClifford circuit from QASM

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

end
