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

function run_mqt_baseline(code_params, network_specs, folder)
   
    stabiliser = stabilizerview(MixedDestabilizer(code_params.qec_code))
    if any(phases(stabiliser) .== 0x02) 
        error("Baseline encoding currently only implemented for all-positive phases")
    end
    hx = stab_to_hx(stabiliser)
    #println("Converted $stabiliser to $hx")


    hx_json = JSON.json(hx) # for passing the check matrix to Python library

    print(@__DIR__)
    print(PyCall.python)
#    mqt_encoding_circ_qasm = readchomp(`/Users/tim/Tim/projects/mqt/qecc/.venv/bin/python3 /Users/tim/Tim/projects/mqt/qecc/scripts/state_encoding.py $hx_json $(code_params.distance)`)

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

    # The MQT library return Hadamard and CNOT gates only

    for instruction in mqt_encoding_circ.data
        #println(verification_circ.data)
        gate = instruction.operation.name
        
        if gate == "h"
            qubits = instruction.qubits
            bit_info = mqt_encoding_circ.find_bit(qubits[1])
            # reg_name = String(bit_info.registers[1][1].name)
            # if reg_name == "q"
            #     continue
            # end
            index = Int(bit_info.index)
            
            #println(gate, index, reg_name)
            # qc_index = 0 
            # if reg_name == "z_anc"
            #     qc_index = num_q + index
            # elseif reg_name == "x_anc"
            #     qc_index = num_q + num_z_anc + index
            # elseif reg_name == "flag"
            #     qc_index = num_q + num_z_anc + num_x_anc + index
            # end

            push!(quantum_clifford_encoding_circ, sHadamard(index+1))

        elseif gate == "cx"
            qubits = instruction.qubits
            control_info = mqt_encoding_circ.find_bit(qubits[1])
            #control_reg_name = String(control_info.registers[1][1].name)
            control = Int(control_info.index)
            target_info = mqt_encoding_circ.find_bit(qubits[2])
            #target_reg_name = String(target_info.registers[1][1].name)
            target = Int(target_info.index)
            push!(quantum_clifford_encoding_circ, sCNOT(control+1,target+1))
        else
            error("Unsupported gate type in MQT Encoding circuit: $gate")
        end
    end

    verification_logical_state = Helper.verify_success(quantum_clifford_encoding_circ, code_params.target_state, network_specs)
    @info "Verification of MQT encoding circuit successful: $verification_logical_state"

    gcounts = Helper.gate_counts(quantum_clifford_encoding_circ, network_specs)

    #println("DQC gate counts for MQT encoding is $gcounts")



    # ----- Data Storage ----------
    dir = joinpath(folder, "mqt_encoding")
    mkpath(dir)

    serialize( joinpath(dir, "mqt_encoding_circuit.jls"), quantum_clifford_encoding_circ )
    #serialize( joinpath(dir, "encoding_circuit_dqc_compiled.jls"), encoding_circ_compiled )

    # open(joinpath(dir, "encoding_gates.txt"), "w") do io
    #     println(io, "# Raw gate sequence (Gottesman encoding circuit) of size $(sum(gate_counts))")
    #     for (i, g) in enumerate(encoding_circ)
    #         println(io, i, "\t", repr(g))
    #     end

    #     println(io, "# Raw gate sequence (DQC compiled version) of size $(sum(gate_counts_compiled))")
    #     for (i, g) in enumerate(encoding_circ_compiled)
    #         println(io, i, "\t", repr(g))
    #     end
    # end

    save_circuit_diagram(quantum_clifford_encoding_circ, dir, "mqt_encoding_circuit.png")
    #save_circuit_diagram(encoding_circ_compiled, dir, "encoding_circuit_dqc_compiled.png")


    # open(joinpath(dir, "summary.txt"), "w") do io
    #     println(io, "# Encoding successful: $verification_logical_state")
    #     println(io, "# Raw gate sequence of size $(sum(gcounts))")
    #     println(io, "# Executable DQC circuit with $(gcounts[1]) single qubit gates, $(gcounts[2]) two qubit gates and $(gcounts[3]) telegates ")
    #     # println(io, "\n")
    #     # println(io, "# DQC Compilation Encoding successful: $verification_logical_state_compiled")
    #     # println(io, "# Raw gate sequence of size $(sum(gate_counts_compiled))")
    #     # println(io, "# Executable DQC circuit with $(gate_counts_compiled[1]) single qubit gates, $(gate_counts_compiled[2]) two qubit gates and $(gate_counts_compiled[3]) telegates ")
    # end

    df = DataFrame(method = ["mqt_encoding"], verified = [verification_logical_state], gate_counts = [gcounts])
    CSV.write(joinpath(dir, "mqt_encoding_stats.csv"), df)

    return mqt_encoding_circ, gcounts

end


function run_qiskit_baseline(code_params, network_specs, folder)
   
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
            @error("Unknown register type encountered in Qiskit circuit")
        end
    end

    verification_logical_state = Helper.verify_success(quantum_clifford_encoding_circ, code_params.target_state, network_specs)
    @info "Verification of Qiskit Greedy Clifford compilation encoding circuit successful: $verification_logical_state"

    gcounts = Helper.gate_counts(quantum_clifford_encoding_circ, network_specs)
    @info "DQC gate counts for Qiskit encoding is $gcounts"


    # ----- Data Storage ----------
    dir = joinpath(folder, "qiskit_encoding")
    mkpath(dir)
    serialize( joinpath(dir, "qiskit_encoding_circuit.jls"), quantum_clifford_encoding_circ )
    save_circuit_diagram(quantum_clifford_encoding_circ, dir, "qiskit_encoding_circuit.png")
    df = DataFrame(method = ["qiskit_encoding"], verified = [verification_logical_state], gate_counts = [gcounts])
    CSV.write(joinpath(dir, "qiskit_encoding_stats.csv"), df)

    return qiskit_encoding_circ, gcounts

end 





end
#using PyCall

# np = pyimport("numpy")
# logging = pyimport("logging")
# logging.basicConfig(
#     level=logging.INFO,
#     format="%(asctime)s %(name)s %(levelname)s: %(message)s"
# )
# logging.getLogger("mqt.qecc").setLevel(logging.INFO)

# CSSCode = pyimport("mqt.qecc").CSSCode
# cs = pyimport("mqt.qecc.circuit_synthesis")
# # sim = pyimport("mqt.qecc.simulation")
# # noise = pyimport("mqt.qecc.noise")

# qiskit = pyimport("qiskit")
# qi = pyimport("qiskit.quantum_info")
# qasm2 = pyimport("qiskit.qasm2")
# synth = pyimport("qiskit.synthesis")

# heuristic_prep_circuit         = cs.heuristic_prep_circuit
# gate_optimal_verification_circuit = cs.gate_optimal_verification_circuit
# VerificationNDFTStatePrepSimulator = cs.VerificationNDFTStatePrepSimulator
# CircuitLevelNoiseIdlingParallel    = cs.CircuitLevelNoiseIdlingParallel


# Define matrices — note PyCall needs explicit int8 type
# hx = np.array([
#     [0,0,1,1,0,0,1,1,0,0,0,0],
#     [1,0,0,0,1,0,0,1,1,0,0,0],
#     [0,1,0,0,0,1,1,0,1,0,0,0],
#     [1,0,0,0,0,1,0,0,0,1,1,0],
#     [0,1,0,1,0,0,0,0,0,0,1,1],
#     [0,0,1,0,1,0,0,0,0,1,0,1],
# ], dtype=np.int8)

# hz = np.array([
#     [1,0,1,0,0,0,0,1,0,1,0,0],
#     [1,1,0,0,0,0,0,0,1,0,1,0],
#     [0,1,1,0,0,0,1,0,0,0,0,1],
#     [0,0,0,1,0,1,1,0,0,0,1,0],
#     [0,0,0,1,1,0,0,1,0,0,0,1],
#     [0,0,0,0,1,1,0,0,1,1,0,0],
# ], dtype=np.int8)

# code = CSSCode(hx,hz,3)

# # code.Lx = np.array([
#     [0,0,0,0,0,1,0,0,1,0,1,0],
#     [0,0,0,0,0,1,1,1,0,1,0,1],
# ])
# code.Lz = np.array([
#     [1,0,1,1,1,0,0,0,0,0,1,0],
#     [0,0,1,1,0,0,0,0,0,0,0,1],
# ])



# println(code.n)

# non_ft_sp = heuristic_prep_circuit(code, zero_state=true)
# println(non_ft_sp)

# plt = pyimport("matplotlib.pyplot")
# #non_ft_sp.draw(output="mpl", initial_state=true, fold=-1, scale=0.4)



# ft_sp = gate_optimal_verification_circuit(non_ft_sp, verify_x_first = true)
# print(ft_sp)
# fig = ft_sp[:draw](
#     output="mpl",
#     initial_state=true,
#     fold=-1,
#     scale=0.4
# )
# plt.show()
# qasm_prog = qasm2.dumps(ft_sp)
# print(qasm_prog)

# p = 1e-3
# noise = CircuitLevelNoiseIdlingParallel(
#     p_tqg=p, p_sqg=p, p_init=p, p_meas=p, p_idle=p/100
# )

# non_ft_sim = VerificationNDFTStatePrepSimulator(non_ft_sp.circ, code=code, zero_state=true)
# ft_sim     = VerificationNDFTStatePrepSimulator(ft_sp,           code=code, zero_state=true)

# pl_non_ft, ra_non_ft = non_ft_sim.logical_error_rate(noise, min_errors=10)[1:2]
# pl_ft,     ra_ft     = ft_sim.logical_error_rate(noise,     min_errors=10)[1:2]

# println("Non-FT logical error rate: $pl_non_ft")
# println("FT logical error rate:     $pl_ft")