""" Credit: https://github.com/QuantumSavory/QuantumClifford.jl/blob/4d524965a11b6d8594d578d4935ccd7bc385f56c/src/ecc/circuits.jl"""

module BaselineEncoding

using ..Types
using ..Helper

using QuantumClifford
using QECCore
using Serialization
using Quantikz: savecircuit
using QuantumClifford: true_success_stat, false_success_stat, continue_stat, failure_stat


export baseline_encoding, baseline_encoding_circuit

"""
This function creates the naive encoding circuit by first obtaining the canonical form of the stabiliser tableau, without re-permuting
    the qubits ( MixedDestabilizer(code, undoperm=false) ), but reporting the permutation operation (see [`canonicalize_gott!`](@ref)). Next, the logical operators are extracted from this canonical form, and the encoding circuit
    is built. Lastly, the qubits (in the encoded circuit) are permuted back to match the creation of state of the original tableau.

    The resulting circuit encodes any state on physical qubits at indices 'n-k+1:n' into a larger logical code by creation of (non fault-tolerant) encoding circuit.
    If we simply start from the all-zero physical register, it encodes the logical all-zero state.

By calling MixedDestabilizer(code, undoperm=true), we can additionally directly extract the standard form of the tableau with its logical operators,
    with the correct qubit indexing (i.e., after repermuting).  The stabiliser generators here might be different as a result of Guassian elimination, 
    which however is not problematic since this new set of generators still stabilises the codespace (columns swaps will be restored correctly).
    (Also, the logical operators match both forms of generators, since those generators generate the same stabiliser group.)
    Using the backtrack argument to MixedDestabilizer, one is also able to restore the entire original tableau, which however is not necessary for our purposes,
    since we work with the final canonical tableau, after column swaps for state preparation).
    
The implementation is based on [cleve1997efficient](@cite) and [gottesman1997stabilizer](@cite), which however is partly erroneous
    (see https://perimeterinstitute.ca/personal/dgottesman/thesis-errata.html), [grassl2002algorithmic](@cite) and [grassl2011variations](@cite) )
"""


function _standard_logical_zero_encoding_circuit(code; undoperm=true, logical_Xs = false)
     
    # Differs from the above by removing the logical X operator CNOTs

    # Creatung the canonical tableau (without re-permuting the columns)
    n = code_n(code)
    k = code_k(code)
    #println("\nFor the given code, we have n=$n, k=$k. \n")
    
    code_standard_form, r, permx, permz = MixedDestabilizer(code, undoperm=false, reportperm=true); # undoperm without returns gives the orignal stabiliser tableau
    #println("Standard form of code is \n$code_standard_form")
    X = logicalxview(code_standard_form)
    Z = logicalzview(code_standard_form)
    
    # Creating te canonical tableau (with re-permuting the columns) for final state specification
    #code_original_with_logicals = MixedDestabilizer(code, undoperm=true);
    # X_cleaned = logicalxview(md_cleaned)
    # println("X is $X_cleaned")
    # Z_cleaned = logicalzview(md_cleaned)
    # println("Z is $Z_cleaned")

    circ = QuantumClifford.AbstractOperation[]
    #push!(circ, Reset(initial_state, [1,2,3,4,5,6,7]))
    
    # Constructing the encoding circuit

    S = stabilizerview(code_standard_form)

    # logical Xs
    if logical_Xs
        for i in 1:k
            for t in 1:n-k
                if X[i,t][1] == true
                    push!(circ, sCNOT(n-k+i, t))
                end
            end
        end
    end
    # projection on codespace
    for i in 1:r
        push!(circ, sHadamard(i))
        if S[i,i][2] == true
            push!(circ, sPhase(i))
        end
        for t in 1:n
            if i!=t
                xz = S[i,t]
                g = if xz == (true, true)  # controlled-Y
                    sZCY
                elseif xz == (true, false) # CNOT (controlled-X)
                    sZCX
                elseif xz == (false, true) && !(i<t<n-k+1) # controlled-Z
                    sZCZ
                end
                isnothing(g) || push!(circ, g(i,t))
            end
        end
    end

    # correct for negative phases in the tableau
    for i in 1:n-k 
        if phases(S)[i]!=0
            if i<=r
                push!(circ, sZ(i))
            else
                push!(circ, sX(i))
            end
        end
    end
    transpositions = nothing
    # undoing the permutations to have the correct circuit for the final (re-permuted) tableau
    if undoperm
        perm = permx[permz]
        transpositions = perm_to_transpositions(perm)
        for (i,j) in transpositions
            push!(circ, sSWAP(i,j))
        end
    end
    circ, transpositions #, code_original_with_logicals, 
end

function baseline_encoding_circuit(qec_code, network_specs; logical_Xs = false)

    @assert qec_code !== nothing 

    baseline_circuit_raw, transpositions = _standard_logical_zero_encoding_circuit(qec_code, logical_Xs = logical_Xs)
    #standard_circuit_length = length(standard_circuit[2])
    permutation = transpositions_to_perm(reverse(transpositions), network_specs.num_data_qubits)

    baseline_gates = Gate[]
    for op in baseline_circuit_raw
        if op isa QuantumClifford.sHadamard
            push!(baseline_gates, HadamardGate(permutation[op.q]))
        elseif op isa QuantumClifford.sPhase
            push!(baseline_gates, SGate(permutation[op.q]))
        elseif op isa QuantumClifford.sX
            push!(baseline_gates, PauliXGate(permutation[op.q]))
        elseif op isa QuantumClifford.sZ
            push!(baseline_gates, PauliZGate(permutation[op.q]))
        elseif op isa QuantumClifford.sZCX
            control, target = Tuple(affectedqubits(op))
            push!(baseline_gates, CX_Gate(permutation[control], permutation[target]))
        elseif op isa QuantumClifford.sCNOT
            control, target = Tuple(affectedqubits(op))
            push!(baseline_gates, CX_Gate(permutation[control], permutation[target]))
        # elseif op isa QuantumClifford.sZCY
        #     control, target = Tuple(affectedqubits(op))
        #     push!(gates, CNOT_Gate(permutation[control], permutation[target]))
        # elseif op isa QuantumClifford.sZCZ
        #     control, target = Tuple(affectedqubits(op))
        #     push!(gates, CNOT_Gate(permutation[control], permutation[target]))
        elseif op isa QuantumClifford.sSWAP
            continue
        else
            error("The baseline encoding circuit uses a gate type that is currently not supported: $(typeof(op))")
        end  
        
    end

    #baseline_raw_circuit = CircuitIndividual(baseline_gates)
    #println(typeof(baseline_raw_circuit))
    #baseline_exec_circuit, num_single_qubit_gates, num_two_qubit_gates, num_telegates = construct_executable_circuit(baseline_gates, gate_set, network_specs, telegate_overhead = true)
    baseline_circuit, gate_counts = gates_to_circuit(baseline_gates, network_specs)
    circuit_sizes = gate_counts
    
    @assert length(baseline_gates) == sum(circuit_sizes) 
   
    #println("Original Code with repermuted logical operators appended: $(code_original_with_logicals)")
    #println("\nRaw size of standard encoding circuit: $(sum(circuit_sizes))\n")
    #println("DQC Size of standard encoding circuit: $(circuit_sizes)\n")
    #println("Single qubit: $num_single_qubit_gates, Two qubit: $num_two_qubit_gates, Telegates: $num_telegates")
    return baseline_gates, baseline_circuit, gate_counts
end

function baseline_encoding(code_params, network_specs, folder)

    ### Baseline comparison of standard encoding in DQC setting
    baseline_gates, baseline_circuit, gate_counts = baseline_encoding_circuit(code_params.qec_code, network_specs)

    @info "Baseline encoding circuit length in DQC setting:: Single-qubit gates: $(gate_counts[1]), Two-qubit gates: $(gate_counts[2]), Telegates: $(gate_counts[3])"

    # ----- Verification ------
    verification_logical_state = verify_success(baseline_circuit, code_params.target_state, network_specs)
    @info "Verification of baseline circuit successful: $verification_logical_state"

    # ----- Data Storage ----------
    dir = joinpath(folder, "baseline_encoding")
    mkpath(dir)

    serialize( joinpath(dir, "baseline_gates.jls"), baseline_gates )

    open(joinpath(dir, "baseline_gates.txt"), "w") do io
        println(io, "# Raw gate sequence of size $(sum(gate_counts))")
        for (i, g) in enumerate(baseline_gates)
            println(io, i, "\t", repr(g))
        end
    end

    save_circuit_diagram(copy(baseline_gates), dir, "baseline_circuit.png")

    open(joinpath(dir, "summary.txt"), "w") do io
        println(io, "# Encoding successful: $verification_logical_state")
        println(io, "# Raw gate sequence of size $(sum(gate_counts))")
        println(io, "# Executable DQC circuit with $(gate_counts[1]) single qubit gates, $(gate_counts[2]) two qubit gates and $(gate_counts[3]) telegates ")
    end

    return baseline_gates
end

end


        #save_circuit_diagram(baseline_exec_circuit, dir, "baseline_exec_circuit__size_$(circuit_sizes).png")

        # open(joinpath(dir, "network_specs.txt"), "w") do io
        #     println(io, "Network Specifications")
        #     for fname in fieldnames(Types.NetworkSpecifications)
        #         println(io, fname, " = ", repr(getfield(network_specs, fname)))
        #     end
        # end

        # open(joinpath(dir, "code_params.txt"), "w") do io
        #     println(io, "Code parameters")
        #     for fname in fieldnames(Types.CodeParameters)
        #         println(io, fname, " = ", repr(getfield(code_params, fname)))
        #     end
        # end

        

        # open(joinpath(dir, "baseline_exec_circuit.txt"), "w") do io
        #     println(io, "# Executable (DQC) circuit operations of sizes $(circuit_sizes) (excl. SWAPS)")
        #     for (i, op) in enumerate(baseline_exec_circuit)
        #         println(io, i, "\t", repr(op))
        #     end
        # end

        

        ###########################################
        ###########################################
        ###########################################


