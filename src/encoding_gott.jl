""" Credit: Almost entirely based on https://github.com/QuantumSavory/QuantumClifford.jl/blob/4d524965a11b6d8594d578d4935ccd7bc385f56c/src/ecc/circuits.jl"""

module EncodingGott

using ..Types
using ..Helper

using QuantumClifford
using QECCore
using Serialization
using Quantikz: savecircuit
using QuantumClifford: true_success_stat, false_success_stat, continue_stat, failure_stat

export encoding_circuit_gott, encoding_gott

"""
This function creates an encoding circuit based on the procedure described by Gottesman. We first obtain the canonical form of the
    stabiliser tableau, without re-permuting the qubits ( MixedDestabilizer(code, undoperm=false) ), while reporting the permutation
    operation (see `canonicalize_gott!`). Next, the logical operators are extracted from this canonical form, and the encoding circuit
    is built. Lastly, the qubits (in the encoded circuit) are permuted back to match the original qubit ordering.

    The resulting circuit encodes any state on physical qubits at indices 'n-k+1:n' into a larger logical code by creation of (non fault-tolerant) encoding circuit.
    If we simply start from the all-zero physical register, it encodes the logical all-zero state.

By calling MixedDestabilizer(code, undoperm=true), we can additionally directly extract the standard form of the tableau with its logical operators,
    with the correct qubit indexing (i.e., after repermuting).  The stabiliser generators here might be different as a result of Guassian elimination, 
    which however is not problematic since this new set of generators still stabilises the codespace (columns swaps will be restored correctly).
    (Also, the logical operators match both forms of generators, since those generators generate the same stabiliser group.)
    Using the backtrack argument to MixedDestabilizer, one is also able to restore the entire original tableau, which however is not necessary for our purposes,
    since we work with the final canonical tableau, after column swaps for state preparation).
    
The implementation is based on [cleve1997efficient] and [gottesman1997stabilizer], which however is partly erroneous
    (see https://perimeterinstitute.ca/personal/dgottesman/thesis-errata.html), [grassl2002algorithmic] and [grassl2011variations] )
"""


function _gottesman_encoding_circuit_raw(code, basis_state)#; reportperm = true)#, undoperm=false)
         
    # Creatung the Gottesman encoding circuit (without re-permuting the columns)
    n = code_n(code)
    k = code_k(code)
    @assert length(basis_state) == k
    
    code_standard_form, r, permx, permz = MixedDestabilizer(code, undoperm=false, reportperm=true); # undoperm without returns gives the orignal stabiliser tableau
  
    T = stabilizerview(code_standard_form)
    X = logicalxview(code_standard_form)  
    
    circ = QuantumClifford.AbstractOperation[]

# ----- Logical Xs -----
    
    for i in 1:k
        if basis_state[i] == 1
            push!(circ, sX(n-k+i))
            for j in 1:n-k-r   # 1:s
                if X[i,r+j][1] == true
                    push!(circ, sCNOT(n-k+i, r+j))
                end
            end
        end
        
    end

# ----- Projection onto codespace -----

    for i in 1:r
        push!(circ, sHadamard(i))
        if T[i,i][2] == true
            push!(circ, sPhase(i))
        end
        for j in 1:n
            if i!=j
                xz = T[i,j]
                g = if xz == (true, true)  # controlled-Y
                    sZCY
                elseif xz == (true, false) # CNOT (controlled-X)
                    sZCX
                elseif xz == (false, true) && !(i<j<n-k+1) # controlled-Z
                    sZCZ
                end
                isnothing(g) || push!(circ, g(i,j))
            end
        end
    end

# ----- Adapting for negative phases in the tableau -----
    for i in 1:n-k  # 1:r+s
        if phases(T)[i]!=0
            if i<=r
                push!(circ, sZ(i))
            else
                push!(circ, sX(i))
            end
        end
    end
    
    gottesman_perm = permz[permx]

"""
For example:
    permx = [1,3,2,4] says: "for x canonicalization, qubit 3 was placed in position 2"
    permz = [1,2,4,3] says: "for z canonicalization, qubit 4 was placed in position 3"
    We performed permx first, then permz.
    Overall, this means [1,2,3,4]>(permx)>[1,3,2,4]>(permz)>[1,3,4,2]
    Now this is precisely permz[permx] = permz[1,3,2,4] = [1,3,4,2]
"""
       
    # if undoperm
    #     transpositions = perm_to_transpositions(gottesman_perm)

    #     for (i,j) in reverse(transpositions) # transpositions capture the permutatation by right application, i.e., from right to left
    #         push!(circ, sSWAP(i,j))
    #     end

    
    return circ, gottesman_perm
end

function encoding_circuit_gott(qec_code, network_specs, basis_state)

    @assert qec_code !== nothing 

    encoding_circuit_gott, perm = _gottesman_encoding_circuit_raw(qec_code, basis_state)
"""
We wish to recover the original ordering from the canonicalization permutation.
The permutation expresses, which index i gets qubit j. Now by the virtue of canonicalization, we need to make
sure that an opertaion performed on index i is actually applied to qubit j. Thus, we relabel the affected qubit of 
an operation from op.q to perm[op.q] (instead of introducing SWAP operations)
"""
    encoding_circuit_orig = QuantumClifford.AbstractOperation[]
    gate_counts = [0,0,0]
    for op in encoding_circuit_gott
        
        T = typeof(op)
        if T <: AbstractSingleQubitOperator
            push!(encoding_circuit_orig, T(perm[op.q]))
            gate_counts += [1,0,0]
        elseif T<: AbstractTwoQubitOperator
            p_control = perm[op.q1]
            p_target = perm[op.q2]
            # perm is the permutation that identifies the qubit to apply the gate to (from permuted encoding circuit),
            # but to identify telegates, we need the qubit mapping (here: the inverse mapping) from the DQC network specifications
            push!(encoding_circuit_orig, T(p_control, p_target))
            if network_specs.register_lookup_array[network_specs.inv_map[p_control]] == network_specs.register_lookup_array[network_specs.inv_map[p_target]]
                # the above is equivalent to network_specs.register_lookup_array[findfirst(==(p_target), network_specs.permutation)]
                gate_counts += [0,1,0]
            else
                gate_counts += [0,0,1]
            end
        end
    end
       
    circuit_sizes = gate_counts
    
    @assert length(encoding_circuit_orig) == sum(circuit_sizes) 
  
    return encoding_circuit_orig, gate_counts
end

function encoding_gott(code_params, network_specs, basis_state, folder)

    encoding_circ, gate_counts = encoding_circuit_gott(code_params.qec_code, network_specs, basis_state)

    @info "Gottesman encoding circuit length in DQC setting:: Single-qubit gates: $(gate_counts[1]), Two-qubit gates: $(gate_counts[2]), Telegates: $(gate_counts[3])"

    # ----- Verification ------
    verification_logical_state = verify_success(encoding_circ, code_params.target_state, network_specs)
    @info "Verification of Gottesman circuit successful: $verification_logical_state"

    # ----- Data Storage ----------
    dir = joinpath(folder, "gottesman_encoding")
    mkpath(dir)

    serialize( joinpath(dir, "encoding_circuit.jls"), encoding_circ )

    open(joinpath(dir, "encoding_gates.txt"), "w") do io
        println(io, "# Raw gate sequence of size $(sum(gate_counts))")
        for (i, g) in enumerate(encoding_circ)
            println(io, i, "\t", repr(g))
        end
    end

    save_circuit_diagram(encoding_circ, dir, "encoding_circuit.png")

    open(joinpath(dir, "summary.txt"), "w") do io
        println(io, "# Encoding successful: $verification_logical_state")
        println(io, "# Raw gate sequence of size $(sum(gate_counts))")
        println(io, "# Executable DQC circuit with $(gate_counts[1]) single qubit gates, $(gate_counts[2]) two qubit gates and $(gate_counts[3]) telegates ")
    end

    return encoding_circ
end

end