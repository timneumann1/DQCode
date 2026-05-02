""" Credit: Almost entirely based on https://github.com/QuantumSavory/QuantumClifford.jl/blob/4d524965a11b6d8594d578d4935ccd7bc385f56c/src/ecc/circuits.jl"""

module EncodingGott

using ..Types
using ..Helper

using QuantumClifford
using QECCore
using Serialization
using Quantikz: savecircuit
using QuantumClifford: true_success_stat, false_success_stat, continue_stat, failure_stat, AbstractOperation, AbstractSingleQubitOperator, AbstractTwoQubitOperator

export encoding_circuit_gott, encoding_gott, overlap_compilation

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
    
    circ = Vector{AbstractOperation}()

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

    encoding_circuit_orig = Vector{AbstractOperation}()
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

    encoding_circ_compiled, gate_counts_compiled = overlap_compilation(encoding_circ, network_specs)

    @info "Encoding Circuit: $encoding_circ"
    @info "Gottesman encoding circuit length in DQC setting:: Single-qubit gates: $(gate_counts[1]), Two-qubit gates: $(gate_counts[2]), Telegates: $(gate_counts[3])"

    @info "Encoding circuit compiled: $encoding_circ_compiled"
    @info "Overlap DQC compilation circuit length in DQC setting:: Single-qubit gates: $(gate_counts_compiled[1]), Two-qubit gates: $(gate_counts_compiled[2]), Telegates: $(gate_counts_compiled[3])"


    # ----- Verification ------
    verification_logical_state = verify_success(encoding_circ, code_params.target_state, network_specs)
    @info "Verification of Gottesman circuit successful: $verification_logical_state"

    verification_logical_state_compiled = verify_success(encoding_circ_compiled, code_params.target_state, network_specs)
    @info "Verification of Compiled circuit successful: $verification_logical_state_compiled"

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


function overlap_compilation(circuit, network_specs)

"""
The overlap is a circuit compilation method used by Paetznick and Reichardt in "Fault-tolerant ancilla preparation and noise threshold lower bounds  for the 23-qubit Golay code" [2013]
    to reduce the number of CNOTs in the encoding circuit of the Steane code. It makes use of overlapping gates to reduce certain sequences of four CNOTs to three, thus effectively reducing
    gate count by 1.

The below code is valid for CSS codes with all-positive phases, which is a very large class of QEC codes. These requirements assure that the circuits consist entirely of H + CX , where 
    a CX will never act on prior qubit based on the Gottesman canonical form, so everything commutes.

Assume that we have a pair of CX(i,k), CX(j,k) that share target qubit k. Further assume that we have a different set of CX(i,m), CX(j,m), which 
    share the same target qubit m, and whose control qubits i and j match those of the first pair of CX gates. For the central insight, if there are no other gates in between the controls that map |0> <-> |1> 
    (such as X and Y, or the target of a CX), i.e., if in between we only find controls for other gates or Z-type gates, for example, then the latter set of controls is active iff the first is. In this case, the target qubit m
    will be flipped iff the target qubit k is flipped, hence we may replace the second set of CX-gates with a CX-gate from k to m.

We use this insight to greedily find such quartets in the Gottesman encoding circuit, reducing the effective gate count. In the DQC setting, however, we don't want to reduce the gate count unconditionally, but only if 
    reducing the gate count does not actually increase the telegate count. The only case in which this happens is when k and m do not share the same core, whereas i, j and m all share the same core (in this case, we
    effectively remove two intra-core two-qubit gates and add one telegate, which is usually not desirable). Since the CX gates commute, however, we should in this case still insert a telegate between k and m, yet with 
    control m and target k, while deleting the telegates CX_ik and CX_jk. 

In order to assure that inserting new gates does not change the overall unitary, it is important to insert the new gate right at the next available spot in the compiled_circ vector. I.e., having identified gates 1-4, we insert
two of them (depending on the DQC criterion), plus the added CX gate right after. Since all CNOT gates in the Gottesman circuit we created commute, this is valid, and it ensure that there is no mixing with other gates.
However, there is a subtlety: While this procedure does the trick for one quartet, if we have another quartet that overlaps with a previous one, the inserted gate will be contaminated from earlier CNOT action that was not present 
in the original circuit. Thus, we must make sure that only those gates can be inserted that have not been before target of a CX (besides from the ones in the quartet itself, of course). This is equivalent to blocking all those insertions
that overlap in at least one qubit with a previous insertion. Since insertions are determined by the values of k and m, we add a blocking mechanism that prevents previously used qubits k and m to be accessed again for insertion. 
 
Also: Network of cleanly separated controls and targets, which we change in the course of compilation, so when grouping gates and adding them to the beginning, we need to be careful: 
    we need to perform all possible such grouping first (also respecting the above point of not reusing inserted gate qubits), and only then append the qubits (commute them to the end) 
    that were not able to participate in any such grouping such that those don't additionally contaminate insertions
"""

    compiled_circ = Vector{AbstractOperation}()

    remaining_gates = [gate for gate in circuit if typeof(gate) <: AbstractTwoQubitOperator ]

    append!(compiled_circ, [gate for gate in circuit if typeof(gate) <: AbstractSingleQubitOperator])

    # In the above construction of encoding circuits, we can safely push the single-qubit gates to the beginning of the circuit, and since we are working with CSS codes,
    # we only need to consider CX gates. Also, controls here are never targets, so we can safely assume that no gates are in between, and disregard order entirely (everything commutes) for the search (not the insertion, see above)
    
    singular_gates = Vector{AbstractOperation}()

    insertion_gate_qubits = Vector{Int64}()

    # we look for two pairs of CNOTs, as descibed above
    @label while_loop
    while !isempty(remaining_gates)
        # Start with the first gate ik
        CX_ik = remaining_gates[1]
        i, k = affectedqubits(CX_ik)
        # if (k ∈ insertion_gate_qubits)
        #     push!(singular_gates, sCNOT(i,k) )
        #     print(i,k, remaining_gates, CX_ik)
        #     filter!(g -> g != CX_ik , remaining_gates)
        #     @goto while_loop
        # end
        
        # Find all potential jk gates that share the same target
        CX_jks = filter(gate -> (affectedqubits(gate)[1] != i && affectedqubits(gate)[2] == k), remaining_gates) #findall(gate, gate.q2 ==k)

     
        # traverse the list
        if !isempty(CX_jks)
            for CX_jk in CX_jks
                # we have a combination ik + jk, now we need im and jm
                j, k2 = affectedqubits(CX_jk)
                @assert k == k2 "Oops, something went wrong in the filtering, the indices k=$k and k2=$k2 should match"

                CX_ims = filter(gate -> (affectedqubits(gate)[1] == i && affectedqubits(gate)[2] != k), remaining_gates) #findall(gate, gate.q1 ==i)

                if isempty(CX_ims)
                    continue
                else 
                    for CX_im in CX_ims
                        i2, m = affectedqubits(CX_im)
                        
                        @assert i == i2 "Oops, something went wrong in the filtering, the indices i=$i and i2=$i2 should match"

                        # if (m ∈ insertion_gate_qubits)
                        #     push!(singular_gates, sCNOT(i,m) )
                        #     filter!(g -> g != CX_im , remaining_gates)
                        #     continue
                        # end
                        # now we have ik + jk + im > only missing jm
                        # there can only be one such CX_km gate by the way the Gottesman circuit is constructed
                        CX_jms = filter(gate -> (affectedqubits(gate)[1] == j && affectedqubits(gate)[2] == m), remaining_gates) #findfirst(gate, gate.q1 ==j, gate.q2 ==m)
                        @assert length(CX_jms) <= 1 "For any valid choice of i=$i, j=$j k=$k and m=$m, there should at most be one gate available"
                        if isempty(CX_jms)
                            continue
                            #push!(compiled_circ, CX_ik)
                            #remaining_gates.pop(CX_ik) 
                        else
                            # now we have a quartet ik + jk + im + jm
                            # ----  DQC criterion  ------
                            # we check the DQC criterion: are k and m in differnet cores?
                            CX_jm = CX_jms[1] # filter returned an array of length 1
                            k_register = network_specs.register_lookup_array[network_specs.inv_map[k]]
                            m_register = network_specs.register_lookup_array[network_specs.inv_map[m]]
            
                            if k_register != m_register
                                i_register = network_specs.register_lookup_array[network_specs.inv_map[i]]
                                j_register = network_specs.register_lookup_array[network_specs.inv_map[j]]
                                if i_register == m_register && j_register == m_register
                                    # in this special case, i and j are in the same register with m, but k is in a different one
                                    # then we want to eliminate the first two gates, while adding a gate m>k (instead of k>m)

                                    # But first, we need to check the contamination criteria again, on control qubit m
                                    if (m ∈ insertion_gate_qubits)
                                        continue
                                    else
                                        push!(compiled_circ, sCNOT(i,m))
                                        push!(compiled_circ, sCNOT(j,m))
                                        push!(compiled_circ, sCNOT(m,k))
                                        filter!(g -> g ∉ (CX_ik, CX_jk, CX_im, CX_jm) , remaining_gates)
                                        @info "Deleted gates $i>$k, $j>$k, $i>$m and $j>$m, added gates $i>$m and $j>$m, $m>$k "
                                        push!(insertion_gate_qubits, m)
                                        push!(insertion_gate_qubits, k)
                                        # remaining_gates.pop(CX_ik) 
                                        # remaining_gates.pop(CX_jk) 
                                        # remaining_gates.pop(CX_im)
                                        # remaining_gates.pop(CX_jm)  
                                        @goto while_loop 
                                    end
                                end
                            end
                            # in this case (we don't add a telegate, or we don't add an additional one at least), we proceed to add the gate and delete the other two
                            # in this case, it is always beneficial to add the new gate, since we are not introducing a telegate
                            # this strategy might not be globally optimally, but it provides a good greedy approach

                            # But first, we need to check the contamination criteria again
                            if (k ∈ insertion_gate_qubits)
                                continue
                            else
                                push!(compiled_circ, sCNOT(i,k))
                                push!(compiled_circ, sCNOT(j,k))
                                push!(compiled_circ, sCNOT(k,m))
                                filter!(g -> g ∉ (CX_ik, CX_jk, CX_im, CX_jm) , remaining_gates)
                                @info "Deleted gates $i>$k, $j>$k, $i>$m and $j>$m, added gates $i>$k and $j>$k, $k>$m "
                                push!(insertion_gate_qubits, k)
                                push!(insertion_gate_qubits, m)
                                # After adding the 3 gates to circuit and deleting 4 gates, we don't traverse the circuit again (one compilation pass)

                                @goto while_loop 
                            end
                        end
                    end
                end
            end
        end
        
        #Now after all the checks and then nothing was found), we remove this gate from the list while adding it to the circuit to maintain its action
        push!(singular_gates, sCNOT(i,k))
        filter!(g -> g != CX_ik, remaining_gates) #remaining_gates.pop(CX_ik) 

    end

    compiled_circ = vcat(compiled_circ, singular_gates)

    gate_counts = [0,0,0]

    for op in compiled_circ
        
        T = typeof(op)
        if T <: AbstractSingleQubitOperator
            gate_counts += [1,0,0]
        elseif T<: AbstractTwoQubitOperator
            control = op.q1
            target = op.q2
            # perm is the permutation that identifies the qubit to apply the gate to (from permuted encoding circuit),
            # but to identify telegates, we need the qubit mapping (here: the inverse mapping) from the DQC network specifications
            if network_specs.register_lookup_array[network_specs.inv_map[control]] == network_specs.register_lookup_array[network_specs.inv_map[target]]
                # the above is equivalent to network_specs.register_lookup_array[findfirst(==(p_target), network_specs.permutation)]
                gate_counts += [0,1,0]
            else
                gate_counts += [0,0,1]
            end
        end
    end
      

        
    # Then verify correctness!


    return compiled_circ, gate_counts


end







end