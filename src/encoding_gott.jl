# encoding_gott.jl

""" 
Encoding circuits for stabiliser codes.

Credit: The implementation of `_gottesman_encoding_circuit_raw` is (up to minor details) identical to the implementation in the 
        QuantumClifford library (https://github.com/QuantumSavory/QuantumClifford.jl, https://arxiv.org/abs/2512.16752), which can
        be inspected at https://github.com/QuantumSavory/QuantumClifford.jl/blob/master/src/ecc/circuits.jl.
"""
module EncodingGott

export encoding_circuit_gott, encoding_gott, overlap_compilation

using ..Types
using ..Helper

using QuantumClifford
using QECCore
using Serialization
using CSV, DataFrames


"""
    _gottesman_encoding_circuit_raw(code::AbstractCSSCode, basis_state::Vector{Int})::Tuple{Vector{QuantumClifford.AbstractOperation}, Vector{Int}}

Accept a CSS code and a computational basis state, and output the Gottesman encoding circuit as originally described by Daniel Gottesman. 

### Input

- `code` -- the CSS code for which the encoding circuit is to be constructed
- `basis_state` -- the computational basis state that is to be encoded

### Output

Tuple of Gottesman encoding circuit based on the permuted qubit space (`circ`) and corresponding permutation of qubits (`gottesman_perm`).

### Notes

In this function, we first obtain the canonical form of the stabiliser tableau, without re-permuting the qubits, while reporting the permutation
operation (see `canonicalize_gott!`). Next, the logical operators are extracted from this canonical form, and the encoding circuit is built. 
While we expose the functionality of encoding an arbitrary computational basis state (labeled `Logical Xs` in the code below), this is not necessary
for encoding the logical zero state of a CSS code. Noteworthy, applying all corresponding `sCNOT` operation irregardless of `basis_state`, we can encode
a general state on the last `k` qubits into the `n` qubit register. However, such encoding is not fault-tolerant. 
    
The Julia implementation of this algorithm, on which the below implementation is based (see Credit in the header), 
is based on [cleve1997efficient] and [gottesman1997stabilizer]. However, the description in these sources was partly erroneous, as noted and 
corrected in [grassl2002algorithmic] and [grassl2011variations] (also cf. https://perimeterinstitute.ca/personal/dgottesman/thesis-errata.html).

### Examples

The `MixedDestabilizer` operation permutes qubits to reach canonical form and returns the permutation. The X- and Z-part of the tableau are canonicalised 
separately, but since they only yield the correct description of the state together, any permutation on one part of the tableau is also applied on the other.
For example, to reach the canonical form of X-part of a 4-qubit system, one might have to apply a column/qubit swap of columns 2 and 3, resulting
in permx = [1,3,2,4] ("qubit 3 is swapped into tableau column 2, and qubit 2 is swapped into tableau column 3"); also, one might have to apply further swaps in 
the Z-part of the tableau, resulting in permz = [1,2,4,3] ("qubit 3 is swapped into tableau column 4, and qubit 4 is swapped into tableau column 3"). Since both 
permutations are applied on the entire tableau, and permx is applied first, this yields the overall permutation [1,2,3,4][permx][permz] = [1,3,2,4][permz] = [1,3,4,2].
Now [1,3,4,2] = [1,2,3,4][permx][permz] = [1,2,3,4][permx[permz]] = [1,2,3,4][perm], so we define gottesman_perm = permx[permz] = [1,3,4,2]. The vector `gottesman_perm` 
then captures the permutation information that qubit 3 has been placed at position 2, qubit 4 at position 3, and qubit 2 at position 4 during the swap operations.
"""
function _gottesman_encoding_circuit_raw(code::AbstractCSSCode, basis_state::Vector{Int})::Tuple{Vector{QuantumClifford.AbstractOperation}, Vector{Int}}
    n = code_n(code)
    k = code_k(code)
    @assert length(basis_state) == k
    code_standard_form, r, permx, permz = MixedDestabilizer(code, undoperm=false, reportperm=true) # setting `undoperm=false` controls that we obtain the canonical Gottesman tableau 
    T = stabilizerview(code_standard_form)
    X = logicalxview(code_standard_form)  
    circ = Vector{QuantumClifford.AbstractOperation}()
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
    gottesman_perm = permx[permz] # Collecting the qubit permutations performed        
    return circ, gottesman_perm
end


"""
    encoding_circuit_gott(qec_code::AbstractCSSCode, network_specs::NetworkSpecifications, basis_state::Vector{Int})::Tuple{Vector{AbstractOperation}, Vector{Int}}

Retrieve the Gottesman encoding circuit on permuted columns/qubits and construct the logical state encoding cicuit using the original
qubit ordering defined in `qec_code`. While constructing the circuit, obtain the gate count for the encoding circuit.

### Input

- `qec_code` -- the CSS code for which the encoding circuit is to be constructed
- `network_specs` -- Type-II DQC networking specifications
- `basis_state` -- the computational basis state that is to be encoded

### Output

Returns the Gottesman encoding circuit using the original qubit ordering as well as the corresponding gate count, collecting 
single-, two-qubit- and telegates.

### Notes

The resulting circuit encodes the specified basis state on physical qubits at indices 'n-k+1:n' into a logical code state using a non-fault-tolerant 
encoding circuit. If we start from the all-zero physical register, the resulting circuit encodes the logical all-zero state.

We wish to recover the original ordering from the canonicalization permutation. From the Notes section of `_gottesman_encoding_circuit_raw`, we know
that `perm` here captures which qubit `j` has been placed at tableau position `i` if the `i`th index of the list is `j`. 
Thus, for a quantum operation `op(i)`, to obtain the original labeling of columns/qubits, we wish to apply the operation on qubit `j`. This can
be achieved by applying `op` on qubit `perm[i] = j`, i.e. to apply op on `perm[op.q]` instead of `op.q`.

In order to retrieve the gate counts, we make use of the `register_lookup_array` as well as the `inv_map` arrays stored in `NetworkSpecifications`.
Since the value `b` at index `a` in `network_specs.mapping` indicates that qubit `b` has been mapped to physical slot `a` in the DQC architecture, and 
`register_lookup_array` collects the registers for each physical slot, we are interested in the value of `register_lookup_array` at index `a` (which 
is populated by qubit `b`). Thus, if an operation is applied to qubit `b` (in the notation of above, `b = perm[op.q]`), we extract its physical slot `a`
by using `register_lookup_array[inv_map[b]]`, where `inv_map` captures the inverse permutation of `mapping`; `inv_map[b]` evaluates to `a`, which is the
correct index to retrieve from `register_lookup_array`.
"""
function encoding_circuit_gott(qec_code::AbstractCSSCode, network_specs::NetworkSpecifications, basis_state::Vector{Int})::Tuple{Vector{QuantumClifford.AbstractOperation}, Vector{Int}}
    @assert qec_code !== nothing 
    encoding_circ_gott, perm = _gottesman_encoding_circuit_raw(qec_code, basis_state)
    encoding_circuit_orig = Vector{QuantumClifford.AbstractOperation}()
    gate_counts = [0,0,0]
    for op in encoding_circ_gott
        T = typeof(op)
        if T <: QuantumClifford.AbstractSingleQubitOperator
            push!(encoding_circuit_orig, T(perm[op.q]))
            gate_counts[1] += 1
        elseif T<: QuantumClifford.AbstractTwoQubitOperator
            p_control = perm[op.q1]
            p_target = perm[op.q2]
            push!(encoding_circuit_orig, T(p_control, p_target))
            if network_specs.register_lookup_array[network_specs.inv_map[p_control]] == network_specs.register_lookup_array[network_specs.inv_map[p_target]]
                gate_counts[2] += 1
            else
                gate_counts[3] += 1
            end
        end
    end
    circuit_sizes = gate_counts
    @assert length(encoding_circuit_orig) == sum(circuit_sizes) 
    return encoding_circuit_orig, gate_counts
end


"""
    overlap_compilation(circuit::Vector{QuantumClifford.AbstractOperation}, 
                        network_specs::NetworkSpecifications)::Tuple{Vector{QuantumClifford.AbstractOperation}, Vector{Int}}

Apply a DQC circuit compilation method based on the overlap method introduced by Paetznick and 
Reichardt (https://arxiv.org/abs/1106.2190) to reduce the number of telegates in the encoding circuit. 

### Input

- `circuit` -- the initial uncompiled Gottesman encoding circuit
- `network_specs` -- Type-II DQC networking specifications

### Output

Tuple of the compiled encoding circuit (`compiled_circ`) and the improved gate counts as a 3-element vector (`gate_counts`).

### Notes

The aforementioned method makes use of overlapping CNOT gates to reduce certain sequences of four CNOTs to three. The compilation 
method is specifically tailored to optimise the logical zero state encoding of CSS codes. These restrictions ensure that control qubits 
in `circuit` are never targets, and they further allow us to assume a H-CNOT template (possibly padded with X- and/or Z-gates for
negative phase correction in the end). Therefore, we append the Hadamard gates in `circuit` first, then perform the compilation on
the CNOT layer, and lastly append the single-qubit gate corrections. 

For details about the compilation method itself, please refer to the publication that this repository accompanies.
"""
function overlap_compilation(circuit::Vector{QuantumClifford.AbstractOperation}, 
                                network_specs::NetworkSpecifications)::Tuple{Vector{QuantumClifford.AbstractOperation}, Vector{Int}}
    compiled_circ = Vector{QuantumClifford.AbstractOperation}()
    remaining_gates = Vector{QuantumClifford.AbstractOperation}()
    phase_correction_gates = Vector{QuantumClifford.AbstractOperation}()
    for gate in circuit
        if typeof(gate) <: QuantumClifford.AbstractSingleQubitOperator
            if gate isa sHadamard
                push!(compiled_circ, gate)
            else
                push!(phase_correction_gates, gate)
            end
        elseif typeof(gate) <: QuantumClifford.AbstractTwoQubitOperator
            if gate isa sCNOT || gate isa sZCX
                push!(remaining_gates, gate)
            else
                @error "oops, this gate was not expected in the encoding circuit -- please double-check your code definition or submit an issue"
            end
        end
    end
    singular_gates = Vector{QuantumClifford.AbstractOperation}() # will accumulate gates that cannot be grouped with others in compilation
    insertion_gate_qubits = Vector{Int64}()
    @label while_loop
    while !isempty(remaining_gates)
        CX_ik = remaining_gates[1] # start with the first gate `ik`
        i, k = affectedqubits(CX_ik)   
        # find all potential `jk` gates that share the same target 
        CX_jks = filter(gate -> (affectedqubits(gate)[1] != i && affectedqubits(gate)[2] == k), remaining_gates) 
        if !isempty(CX_jks)
            for CX_jk in CX_jks
                j, k2 = affectedqubits(CX_jk) # for `{ik, jk}`, we need to find `im` and `jm`
                @assert k == k2 "Oops, something went wrong in the filtering, the indices k=$k and k2=$k2 should match"
                CX_ims = filter(gate -> (affectedqubits(gate)[1] == i && affectedqubits(gate)[2] != k), remaining_gates)
                if isempty(CX_ims)
                    continue
                else 
                    for CX_im in CX_ims
                        i2, m = affectedqubits(CX_im)
                        @assert i == i2 "Oops, something went wrong in the filtering, the indices i=$i and i2=$i2 should match"
                        # for `{ik, jk, im}`, we need to find unique `jm` 
                        CX_jms = filter(gate -> (affectedqubits(gate)[1] == j && affectedqubits(gate)[2] == m), remaining_gates)
                        @assert length(CX_jms) <= 1 "for any valid choice of i=$i, j=$j k=$k and m=$m, the circuit can contain at most one CNOT_{jm}"
                        if isempty(CX_jms)
                            continue
                        else
                            # -------  DQC compilation for quartet `{ik, jk, im, jm}` ---------
                            CX_jm = CX_jms[1] # filter() returned an array of length 1
                            k_register = network_specs.register_lookup_array[network_specs.inv_map[k]]
                            m_register = network_specs.register_lookup_array[network_specs.inv_map[m]]
                            # determine whether k and m are in different cores
                            if k_register != m_register
                                i_register = network_specs.register_lookup_array[network_specs.inv_map[i]]
                                j_register = network_specs.register_lookup_array[network_specs.inv_map[j]]
                                if i_register == m_register && j_register == m_register
                                    if (m ∈ insertion_gate_qubits) # verify contamination criteria on m
                                        continue
                                    else
                                        # `{i,j,m}` share a core disjoint from `k` -> eliminate `ik` and `jk` and add `mk`
                                        push!(compiled_circ, sCNOT(i,m))
                                        push!(compiled_circ, sCNOT(j,m))
                                        push!(compiled_circ, sCNOT(m,k))
                                        filter!(g -> g ∉ (CX_ik, CX_jk, CX_im, CX_jm) , remaining_gates)
                                        push!(insertion_gate_qubits, m)
                                        push!(insertion_gate_qubits, k) 
                                        @goto while_loop 
                                    end
                                end
                            end
                            if (k ∈ insertion_gate_qubits) # verify contamination criteria on k
                                continue
                            else
                                # eliminate `im` and `jm` and add `km`
                                push!(compiled_circ, sCNOT(i,k))
                                push!(compiled_circ, sCNOT(j,k))
                                push!(compiled_circ, sCNOT(k,m))
                                filter!(g -> g ∉ (CX_ik, CX_jk, CX_im, CX_jm) , remaining_gates)
                                push!(insertion_gate_qubits, k)
                                push!(insertion_gate_qubits, m)
                                @goto while_loop 
                            end
                        end
                    end
                end
            end
        end        
        push!(singular_gates, sCNOT(i,k))
        filter!(g -> g != CX_ik, remaining_gates)
    end
    compiled_circ = vcat(compiled_circ, singular_gates,phase_correction_gates)
    # retrieve the update gate counts
    gate_counts = [0,0,0]
    for op in compiled_circ
        T = typeof(op)
        if T <: QuantumClifford.AbstractSingleQubitOperator
            gate_counts += [1,0,0]
        elseif T<: QuantumClifford.AbstractTwoQubitOperator
            control = op.q1
            target = op.q2
            if network_specs.register_lookup_array[network_specs.inv_map[control]] == network_specs.register_lookup_array[network_specs.inv_map[target]]
                gate_counts += [0,1,0]
            else
                gate_counts += [0,0,1]
            end
        end
    end
    return compiled_circ, gate_counts
end


"""
    encoding_gott(code_params::CodeParameters, network_specs::NetworkSpecifications, 
                    basis_state::Vector{Int})::Tuple{Vector{AbstractOperation}, Vector{AbstractOperation}, Bool, Bool, Vector{Int}, Vector{Int}}

Entry point to module used to initialise Gottesman circuit generation and DQC compilation.

### Input

- `code_params` -- code parameters
- `network_specs` -- Type-II DQC networking specifications
- `basis_state` -- the computational basis state to be encoded

### Output

Returns a 6-element Tuple containin the Gottesman encoding circuit (on original qubit ordering), the overlap-compiled
encoding circuit, the verification booleans for both circuits, and the gate counts for both circuits.
"""
function encoding_gott(code_params::CodeParameters, network_specs::NetworkSpecifications, 
                        basis_state::Vector{Int})::Tuple{Vector{QuantumClifford.AbstractOperation}, Vector{QuantumClifford.AbstractOperation}, Bool, Bool, Vector{Int}, Vector{Int}}
    encoding_circ, gate_counts = encoding_circuit_gott(code_params.qec_code, network_specs, basis_state)
    encoding_circ_compiled, gate_counts_compiled = overlap_compilation(encoding_circ, network_specs)
    verification_logical_state = verify_success(encoding_circ, code_params.target_state, network_specs)
    @info "Verification of Gottesman circuit successful: $verification_logical_state, Gate count: $gate_counts"
    verification_logical_state_compiled = verify_success(encoding_circ_compiled, code_params.target_state, network_specs)
    @info "Verification of Compiled circuit successful: $verification_logical_state_compiled, Gate count: $gate_counts_compiled"
    return encoding_circ, encoding_circ_compiled, verification_logical_state, 
            verification_logical_state_compiled, gate_counts, gate_counts_compiled
end


end