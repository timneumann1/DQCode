# helper.jl

"""
Helper functions that performs networking configuration computations,
gate counting, data type conversions and circuit execution.
"""
module Helper

export create_lookup_array, perm_to_transpositions, transpositions_to_perm, data_qubit_partitioning
export tableau_distance, tableau_to_bitmatrix
export gate_counts, verify_success
export qc_circuit_to_qasm
export save_circuit_diagram, code_dirname, save_txt

using ..Types

using QuantumClifford
using QECCore
using KaHyPar
using SparseArrays
using Quantikz: savecircuit, @with, classicalbitslayout
using CairoMakie

import Quantikz: QuantikzOp, ClassicalDecision


"""
    create_lookup_array(num_data_qubits_per_register::Vector{Int})::Tuple{Vector{Int}, Int, Int}

Build a register lookup array mapping each physical slot on the DQC architecture
to its assigned QPU register.

### Input

- `qpu_sizes` -- vector whose `i`-th entry indicates the number of available data qubit slots in register `i`

### Output

Returns a 3-element tuple containing the register lookup array the total number of data qubits,
and the number of communication qubits per register.

### Examples
For the Shor-[[9,1,3]] code, we can define three QPUs with three data qubits each. Then 
`create_lookup_array` returns `([1, 1, 1, 2, 2, 2, 3, 3, 3], 9, 2)`, where the first, middle 
and last three slots in the Type-II architecture belong to the first, second and third core,
respectively, and where have nine memory qubits (on which the QEC code is defined) as well as 
two communication qubits (one for each register).
"""
function create_lookup_array(qpu_sizes::Vector{Int})::Tuple{Vector{Int}, Int, Int}
    num_data_qubits = sum(qpu_sizes)
    num_registers = length(qpu_sizes)
    num_comm_qubits_per_register = num_registers-1
    register_lookup_array = Vector{Int}(undef, num_data_qubits)
    register = 1
    register_start_index = 1
    for j in qpu_sizes
        register_lookup_array[register_start_index:register_start_index+j-1] .= register 
        register_start_index+=j
        register +=1
    end
    return register_lookup_array, num_data_qubits, num_comm_qubits_per_register
end


"""
    perm_to_transpositions(perm::Vector{Int})::Vector{Tuple{Int, Int}}

Transform a permutation to an equivalent sequence of transpositions (two-element swaps).

### Input

- `perm` -- a permutation vector of length `n` (number of data qubits)

### Output

Returns a vector of index pairs `(i, j)` such that applying the transpositions as SWAPS from 
right to left (in the spirit of mathematical notation) performs the permutation captured by `perm`.

### Examples

We assume that `perm = [4, 1, 3, 2]` encodes that qubit `4` is permuted to slot `1`. Calling 
`perm_to_transpositions` on `perm` results in the vector `[(1,2), (2,4)]`. Applying SWAPS on 
an ordered 4-qubit register from right to left maps qubit `1` to slot `2`, qubit `2` to slot `4`,
qubit `3` to slot `3`, and qubit `4` to slot `1`, as desired.
"""
function perm_to_transpositions(perm::Vector{Int})::Vector{Tuple{Int, Int}}
    n = length(perm)
    transpositions = Tuple{Int, Int}[]
    for i in 1:n
        if perm[i]!=i
            j = findlast(==(i), perm)
            @assert !isnothing(j)
            push!(transpositions, (i, j))
            perm[j] = perm[i]
        end
    end
    return transpositions
end


"""
    transpositions_to_perm(transpositions::Vector{Tuple{Int,Int}}, n::Int)::Vector{Int}

Reconstruct a permutation vector from a sequence of transpositions.

### Input

- `transpositions` -- ordered sequence of index pairs `(i, j)` representing two-element swaps
- `n` -- length of the target permutation

### Output

Returns the permutation vector of length `n` obtained by applying the transpositions from 
right to left to the identity permutation, or equivalently, by applying SWAPS corresponding 
to the transpositions from right to left to an ordered `n` qubit register. 

### Examples

Composing the functions `perm_to_transpositions` and `transpositions_to_perm` yields the input back:
`transpositions_to_perm( perm_to_transpositions( [4,1,3,2] ), 4) = [4,1,3,2]`
"""
function transpositions_to_perm(transpositions::Vector{Tuple{Int,Int}}, n::Int)::Vector{Int}
    perm = collect(1:n)
    for (i, j) in reverse(transpositions)
        perm[i], perm[j] = perm[j], perm[i]
    end
    return perm
end


"""
    data_qubit_partitioning(capacities::Vector{Int}, stabilizers::Matrix)::Vector{Int}

Partition/Map data qubits across registers using hypergraph partitioning to minimise inter-register
communication in distributed QEC cycles, subject to per-register capacity constraints, using the 
Karlsruhe Hypergraph Partitioning (KaHyPar) tool (https://doi.org/10.1137/1.9781611974768.3).

### Input

- `capacities` -- vector of length `k` specifying the number of data qubits each register can hold
- `stabilizers` -- tableau of stabilizer group generators

### Output

Returns a mapping vector of length `num_qubits` where value `j` at index `i` indicates that qubit `j`
is mapped to register slot `i` in a flattened register layout (register 1 first, then register 2, etc.).

### Notes

Partitioning is performed via `KaHyPar` using a connectivity objective. The precise optimisation 
configuration can be accessed in `km1_kKaHyPar_sea20.ini`, which encodes "direct k-way partitioning
optimizing the (connectivity - 1) objective" (cf. https://github.com/kahypar/kahypar). 

### Examples

Consider the stabiliser generators for the Steane-[[7,1,3]] code:
    + ___XXXX
    + _XX__XX
    + X_X_X_X
    + ___ZZZZ
    + _ZZ__ZZ
    + Z_Z_Z_Z

Traversing the stabilisers (rows) of this tableau, we build the incidence vectors `I` and `J`

`I = [4, 5, 6, 7, 2, 3, 6, 7, 1, 3, 5, 7, 4, 5, 6, 7, 2, 3, 6, 7, 1, 3, 5, 7]` 
`J = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6]`, 

where `I` collects the qubits involved in the `p`th hyperedge, where `p` is the counter value 
stored at the respective indices in `J`.

The KaHyPar algorithm then partitions the nodes/qubits based on this information, and returns the 
core assignment `[1, 0, 0, 1, 1, 0, 0]`.

This captures the information that qubits `1, 4, 5` are mapped to core `2`, and qubits `2, 3, 6, 7'
are mapped to core `1`. Since we assume full connectivity of cores, these qubits are then ordered
in the `mapping` vector `mapping = [2, 3, 6, 7, 1, 4, 5]`, such that value `j` at index `i` indicates
that qubit `j` is mapped to physical slot `i`. 
Calling `register_lookup_array` on `i` will then return the core of the slot, which can be used to 
evaluate telegate counts etc.
"""
function data_qubit_partitioning(capacities::Vector{Int}, stabilizers::Stabilizer)::Vector{Int}
    k = length(capacities)
    nqubits = size(stabilizers, 2)
    @assert sum(capacities) == nqubits "Register capacities must sum to code length: $(sum(capacities)) vs. $nqubits"  
    # Build incidence matrix: rows=qubits (vertices), cols=stabilizers (hyperedges)
    I = Int[]
    J = Int[]
    e = 0
    for r in 1:size(stabilizers, 1) # traverse through rows in stabilisers and build up the connectivity graph
        support = Int[]
        for q in 1:nqubits
            x, z = stabilizers[r, q]
            if x || z
                push!(support, q)
            end
        end
        if length(support) >= 2 # ignore trivial or weight-1 stabilisers for the mapping
            e += 1
            for q in support
                push!(I, q)
                push!(J, e)
            end
        end
    end
    if e == 0
        @warn "no stabilizer supports found; returning identity mapping"
        return collect(1:nqubits)
    end
    A = sparse(I, J, ones(Int, length(I)), nqubits, e)
    h = KaHyPar.HyperGraph(A) # apply KaHyPa to partition the graph
    cfg = joinpath(@__DIR__, "km1_kKaHyPar_sea20.ini")
    partition = KaHyPar.partition(h, k; configuration=cfg)
    assignments = partition.+1 # KaHyPar is a Python-based optimiser based on 0 indexing
    block_sizes = [count(==(b), assignments) for b in 1:k]
    if block_sizes != capacities
        error("`KaHyPar`` failed to satisfy the capacity constraint with correct block sizes -- please change the capacities or alter the mapping manually.")
    end
    mapping = Int[]
    for b in 1:k
        append!(mapping, sort(findall(==(b), assignments)))
    end
    @assert length(mapping) == nqubits
    @assert sort(mapping) == collect(1:nqubits)
    @info "DQC qubit mapping determined by `KaHyPar``: $mapping"
    return mapping
end


"""
    qc_circuit_to_qasm(circ::Vector{QuantumClifford.AbstractOperation})::String

Convert a `QuantumClifford` circuit to an OpenQASM 2.0 string.

### Input

- `circ` -- sequence of quantum gates

### Output

Returns a valid OpenQASM 2.0 string with a `qreg` sized to the largest
qubit index found in `circ`.

### Notes
Qubit indices are converted from 1-based (Julia) to 0-based (QASM).
Currently, only the conversion of Hadamard and CNOT gates is supported.

### Examples

The `QuantumClifford.AbstractOperation` circuit 

    `QuantumClifford.AbstractOperation[sHadamard(1), sCNOT(1,2), sHadamard(3), sCNOT(4,2)]`

is converted to the QASM 2.0 string
    
    OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[4];
    h q[0];
    cx q[0],q[1];
    h q[2];
    cx q[3],q[1];
"""
function qc_circuit_to_qasm(circ::Vector{QuantumClifford.AbstractOperation})::String
    max_q = 1
    for g in circ
        if g isa AbstractSingleQubitOperator
            @assert g isa sHadamard "the conversion from QuantumClifford circuits to qasm is currently only supported for Hadamard and CNOT gates"
            if g.q > max_q 
                max_q = g.q
            end
        elseif g isa AbstractTwoQubitOperator
            @assert g isa sCNOT || g isa sZCX
            if max(g.q1, g.q2) > max_q 
                max_q = max(g.q1, g.q2)
            end
        end
    end
    qasm_code = String[ "OPENQASM 2.0;", "include \"qelib1.inc\";" , "qreg q[$max_q];"    ]
    for g in circ
        if g isa AbstractSingleQubitOperator
            push!(qasm_code, "h q[$(g.q-1)];")
        elseif g isa AbstractTwoQubitOperator
            push!(qasm_code, "cx q[$(g.q1-1)],q[$(g.q2-1)];")
        end
    end
    return join(qasm_code, "\n")
end


"""
    tableau_to_bitmatrix(tableau::QuantumClifford.Tableau)::Matrix{Int}

Convert a stabilizer `Tableau` to a matrix representation, where integers encode the Pauli operators.

### Input

- `tableau` -- a `QuantumClifford.Tableau` 

### Output

Returns a matrix of size `(# stabilizers \times # qubits)`, where each Pauli entry is encoded as
`0` (I), `1` (X), `2` (Z), or `3` (Y), and the final column stores the phase of each stabilizer 
row as `0` (+) or `1` (-).

### Examples

The tableau

    + ___ZZZZ
    + ___XXXX
    + __Z_ZZ_
    + _Z__Z_Z
    + _XX__XX
    + Z____ZZ
    + X_X_X_X

is converted to the bitmatrix

    0  0  0  2  2  2  2  0
    0  0  0  1  1  1  1  0
    0  0  2  0  2  2  0  0
    0  2  0  0  2  0  2  0
    0  1  1  0  0  1  1  0
    2  0  0  0  0  2  2  0
    1  0  1  0  1  0  1  0

### Notes

The resulting matrix is evidently not a 'bitmatrix'; yet, we stick to this terminology, alluding to 
https://doi.org/10.1103/gqpr-dgz7, in which a tableau has `2n` columns (separating X and Z part), rather 
than `n` columns, as in the QuantumClifford representation.
"""
function tableau_to_bitmatrix(tableau::QuantumClifford.Tableau{<:AbstractVector{UInt8}, <:AbstractMatrix{<:Unsigned}})::Matrix{Int}
    rows, cols = size(tableau)
    bits = Matrix{Int}(undef, rows, cols+1)
    @inbounds for r in 1:rows
        for c in 1:cols
            x, z = tableau[r, c]      
            bits[r, c] = x + 2*z
        end
        bits[r, cols + 1] = tableau.phases[r] 
    end
    return bits
end


"""
    tableau_distance(matrix::Matrix{Int}, target_matrix::Matrix{Int}; metric::String = "jaccard")::Float64

Compute the normalised distance between two tableau bit-matrices.

### Input

- `matrix` -- the current tableau bit-matrix, as produced by `tableau_to_bitmatrix`
- `target_matrix` -- the target tableau bit-matrix to compare against
- `metric` -- (optional, default: `jaccard`) distance metric to use; either `"hamming"` (fraction of differing entries) 
                            or `"jaccard"` (fraction of differing entries restricted to the joint non-identity support)

### Output

Returns a scalar in `[0.0, 1.0]`, where `0.0` indicates identical tableaus.

### Notes

There is an implicit guardrail against a division-by-zero error for the `jaccard`-distance; it can never
occur, since otherwise both `matrix` and `target_matrix` would have to be trivial, and thus not define a
valid quantum state.
"""
function tableau_distance(matrix::Matrix{Int}, target_matrix::Matrix{Int}; metric = "jaccard")::Float64
    @assert size(matrix) == size(target_matrix) "Mismatch in dimensions"
    difference_mask = matrix .!= target_matrix
    if metric == "hamming"
        return count(difference_mask) / length(matrix) 
    elseif metric == "jaccard"
        support_mask = (matrix .!= 0) .| (target_matrix .!= 0) 
        return count(difference_mask .& support_mask) / count(support_mask) # Jaccard distance = 1-intersection / union
    end
end


"""
    gate_counts(circuit::Vector{AbstractOperation}, n::NetworkSpecifications)::Vector{Int}

Count single-qubit gates, intra-regiser/local two-qubit gates, and
inter-register/tele- two-qubit gates in an encoding circuit.

### Input

- `circuit` -- sequence of `AbstractOperation` gates constituting the encoding circuit
- `n` -- network specification object providing `inv_map` and `register_lookup_array`

### Output

Returns a 3-element integer vector `[num_single, num_local_two_qubit, num_telegates]`.

### Notes

`sSWAP` gates are decomposed into three CNOTs when counting, since we can decompose 
`sSWAP(i,j) = sCNOT(i,j)sCNOT(j,i)sCNOT(i,j)`.
"""
function gate_counts(circuit::Vector{QuantumClifford.AbstractOperation}, n::NetworkSpecifications)::Vector{Int}    
    mapping = copy(n.inv_map)
    gate_counts = [0,0,0]
    for op in circuit
        T = typeof(op)
        if T <: AbstractSingleQubitOperator
            gate_counts += [1,0,0]
        elseif T<: AbstractTwoQubitOperator
            control_register = n.register_lookup_array[mapping[op.q1]] 
            target_register = n.register_lookup_array[mapping[op.q2]]
            if control_register == target_register
                if T == sSWAP
                    gate_counts += [0,3,0] # decompose SWAP into three CNOT gates
                else 
                    gate_counts += [0,1,0]
                end
            else
                if T == sSWAP
                    gate_counts += [0,0,3] # decompose SWAP into three non-local CNOT gates
                else 
                    gate_counts += [0,0,1]
                end
            end
        else
            error("Unsupported gate type in gates_to_circuit: $T")
        end
    end
    return gate_counts
end


"""
    verify_success(circuit::Vector{AbstractOperation}, target_state::Stabilizer, n::NetworkSpecifications)::Bool

Verify whether a given circuit successfully prepares the target logical zero state.

### Input

- `circuit` -- quantum circuit to verify
- `target_state` -- the target stabilizer state to check against
- `n` -- network specification object 

### Output

Returns `true` if the circuit encodes the target logical zero state, `false` otherwise.

### Notes

Verification appends a `VerifyOp` to the circuit and uses `mctrajectory!` for simulation, and
throws an `ErrorException` if the simulation terminates with an unexpected status.
"""
function verify_success(circuit::Vector{QuantumClifford.AbstractOperation}, target_state::Stabilizer, n::NetworkSpecifications)::Bool
    verification_circuit = copy(circuit)
    push!(verification_circuit, VerifyOp(target_state, n.data_qubits))
    initial_state = Register(one(MixedDestabilizer, n.num_data_and_comm_qubits),n.num_registers*(n.num_registers-1))
    state, stat = mctrajectory!(initial_state, verification_circuit)
    if stat == true_success_stat
        return true
    elseif stat == false_success_stat
        return false
    else
        throw(ErrorException("Run was invalid due to status: $stat"))
    end
end


"""
    code_dirname(code)::String

Extract a normalised base name for a QEC code for usage as directory name.

### Input

- `code` -- a QEC code object whose type name encodes the code family

### Output

Returns the directory name string corresponding to the code.

### Examples
We remove trailing digits and any `Via...` suffix in the code name, e.g., 
`Steane7` becomes `Steane`, and `TrivariateBicycleViaCirculantMat` becomes `TrivariateBicycle`.
"""
function code_dirname(code::AbstractCSSCode)::String
    name = string(nameof(typeof(code))) 
    name = replace(name, r"Via.*$" => "") 
    name = replace(name, r"\d+$" => "")  
    return name
end


"""
    save_txt(folder::String, title::String, obj::Any)

Serialise all fields of a struct to a plain-text file.

### Input

- `folder` -- path to the output directory
- `title` -- filename for the output file
- `obj` -- any Julia struct whose fields are to be written

### Notes

Each line of the output file has the form `fieldname = repr(value)`.
"""
function save_txt(folder::String, title::String, obj::Any)
    open(joinpath(folder, title), "w") do io
        for fn in fieldnames(typeof(obj))
            println(io, fn, " = ", repr(getfield(obj, fn)))
        end
    end
end


"""
    save_circuit_diagram(circuit::Vector{QuantumClifford.AbstractOperation}, directory::String, label::String)

Render and save a circuit diagram to disk using the `Quantikz`
library (https://arxiv.org/abs/1809.03842, https://github.com/QuantumSavory/Quantikz.jl).

### Input

- `circuit` -- the quantum circuit to render
- `directory` -- path to the output directory
- `label` -- filename (without path) for the saved diagram

### Notes

For moderate register sizes (~ < 15 qubits), the saving works well, yet for larger registers we Quantikz 
memory capacity might be exceeded, in which case we skip the saving.
"""
function save_circuit_diagram(circuit::Vector{QuantumClifford.AbstractOperation}, directory::String, label::String)
    try
        @with classicalbitslayout => :expanded begin   
            savecircuit(
                circuit,
                joinpath(directory, label);
                scale = 1
            )
        end
    catch e
        @warn "saving circuit in $label failed" exception=e
    end
end


function QuantikzOp(op::ConditionalGate)
    # Maps a ConditionalGate to a Quantikz ClassicalDecision for circuit diagram rendering.
    targets = collect(affectedqubits(op.truegate))
    label = _conditional_gate_label(op.truegate)
    return ClassicalDecision(label, targets, op.controlbit)
end


function _conditional_gate_label(g::QuantumClifford.AbstractOperation)
    # Returns the single-character Pauli/gate label for a conditional gate operand.
    repr_g = string(g)
    if occursin("sX", repr_g)
        return "X"
    elseif occursin("sY", repr_g)
        return "Y"
    elseif occursin("sZ", repr_g)
        return "Z"
    elseif occursin("sHadamard", repr_g)
        return "H"
    elseif occursin("sId1", repr_g)
        return "I"
    end
    return "U"
end



end