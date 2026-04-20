module QECTools

# from https://github.com/QuantumSavory/QuantumClifford.jl/blob/master/src/ecc/decoder_pipeline.jl

using ..Types
using ..Helper

using QuantumClifford
using QuantumClifford: AbstractOperation
using QECCore
using QuantumClifford.ECC: AbstractSyndromeDecoder, faults_matrix
using Combinatorics: combinations

import QuantumClifford.ECC: parity_checks
export parity_checks

export ClassicalTableDecoder, CSSTableDecoder, decode

function perfect_ancillary_paulimeasurement(p::PauliOperator, ancillary_index, bit_index, network_specs)
    circuit = AbstractOperation[]
    num_data_qubits = nqubits(p)
    @assert num_data_qubits == network_specs.num_data_qubits
    for qubit in 1:num_data_qubits
        if p[qubit] == (1,0)
            push!(circuit, sXCX(network_specs.comm_idx[qubit], ancillary_index)) # X-controlled X
        elseif p[qubit] == (0,1)
            push!(circuit, sCNOT(network_specs.comm_idx[qubit], ancillary_index)) # Z-controlled X
        elseif p[qubit] == (1,1)
            push!(circuit, sYCX(network_specs.comm_idx[qubit], ancillary_index)) # Y-controlled X
        end
    end
    p.phase[] == 0 || push!(circuit, sX(ancillary_index))
    mz = sMRZ(ancillary_index, bit_index)
    push!(circuit, mz)

    return circuit
end

function syndrome_circuit(parity_check_tableau, ancillary_index, bit_index, network_specs)
    syndrome_circ = AbstractOperation[]
    ancillaries = 0
    bits = 0
    for check in parity_check_tableau
        append!(syndrome_circ, perfect_ancillary_paulimeasurement(check, ancillary_index+ancillaries, bit_index+bits, network_specs))
        ancillaries +=1
        bits +=1
    end

    return syndrome_circ, ancillaries, bit_index:bit_index+bits-1
end

# function physical_ECC_circuit(H, setup::NaiveSyndromeECCSetup)
#     syndrome_circ, n_anc, syndrome_bits = syndrome_circuit(H)
#     return syndrome_circ, syndrome_bits, n_anc
# end

# function evaluate_decoder(d::AbstractSyndromeDecoder, setup::AbstractECCSetup, nsamples::Int)
#     H = parity_checks(d)
#     n = code_n(H)
#     k = code_k(H)
#     O = faults_matrix(H)

#     physical_noisy_circ, syndrome_bits, n_anc = physical_ECC_circuit(H, setup)
#     encoding_circ = naive_encoding_circuit(H)
#     preX = sHadamard[sHadamard(i) for i in n-k+1:n]

#     mdH = MixedDestabilizer(H)
#     logX_circ, _, logX_bits = syndrome_circuit(logicalxview(mdH), n_anc+1, last(syndrome_bits)+1)
#     logZ_circ, _, logZ_bits = syndrome_circuit(logicalzview(mdH), n_anc+1, last(syndrome_bits)+1)

#     # Evaluate the probability for X logical error (the Z-observable part of the faults matrix is used)
#     X_error = evaluate_decoder(
#         d, nsamples,
#         vcat(encoding_circ, physical_noisy_circ, logZ_circ),
#         syndrome_bits, logZ_bits, O[end÷2+1:end,:])
#     # Evaluate the probability for Z logical error (the X-observable part of the faults matrix is used)
#     Z_error = evaluate_decoder(
#         d, nsamples,
#         vcat(preX, encoding_circ, physical_noisy_circ, logX_circ),
#         syndrome_bits, logX_bits, O[1:end÷2,:])
#     return (X_error, Z_error)
# end


struct ClassicalTableDecoder <: AbstractSyndromeDecoder
    """Parity check matrix defining the code"""
    H::Matrix{Bool}
    """The number of bits in the code"""
    n::Int
    """The number of parity checks"""
    s::Int
    """The maximum weight of errors in the lookup table"""
    error_weight::Int
    """The lookup table corresponding to the code"""
    lookup_table::Dict{Vector{Bool},Vector{Bool}}
    lookup_buffer::Vector{Bool}
    ClassicalTableDecoder(H, n, s, error_weight, lookup_table) = new(H, n, s, error_weight, lookup_table, fill(false, s))
end

function ClassicalTableDecoder(H::Matrix{Bool}; error_weight=1)
    s, n = size(H)
    lookup_table = create_lookup_table(H; error_weight)
    return ClassicalTableDecoder(H, n, s, error_weight, lookup_table)
end

function create_lookup_table(H::Matrix{Bool}; error_weight=1) # TODO there is inefficient casting between Bool vectors and bitvectors here
    lookup_table = Dict{Vector{Bool},Vector{Bool}}()
    s, n = size(H)
    # Process errors from highest weight to lowest
    # so that lower-weight errors overwrite higher-weight ones
    # (lower-weight errors are more probable)
    for w in error_weight:-1:1
        for positions in combinations(1:n, w)
            error = falses(n)
            for pos in positions
                error[pos] = true
            end
            # Calculate syndrome: s = H * e (mod 2)
            syndrome = Bool[(sum(H[row, pos] for pos in positions) % 2) == 1 for row in 1:s]
            lookup_table[syndrome] = error
        end
    end
    # In the case of no errors
    lookup_table[falses(s)] = falses(n)
    lookup_table
end

function decode(d::ClassicalTableDecoder, syndrome_sample)
    d.lookup_buffer .= syndrome_sample
    return get(d.lookup_table, d.lookup_buffer, nothing)
end


struct CSSTableDecoder <: AbstractSyndromeDecoder
    """Stabilizer tableau defining the code"""
    H
    """Faults matrix corresponding to the code"""
    faults_matrix
    """The number of qubits in the code"""
    n::Int
    """The number of parity checks"""
    s::Int
    """The number of encoded qubits"""
    k::Int
    """The number of X checks"""
    cx::Int
    """The number of Z checks"""
    cz::Int
    """The maximum weight of errors in the lookup table"""
    error_weight::Int
    """Classical decoder for X errors (decodes Z syndrome)"""
    tabledecoderx::ClassicalTableDecoder
    """Classical decoder for Z errors (decodes X syndrome)"""
    tabledecoderz::ClassicalTableDecoder
end

function CSSTableDecoder(c; error_weight=1)
    Hx = parity_matrix_x(c)
    Hz = parity_matrix_z(c)
    H = parity_checks(c) 
    # is this not just Stabilizer(c)?
    #H = Stabilizer(parity_matrix(c))
    s, n = size(H)
    _, _, r = canonicalize!(copy(H), ranks=true)
    k = n - r
    cx = size(Hx, 1)
    cz = size(Hz, 1)
    fm = faults_matrix(H)
    tabledecoderx = ClassicalTableDecoder(Matrix{Bool}(Hx); error_weight)
    tabledecoderz = ClassicalTableDecoder(Matrix{Bool}(Hz); error_weight)
    return CSSTableDecoder(H, fm, n, s, k, cx, cz, error_weight, tabledecoderx, tabledecoderz)
end

parity_checks(d::CSSTableDecoder) = d.H

function decode(d::CSSTableDecoder, syndrome_sample)
    row_x = @view syndrome_sample[1:d.cx]
    row_z = @view syndrome_sample[d.cx+1:d.cx+d.cz]
    guess_z = decode(d.tabledecoderx, row_x)
    guess_x = decode(d.tabledecoderz, row_z)
    return isnothing(guess_x) || isnothing(guess_z) ? nothing : vcat(guess_x, guess_z)
end



end