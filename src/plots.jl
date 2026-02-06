""" Credit: https://github.com/QuantumSavory/QuantumClifford.jl/blob/4d524965a11b6d8594d578d4935ccd7bc385f56c/src/ecc/circuits.jl"""

using QuantumClifford
using QECCore
using Quantikz




function naive_encoding_circuit(code, initial_state; undoperm=true)

    # Creatung the canonical tableau (without re-permuting the columns)
    n = code_n(code)
    k = code_k(code)
    #println("\nFor the given code, we have n=$n, k=$k. \n")
    
    code_standard_form, r, permx, permz = MixedDestabilizer(code, undoperm=false, reportperm=true); # undoperm without returns gives the orignal stabiliser tableau
    #println("Standard form of code is \n$code_standard_form")
    X = logicalxview(code_standard_form)
    Z = logicalzview(code_standard_form)
    
    # Creating te canonical tableau (with re-permuting the columns) for final state specification
    code_original_with_logicals = MixedDestabilizer(code, undoperm=true);
    # X_cleaned = logicalxview(md_cleaned)
    # println("X is $X_cleaned")
    # Z_cleaned = logicalzview(md_cleaned)
    # println("Z is $Z_cleaned")

    circ = QuantumClifford.AbstractOperation[]
    push!(circ, Reset(initial_state, [1,2,3,4,5,6,7]))
    
    # Constructing the encoding circuit

    S = stabilizerview(code_standard_form)

    # logical Xs
    for i in 1:k
        for t in 1:n-k
            if X[i,t][1] == true
                push!(circ, sCNOT(n-k+i, t))
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
                g = if xz == (true, true)  # Y
                    sZCY
                elseif xz == (true, false) # X
                    sZCX
                elseif xz == (false, true) && !(i<t<n-k+1) # Z
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

    # undoing the permutations to have the correct circuit for the final (re-permuted) tableau
    if undoperm
        perm = permx[permz]
        transpositions = perm_to_transpositions(perm)
        for (i,j) in transpositions
            push!(circ, sSWAP(i,j))
        end
    end
    code_original_with_logicals, circ
end

function perm_to_transpositions(perm)
    n = length(perm)
    transpositions = Tuple{Int, Int}[]
    for i in n:-1:1
        if perm[i]!=i
            j = findfirst(==(i), perm)
            @assert !isnothing(j)
            push!(transpositions, (i, j))
            perm[j] = perm[i]
        end
    end
    return transpositions
end

#####################################################################################
#####################################################################################
#####################################################################################

function golay_encoding_circuit()
    """[Fault-tolerant ancilla preparation and noise threshold lower bounds  for the 23-qubit Golay code, 2013]"""
    circ = QuantumClifford.AbstractOperation[]
    push!(circ, sHadamard(1))
    push!(circ, sHadamard(2))
    push!(circ, sHadamard(4))
    push!(circ, sCNOT(1,7))
    push!(circ, sCNOT(2,3))
    push!(circ, sCNOT(4,6))
    push!(circ, sCNOT(1,5))
    push!(circ, sCNOT(2,6))
    push!(circ, sCNOT(1,3))
    push!(circ, sCNOT(4,5))
    push!(circ, sCNOT(6,7))
    circ
end

#####################################################################################
#####################################################################################
#####################################################################################

function plot_zero_to_zero()

    # State to be encoded must be at indices 'n-k+1:n', so for Steane-7: |0>^6 \otimes |psi> at index 7
    initial_zero = S"IIIIIZI IIIIZII IIIZIII IIZIIII IZIIIII ZIIIIII IIIIIIZ" 
    initial_one =  S"IIIIIZI IIIIZII IIIZIII IIZIIII IZIIIII ZIIIIII -IIIIIIZ" 
    initial_plus = S"IIIIIZI IIIIZII IIIZIII IIZIIII IZIIIII ZIIIIII IIIIIIX" 


    encoded_zero = S"XIXIXIX IXXIIXX IIIXXXX ZIZZIZI  ZZIIZZI ZZIZIIZ IZIZIZI" # Z_L = Z_2 Z_4 Z_6, matching the logical operator retrieved with MixedDestabiliser
    encoded_one = S"XIXIXIX IXXIIXX IIIXXXX ZIZZIZI  ZZIIZZI ZZIZIIZ -ZZZZZZZ"
    encoded_plus = S"XIXIXIX IXXIIXX IIIXXXX ZIZZIZI  ZZIIZZI ZZIZIIZ XXXXXXX"

    _, circuit = naive_encoding_circuit(Steane7(), initial_zero)
    push!(circuit, VerifyOp(encoded_zero , [1,2,3,4,5,6,7]))
    savecircuit(circuit, "src/plots/circuits/circuit_encoding_zero.png")
end

function plot_golay()
    savecircuit(golay_encoding_circuit(), "src/plots/circuits/golay_encoding_zero.png")
end

plot_golay()
