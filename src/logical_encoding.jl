""" Credit: https://github.com/QuantumSavory/QuantumClifford.jl/blob/4d524965a11b6d8594d578d4935ccd7bc385f56c/src/ecc/circuits.jl"""

module LogicalEnc

using ..Helper

using QuantumClifford
using QECCore
using Quantikz: savecircuit
using QuantumSavory: H, CNOT, X, Y, Z, stateof
using QuantumClifford: true_success_stat, false_success_stat, continue_stat, failure_stat

export naive_encoding_circuit
#export golay_encoding_circuit
export run_tests

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


function naive_encoding_circuit(code; undoperm=true)

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
    #push!(circ, Reset(initial_state, [1,2,3,4,5,6,7]))
    
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

#####################################################################################
#####################################################################################
#####################################################################################

function golay_encoding_circuit()
    """ From [Fault-tolerant ancilla preparation and noise threshold lower bounds  for the 23-qubit Golay code, 2013], Figure 5a"""
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

function test(;initial_state, circuit_type, hadamards, final_state_verify, plotting, label)
    if circuit_type == "naive"
        # N&C version: S"ZIZIZIZ XIXIXIX IZZIIZZ IXXIIXX IIIZZZZ IIIXXXX" 
        # Steane7(): S"XIXIXIX IXXIIXX IIIXXXX ZIZZIZI  ZZIIZZI ZZIZIIZ" 
        code_original_with_logicals, circuit = naive_encoding_circuit(Steane7())
        #println("Final state tableau: \n$code_original_with_logicals")

    end
    if circuit_type == "golay"
        circuit = golay_encoding_circuit()
    end
    if hadamards == true
        for i in 1:7
            push!(circuit, sHadamard(i))
        end
    end
    push!(circuit, VerifyOp(final_state_verify , [1,2,3,4,5,6,7]))
    if plotting
        savecircuit(circuit, "src/plots/circuits/circuit_$(label).png") # plotting is performed by enabling the reset function
    end
    petrajectories(initial_state, circuit)
end


### TESTS ###

function run_tests()

    # State to be encoded must be at indices 'n-k+1:n', so for Steane-7: |0>^6 \otimes |psi> at index 7
    initial_zero = S"IIIIIZI IIIIZII IIIZIII IIZIIII IZIIIII ZIIIIII IIIIIIZ" 
    initial_one =  S"IIIIIZI IIIIZII IIIZIII IIZIIII IZIIIII ZIIIIII -IIIIIIZ" 
    initial_plus = S"IIIIIZI IIIIZII IIIZIII IIZIIII IZIIIII ZIIIIII IIIIIIX" 


    encoded_zero = S"XIXIXIX IXXIIXX IIIXXXX ZIZZIZI  ZZIIZZI ZZIZIIZ IZIZIZI" # Z_L = Z_2 Z_4 Z_6, but works with ZZZZZZZ as well 
    encoded_one = S"XIXIXIX IXXIIXX IIIXXXX ZIZZIZI  ZZIIZZI ZZIZIIZ -IZIZIZI"
    encoded_plus = S"XIXIXIX IXXIIXX IIIXXXX ZIZZIZI  ZZIIZZI ZZIZIIZ IIXIXXI" # works with XXXXXXX

    test1 = test(initial_state=initial_zero, circuit_type="naive", hadamards=false, final_state_verify=encoded_zero, plotting=true, label="zero_to_logical_zero")
    test2 = test(initial_state=initial_one, circuit_type="naive", hadamards=false, final_state_verify=encoded_zero,  plotting =false, label="")
    test3 = test(initial_state=initial_one, circuit_type="naive", hadamards=false, final_state_verify=encoded_one,  plotting =true, label="one_to_logical_one")
    test4 = test(initial_state=initial_zero, circuit_type="naive", hadamards=true, final_state_verify=encoded_plus, plotting =true, label="zero_to_logical_plus")
    test5 = test(initial_state=initial_zero, circuit_type="naive", hadamards=true, final_state_verify=encoded_zero, plotting =false, label="")

    """
    Expected Test Results:

    Test1: 1.0
    Test2: 0.0
    Test3: 1.0
    Test4: 1.0
    Test5: 0.0
    """

    test6 = test(initial_state=initial_zero, circuit_type="golay", hadamards=false, final_state_verify=encoded_zero, plotting=true, label="zero_to_logical_zero_golay")
    test7 = test(initial_state=initial_one, circuit_type="golay", hadamards=false, final_state_verify=encoded_zero,  plotting =false, label="")
    test8 = test(initial_state=initial_one, circuit_type="golay", hadamards=false, final_state_verify=encoded_one,  plotting =true, label="one_to_logical_one_golay")
    test9 = test(initial_state=initial_zero, circuit_type="golay", hadamards=true, final_state_verify=encoded_plus, plotting =true, label="zero_to_logical_plus_golay")
    test10 = test(initial_state=initial_zero, circuit_type="golay", hadamards=true, final_state_verify=encoded_zero, plotting =false, label="")

    """
    Expected Test Results:
    Test6: 1.0
    Test7: 0.0
    Test8: 0.0
    Test9: 1.0
    Test10: 0.0
    """

    test_superposition = test(initial_state=initial_plus, circuit_type="naive", hadamards=false, final_state_verify=encoded_plus, plotting =true, label="superposition_circuit")

    """
    Expected Test Results: 1.0
    """

    println()
    print("Naive Encoding Test\n")
        
    print("Test 1 successful: $( test1[ true_success_stat ]  )\n") # prints the result of true_success (verification was successful)
    print("Test 2 successful: $( test2[ true_success_stat]  )\n") 
    print("Test 3 successful: $( test3[ true_success_stat ]  )\n") 
    print("Test 4 successful: $( test4[ true_success_stat ]  )\n") 
    print("Test 5 successful: $( test5[ true_success_stat]  )\n") 


    println()
    print("Golay Encoding Test\n")


    print("Test 6 successful: $( test6[ true_success_stat ]  )\n") # prints the result of true_success (verification was successful)
    print("Test 7 successful: $( test7[ true_success_stat ]  )\n") 
    print("Test 8 successful: $( test8[ true_success_stat ]  )\n") 
    print("Test 9 successful: $( test9[ true_success_stat ]  )\n") 
    print("Test 10 successful: $( test10[ true_success_stat ] )\n") 


    # The compressed circuit in [Fault-tolerant ancilla preparation and noise threshold lower bounds  for the 23-qubit Golay code] is specific to the initial all-zero
    # state and therefore does not work for encoding arbitrary input states.

    println()
    print("Superposition Encoding Test\n")

    print("Test successful: $( test_superposition[ true_success_stat ]  )\n")
end

end