""" Credit: https://github.com/QuantumSavory/QuantumClifford.jl/blob/4d524965a11b6d8594d578d4935ccd7bc385f56c/src/ecc/circuits.jl"""

module LogicalEnc

using ..Types
using ..Helper
#using ..CircuitSimulator


using QuantumClifford
using QECCore
using Quantikz: savecircuit
#using QuantumSavory: H, CNOT, X, Y, Z, stateof
using QuantumClifford: true_success_stat, false_success_stat, continue_stat, failure_stat


export standard_logical_zero_encoding_circuit
#export golay_encoding_circuit
export baseline_encoding, baseline_encoding_circuit, run_tests

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


function standard_logical_zero_encoding_circuit(code; undoperm=true, logical_Xs = false)
     
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
    transpositions = nothing
    # undoing the permutations to have the correct circuit for the final (re-permuted) tableau
    if undoperm
        perm = permx[permz]
        transpositions = perm_to_transpositions(perm)
        for (i,j) in transpositions
            push!(circ, sSWAP(i,j))
        end
    end
    code_original_with_logicals, circ, transpositions
end

function baseline_encoding_circuit(qec_code, network_specs, gate_set; logical_Xs = false)

    @assert qec_code !== nothing 

    code_original_with_logicals, baseline_circuit, transpositions = standard_logical_zero_encoding_circuit(qec_code, logical_Xs = logical_Xs)
    #standard_circuit_length = length(standard_circuit[2])
    permutation = transpositions_to_perm(reverse(transpositions), network_specs.num_data_qubits)

    baseline_gates = Gate[]
    for op in baseline_circuit
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
            error("Unsupported gate type: $(typeof(op))")
        end  
        
    end

    #baseline_raw_circuit = CircuitIndividual(baseline_gates)
    #println(typeof(baseline_raw_circuit))
    #baseline_exec_circuit, num_single_qubit_gates, num_two_qubit_gates, num_telegates = construct_executable_circuit(baseline_gates, gate_set, network_specs, telegate_overhead = true)
    baseline_circuit, gate_counts = gates_to_circuit(baseline_gates, network_specs)
    circuit_sizes = gate_counts
    # print(baseline_circuit)
    # print(circuit_sizes)
    # print(circuit_size(baseline_exec_circuit))
    # print(sum(circuit_sizes))
    #print(baseline_gates)
    #println(circuit_sizes)
    @assert length(baseline_gates) == sum(circuit_sizes) 
    #println(typeof(baseline_exec_circuit))
    #println(typeof(baseline_raw_circuit.gates))
    println("Original Code with repermuted logical operators appended: $(code_original_with_logicals)")
    println("\nRaw size of standard encoding circuit: $(sum(circuit_sizes))\n")
    println("DQC Size of standard encoding circuit: $(circuit_sizes)\n")
    #println("Single qubit: $num_single_qubit_gates, Two qubit: $num_two_qubit_gates, Telegates: $num_telegates")
    return baseline_gates, baseline_circuit, gate_counts
end

function baseline_encoding(code_params, network_specs, gate_set; data_storage = true)

    ### Baseline comparison of standard encoding in DQC setting
    baseline_gates, baseline_circuit, gate_counts = baseline_encoding_circuit(code_params.qec_code, network_specs, gate_set)

    verification_logical_state = verify_success(baseline_circuit, code_params.target_state, network_specs)

    println("\nVerification of baseline circuit successful (target state fidelity; only expressive (binary) in noiseless setting): $verification_logical_state\n")

    if data_storage
       
        ###########################################
        ############# DATA STORAGE ################
        ###########################################

        dir = joinpath(@__DIR__, "results", string(code_dirname(code_params.qec_code)), "baseline_encoding")
        mkpath(dir)

        println("Saving results to $(dir)")

        save_circuit_diagram(baseline_gates, dir, "baseline_raw_circuit__size_$(sum(gate_counts)).png")
        #save_circuit_diagram(baseline_exec_circuit, dir, "baseline_exec_circuit__size_$(circuit_sizes).png")

        open(joinpath(dir, "network_specs.txt"), "w") do io
            println(io, "Network Specifications")
            for fname in fieldnames(Types.NetworkSpecifications)
                println(io, fname, " = ", repr(getfield(network_specs, fname)))
            end
        end

        open(joinpath(dir, "code_params.txt"), "w") do io
            println(io, "Code parameters")
            for fname in fieldnames(Types.CodeParameters)
                println(io, fname, " = ", repr(getfield(code_params, fname)))
            end
        end

        open(joinpath(dir, "baseline_raw.txt"), "w") do io
            println(io, "# Raw gate sequence of size $(sum(gate_counts))")
            for (i, g) in enumerate(baseline_gates)
                println(io, i, "\t", repr(g))
            end
        end

        # open(joinpath(dir, "baseline_exec_circuit.txt"), "w") do io
        #     println(io, "# Executable (DQC) circuit operations of sizes $(circuit_sizes) (excl. SWAPS)")
        #     for (i, op) in enumerate(baseline_exec_circuit)
        #         println(io, i, "\t", repr(op))
        #     end
        # end

        open(joinpath(dir, "summary.txt"), "w") do io
            println(io, "# Encoding successful: $verification_logical_state")
            println(io, "# Raw gate sequence of size $(sum(gate_counts))")
            println(io, "# Executable (DQC) circuit operations of size $(gate_counts) (excl. SWAPS)")
        end
    end
        ###########################################
        ###########################################
        ###########################################

    return baseline_gates
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
        _, circuit, _ = standard_logical_zero_encoding_circuit(Steane7())
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
    encoded_plus = S"XIXIXIX IXXIIXX IIIXXXX ZIZZIZI  ZZIIZZI ZZIZIIZ IIXIXXI" # works with XXXXXXX as logical X as well

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