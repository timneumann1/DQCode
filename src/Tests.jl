
module UnitTests

include("types.jl")
include("trivariate_bicycle_code.jl")
include("helper.jl")
include("experiment_config.jl")
#include("circsim.jl")
include("encoding_gott.jl")

using .Helper: perm_to_transpositions, transpositions_to_perm, compare_states
using .EncodingGott: encoding_circuit_gott
using QuantumClifford: Register, MixedDestabilizer, mctrajectory!, stabilizerview, logicalzview, VerifyOp, true_success_stat, sX

using Serialization
# Permutations

function permutation_test()

    # ---------- TEST 1 ----------

    # Our convention is: mapping contains, at index i, the qubit j that will be mapped to this position in the DQC setting, "position i gets qubit j"
    mapping = [4,1,2,5,3] # qubit 4 is mapped to 1, qubit 1 is mapped to 2, qubit 2 is mapped to 3, qubit 5 is mapped to 4, qubit 3 is mapped to 5
    # The corresponding product of transpositions should implement the above mapping, when traversed from right to left.
    @assert perm_to_transpositions(mapping) == [(1,2),(2,3),(3,5),(4,5)] "Permutation Test 1 failed" # This transposition implements "qubit i > index j": 1>2, 2>3, 3>5, 4>1 and 5>4, as desired 
    
    # ---------- TEST 2 ----------

    transpositions = [(1,2), (2,3), (3,4)] # This transposition implements "qubit i > index j": 1>2, 2>3, 3>4, 4>1 
    @assert transpositions_to_perm(transpositions, 4) == [4,1,2,3] "Permutation Test 2 failed" # The corresponding permutation is [4,1,2,3]

    # ---------- TEST 3 ----------
    # transpositions_to_perm and perm_to_transpositions should yield the same vector when concatenated
    mapping = [5,2,1,3,6,9,4,7,8]
    @assert transpositions_to_perm(perm_to_transpositions(copy(mapping)),9) == mapping "Permutation Test 3 Failed"
    
    @info "Permutation-Transposition Tests: PASSED"
end

function gottesman_encoding()

    # ---------- TEST 1 ----------

    # --- Encoding logical zero ---
    network_specs=deserialize(joinpath("/Users/tim/Tim/projects/thesis/src/results/Steane/[4, 3]/network_specs.jls"))
    code_params=deserialize(joinpath("/Users/tim/Tim/projects/thesis/src/results/Steane/[4, 3]/code_params.jls"))
    basis_state = [0] # Encoding the logical |0> state

    # ---- Executing the encoding circuit ----
    encoding_circuit, _ = encoding_circuit_gott(code_params.qec_code, network_specs, basis_state)
    initial_state = Register(one(MixedDestabilizer, 7))
    encoded_state, stat = mctrajectory!(initial_state, encoding_circuit)
   
    # ---- Constructing Target State ----
    code = MixedDestabilizer(code_params.qec_code)
    target_state = vcat(stabilizerview(code), logicalzview(code))

    # ---- Verifying Equivalence ----
    state, stat = mctrajectory!(encoded_state, [VerifyOp(target_state,collect(1:code_params.n))])
    @assert stat == true_success_stat "Encoding Circuit Test 1 Failed"

    # ---------- TEST 2 ----------

    # --- Encoding logical one ---
    network_specs=deserialize(joinpath("/Users/tim/Tim/projects/thesis/src/results/Steane/[4, 3]/network_specs.jls"))
    code_params=deserialize(joinpath("/Users/tim/Tim/projects/thesis/src/results/Steane/[4, 3]/code_params.jls"))
    basis_state = [1] # encoding the logical |1> state

    # ---- Executing the encoding circuit ----
    encoding_circuit, _ = encoding_circuit_gott(code_params.qec_code, network_specs, basis_state)
    initial_state = Register(one(MixedDestabilizer, code_params.n))
    encoded_state, stat = mctrajectory!(initial_state, encoding_circuit)
   
    # ---- Constructing Target State ----
    code = MixedDestabilizer(code_params.qec_code)
    inital_state = vcat(stabilizerview(code), logicalzview(code))
    # Logical X operator: + __X_XX_ ... = X_3 X_5 X_6
    target_state, _ = mctrajectory!(inital_state, [sX(3), sX(5), sX(6)] ) # applying logical X_L operator to the logical zero state
    
    # ---- Verifying Equivalence----
    state, stat = mctrajectory!(encoded_state, [VerifyOp(target_state, collect(1:code_params.n))])
    @assert stat == true_success_stat "Encoding Circuit Test 2 Failed"

    # ---------- TEST 3 ----------

    # --- Encoding logical |0010> ---
    network_specs=deserialize(joinpath("/Users/tim/Tim/projects/thesis/src/results/BivariateBicycle/[3, 3, 3, 3, 3, 3]/network_specs.jls"))
    code_params=deserialize(joinpath("/Users/tim/Tim/projects/thesis/src/results/BivariateBicycle/[3, 3, 3, 3, 3, 3]/code_params.jls"))
    basis_state = [0,0,1,0] # encoding the logical |0010> state

    # ---- Executing the encoding circuit ----
    encoding_circuit, _ = encoding_circuit_gott(code_params.qec_code, network_specs, basis_state)
    initial_state = Register(one(MixedDestabilizer, code_params.n))
    encoded_state, stat = mctrajectory!(initial_state, encoding_circuit)
   
    # ---- Constructing Target State ----
    code = MixedDestabilizer(code_params.qec_code)
    inital_state = vcat(stabilizerview(code), logicalzview(code))
    # Logical X operators: + _______X_XXXXX____
                         # + _______XX___X_X___
                         # + _______XX_XX___XX_   ... X_{L_3} = X_8 X_9 X_11 X_12 X_16 X_17 
                         # + ________X__X___X_X
    target_state, _ = mctrajectory!(inital_state, [sX(8), sX(9), sX(11), sX(12), sX(16), sX(17)] ) # applying logical X_{L_3} operator to the logical zero state
    
    # ---- Verifying Equivalence----
    state, stat = mctrajectory!(encoded_state, [VerifyOp(target_state, collect(1:code_params.n))])
    @assert stat == true_success_stat "Encoding Circuit Test 3 Failed"

    # ---------- TEST 4 ----------

    # --- Encoding logical |11> ---
    network_specs=deserialize(joinpath("/Users/tim/Tim/projects/thesis/src/results/TrivariateBicycle/[3, 3, 3, 3]/network_specs.jls"))
    code_params=deserialize(joinpath("/Users/tim/Tim/projects/thesis/src/results/TrivariateBicycle/[3, 3, 3, 3]/code_params.jls"))
    basis_state = [1,1] # encoding the logical |11> state

    # ---- Executing the encoding circuit ----
    encoding_circuit, _ = encoding_circuit_gott(code_params.qec_code, network_specs, basis_state)
    initial_state = Register(one(MixedDestabilizer, code_params.n))
    encoded_state, stat = mctrajectory!(initial_state, encoding_circuit)
   
    # ---- Constructing Target State ----
    code = MixedDestabilizer(code_params.qec_code)
    inital_state = vcat(stabilizerview(code), logicalzview(code))
    # Logical X operators: + _____X__X_X_ ... X_{L_1} = X_6 X_9 X_11
                        #  + _____XXX_X_X ... X_{L_2} = X_6 X_7 X_8 X_10 X_12
    target_state, _ = mctrajectory!(inital_state, [sX(7), sX(8), sX(9), sX(10), sX(11), sX(12)] ) # applying logical X_{L_1}X_{L_2} operator to the logical zero state
    
    # ---- Verifying Equivalence----
    state, stat = mctrajectory!(encoded_state, [VerifyOp(target_state, collect(1:code_params.n))])
    @assert stat == true_success_stat "Encoding Circuit Test 4 Failed"
    
    @info "Encoding Circuit Tests: PASSED"
end

end