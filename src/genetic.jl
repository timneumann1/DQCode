module Genetic

using ..Types
using ..CircuitSimulator
using ..Helper
using ..LogicalEnc

using Quantikz: savecircuit, @with, classicalbitslayout
using QECCore: Steane7
using QuantumClifford: MixedDestabilizer, sHadamard, sCNOT, sSWAP, @S_str, true_success_stat, false_success_stat, continue_stat, failure_stat, PauliMeasurement, VerifyOp
using BenchmarkTools

export run_genetic_search


function define_parameters()
    params = GeneticParameters(
        [3,4], #register sizes
        0.01, # depolarising_prob 
        0.01, # gate_noise_prob 
    )
    return params
end


function build_start_circuit(num_data_qubits)

    # could start with circuit from logical_encoding.jl here
    # (need to be recast from Vector{AbstractOperation} to tensor)

    circuit = Circuit(num_data_qubits, 12) 
    
    circuit.gates[1,1] = HadamardGate()
    circuit.gates[2,1] = HadamardGate()
    circuit.gates[3,1] = HadamardGate()

    circuit.gates[7,2] = circuit.gates[4,2] = CNOT_Gate(7,4)

    circuit.gates[1,3] = circuit.gates[4,3] = CNOT_Gate(1,4)
    circuit.gates[7,3] = circuit.gates[5,3] = CNOT_Gate(7,5)

    circuit.gates[1,4] = circuit.gates[5,4] = CNOT_Gate(1,5)

    circuit.gates[1,5] = circuit.gates[6,5] = CNOT_Gate(1,6)

    circuit.gates[2,6] = circuit.gates[4,6] = CNOT_Gate(2,4)

    circuit.gates[2,7] = circuit.gates[6,7] = CNOT_Gate(2,6)

    circuit.gates[2,8] = circuit.gates[7,8] = CNOT_Gate(2,7)

    circuit.gates[3,9] = circuit.gates[5,9] = CNOT_Gate(3,5)

    circuit.gates[3,10] = circuit.gates[6,10] = CNOT_Gate(3,6)

    circuit.gates[3,11] = circuit.gates[7,11] = CNOT_Gate(3,7)

    circuit.gates[3,12] = circuit.gates[4,12] = SWAP_Gate(3,4)
    circuit.gates[6,12] = circuit.gates[7,12] = SWAP_Gate(6,7)

    #savecircuit(circuit, "src/plots/circuit_sim/circuit.png") # plotting is performed by enabling the reset function

    return circuit.gates
end


function run_genetic_search()

    params = define_parameters()                             # retrieve parameters

    # TODO: Mapping stage -> use dictionary to map indices to one another
    # As extracted from Hypergraph Partitoning
    permutation = [1,7,4,2,3,5,6]
    inv_perm = invperm(permutation)
    mapping = perm_to_transpositions(deepcopy(permutation)) # careful: without deepcopy, this does in-place substitution of permutation
    
    # NOTE: When generating the infromation for hypergraph part., we need to consult the naive encoding function in the logical encoding script to obtain the logical oeprators.
    # For the inversion of the circuit, we have a custoim function in circsim.jl since this requries applicaiton of correct indices, accounting for communication qubits.

    data_qubit_capacities = params.register_sizes
    num_registers  = length(data_qubit_capacities)
    register_lookup_array, data_qubits, num_data_qubits = create_lookup_array_cliff(data_qubit_capacities)      # create lookup array
    num_comm_qubits_per_register = num_registers-1
    num_qubits = num_data_qubits + num_comm_qubits_per_register*(num_registers) # one verification qubit
    println("number of qubits is $num_qubits, number of comm. qubits per register is $num_comm_qubits_per_register")
    println("Lookup Array: $register_lookup_array")
    println("Data qubits: $data_qubits")
    
    circuit_tensor = build_start_circuit(num_data_qubits)                  # build initial circuit
    target_state = S"XIXIXIX IXXIIXX IIIXXXX ZIZZIZI ZZIIZZI ZZIZIIZ IZIZIZI"
    code = Steane7()

    circuit = tensor_to_circuit(code, params.depolarising_noise, params.gate_noise, circuit_tensor, mapping, inv_perm, register_lookup_array, data_qubits, num_comm_qubits_per_register, num_qubits, target_state, data_qubit_capacities)
    #@btime tensor_to_circuit($code, $params.depolarising_noise, $params.gate_noise, $circuit_tensor, $mapping, $inv_perm, $register_lookup_array, $data_qubits, $num_comm_qubits_per_register, $num_qubits, $target_state, $data_qubit_capacities)

     # Two methods of verifying the creation of the encoded state (Method 1 is preferable since simpler)

    # 1. VerifyOp 
    push!(circuit, VerifyOp(target_state, data_qubits)) 

    """
    # 2. Apply inverse unitary and measure physical all-zero state
    #    Assumes the existence of an extra qubit for verification
    
    encoding_circuit = naive_encoding_circuit_mapping(code, num_comm_qubits_per_register, register_lookup_array)

    for gate in reverse(encoding_circuit)
        push!(circuit, gate)
    end
   
   
    # push!(circuit, VerifyOp(S"ZIIIIII IZIIIII IIZIIII IIIZIII IIIIZII IIIIIZI IIIIIIZ", data_qubits)) 
   
    # or 

    circuit, pauli_string = measure_zero(circuit, data_qubits, num_qubits) # one verification qubit is appended

    #  The measurment of zero yields true in case the eigenvalue of the state is -1. That is, when the zero state actually lives in the register, 
    #  the measurement outcome will be +1, so the boolean will be False. In that case (when we measured zero), we apply the flip 
    push!(circuit,ConditionalGate(sId1(num_qubits),sX(num_qubits), pauli_string.bit))
    println("circuit after conditional: $circuit")
    # If we have measured zero and thus applied the flip, the qubit at index num_qubits will be in the |1> state, which is stabilise by -Z
    push!(circuit, VerifyOp(S"-Z", [num_qubits])) # 
    """
        
    max_classical_bit = maximum([op.bit for op in circuit if op isa PauliMeasurement] ; init=0)
    println("No. of classical bits: $max_classical_bit")

    # @with classicalbitslayout => :expanded begin
    #     savecircuit(circuit, "src/plots/circuit_sim/circuit_noise.png")
    # end

    #savecircuit(circuit, "src/plots/circuit_sim/circuit_noise.png", nclassical = max_classical_bit+1) # plotting is performed by enabling the reset function

    num_traj = 10000
    mc_result = execute_circuit(circuit, num_qubits, num_registers; num_traj = num_traj) # if specifying num_traj, we use MC sampling, otherwise perturbation.
    # for perturbative expansion, only the leading order is kept, so probabilities can be smaller than 1, 
    # also, PauliMeasurement don't work with pert. expansion currently
    println()
    @btime execute_circuit($circuit, $num_qubits, $num_registers, num_traj = $num_traj)

    println("\nFinal Steane-7 dict: $(mc_result) \n")
    if (mc_result[true_success_stat]  + mc_result[false_success_stat]) != num_traj
        throw(ErrorException("Some runs were invalid"))
    end

    fidelity = (round(mc_result[true_success_stat] / (mc_result[true_success_stat]+mc_result[false_success_stat]),digits=3))
    println("This is a fidelity of $fidelity")
end

end
