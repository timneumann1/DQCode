module Genetic

using ..Types
using ..CircuitSimulator
#using ..Helper
using ..LogicalEnc

using Quantikz: savecircuit
using QECCore: Steane7
using QuantumClifford: MixedDestabilizer, sHadamard, sCNOT, sSWAP, @S_str, true_success_stat, false_success_stat, continue_stat, failure_stat
using BenchmarkTools


export run_genetic_search


function define_parameters()
    params = SimulationParameters(
        [3,4], #register sizes
        #12000.0,#T1
        10, #4200 depolarising noise time
        20e-6, # Execution Time of a single qubit gate   #20e^-6
        200e-6,   # Two-qubit gates 200e^-6
        1e-5,  #1e^-5  projective measurement time
        1e-2,  # 1e^-2 classical comm time
        1,#0.9689, # Bell state fidelity  #TODO: make this smaller than 1, should be fixed upon MixedDestabilizer fidelity
        1.41e-4, # Bell state generation,from [Main, 2025]    success probability 
        1.168e-9 #1.168e-9  # attempt time

    #TODO: Define units and insert realistic values
    #TODO: Add state preparation fidelity and single-shot readout of 99.93% [Harty], single-qubit gate fidelity of 99.99916%, two-qubit fidelity of 99.97% [Löschnauer]
    # characteristic_time = 1000
    # p = 1-exp(-1/characteristic_time) # define probability for Pauli Noise application (Poisson point process)
    #TODO: Define single qubit error rate
    )
    return params
end


function build_start_circuit(num_qubits)
    #_, circuit = naive_encoding_circuit(Steane7())


    circuit = Circuit(num_qubits, 12)   # params.register_sizes rows (qubits) and 8 columns (time steps)
    
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


    #print(typeof(circuit))
    #savecircuit(circuit, "src/plots/circuit_sim/circuit.png") # plotting is performed by enabling the reset function

    return circuit.gates
    # could start with circuit from logical_encoding.jl here
end


function run_genetic_search()

    params = define_parameters()                             # retrieve parameters

    # TODO: Mapping stage -> use dictionary to map indices to one another
    permutation = [1,7,4,2,3,5,6]
    inv_perm = invperm(permutation)
    mapping = perm_to_transpositions(deepcopy(permutation)) # careful: this does in-place substitution of permutation
    # as extracted from Hypergraph Partitoning DO I WANT TO DO THIS HERE ONCE AND ALWAYS JSUT PASS IT?

    register_lookup_array, data_qubits, num_data_qubits = create_lookup_array_cliff(params.register_sizes)      # create lookup array
    num_comm_qubits_per_register = length(params.register_sizes)-1
    num_qubits = num_data_qubits + num_comm_qubits_per_register*(length(params.register_sizes)) +1 # one verification qubit
    print("number of qubits is $num_qubits, $num_comm_qubits_per_register")
    #mapping = [(7,2),(6,2),(5,2),(4,3),(3,2)]  #this mapping is an update of the oroginal transpoitions, taking into account that we inserted comm qubits
    # Make array of data qubits
    println("Lookup Array: $register_lookup_array")
    println("Data qubits: $data_qubits")
    
    circuit_tensor = build_start_circuit(num_qubits)                  # build initial circuit
    #target_state = S"XIXIXIX IXXIIXX IIIXXXX ZIZZIZI ZZIIZZI ZZIZIIZ IZIZIZI"
    # convert tensor of DATA QUBITS to QS circuit
    circuit = tensor_to_circuit(circuit_tensor, mapping, inv_perm, register_lookup_array, data_qubits, num_comm_qubits_per_register, num_qubits)
    
    #@btime tensor_to_circuit($circuit_tensor, $mapping, $inv_perm, $register_lookup_array, $data_qubits, $num_comm_qubits_per_register, $target_state)

    
    
    #circuit = add_verification(circuit, target_state, data_qubits)
    
    #savecircuit(circuit, "src/plots/circuit_sim/circuit_noise.png") # plotting is performed by enabling the reset function



    #TODO: block all communication qubit layers! Can be done via row check != comm_qubits,
    #TODO: Include check for no overlaps within one layer

    # Pauli measuremnt with project!
    num_traj = 50000
    mc_result = execute_circuit(circuit, num_qubits; num_traj = num_traj) # if specificg num_traj = 100000, we use mc sampling, otherwise pert.
    # for perturbative expansion, only the leading order is kept, so probabilies can be smaller than 1
    println()
    @btime execute_circuit($circuit, $num_qubits, num_traj = $num_traj)


    println("\nFinal Steane-7 dict: $(mc_result) \n")
    if (mc_result[true_success_stat]  + mc_result[false_success_stat]) != num_traj
        throw(ErrorException("Some runs were invalid"))
    end

    fidelity = (round(mc_result[true_success_stat] / (mc_result[true_success_stat]+mc_result[false_success_stat]),digits=3))
    println("This is a fidelity of $fidelity")
end

end
