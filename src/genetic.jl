module Genetic

using ..Types
using ..Simulation

export run_genetic_search

function define_parameters()
    params = SimulationParameters(
        [5,4], #register sizes
        12000.0,#T1
        100000000.0, #4200 T2
        20e-6, # Execution Time of a single qubit gate   #20e^-6
        200e-6,   # Two-qubit gates 200e^-6
        1e-5,  #1e^-5  projective measurement time
        1e-2,  # 1e^-2 classical comm time
        1,#0.9689, # Bell state fidelity
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

# Need to be run only once before the EA
#TODO: Take care of scoping of register and register_start_index
function create_lookup_array(params)
    register_lookup_array = Vector{Int}(undef, sum(params.register_sizes))
    register = 1
    register_start_index = 1
    for i in eachindex(params.register_sizes)
        j = params.register_sizes[i]
        register_lookup_array[register_start_index:register_start_index+j-1] .= register
        register_start_index+=j
        register +=1
    end
    return register_lookup_array
end

function comm_qubits_array(params)
    comm_qubits_array = Vector{Int}(undef,length(params.register_sizes))
    index = 1
    for i in eachindex(params.register_sizes)
        j = params.register_sizes[i]
        comm_qubits_array[i] = index
        index += j
    end
    return comm_qubits_array
end

# register_lookup_array = Int[]
# for (register, size) in enumerate(params.register_sizes)
#     append!(register_lookup_array, fill(register, size))
# end

function build_start_circuit(params)
    circuit = Circuit(sum(params.register_sizes), 8)   # params.register_sizes rows (qubits) and 8 columns (time steps)
    circuit.gates[2,1] = HadamardGate()
    circuit.gates[3,1] = HadamardGate()
    circuit.gates[5,1] = HadamardGate()

    circuit.gates[2,2] = circuit.gates[4,2] = CNOT_Gate(2,4)
    circuit.gates[5,2] = circuit.gates[8,2] = CNOT_Gate(5,8)

    circuit.gates[3,3] = circuit.gates[9,3] = CNOT_Gate(3,9)

    circuit.gates[2,4] = circuit.gates[7,4] = CNOT_Gate(2,7)

    circuit.gates[5,5] = circuit.gates[9,5] = CNOT_Gate(5,9)

    circuit.gates[3,6] = circuit.gates[8,6] = CNOT_Gate(3,8)

    circuit.gates[2,7] = circuit.gates[9,7] = CNOT_Gate(2,9)

    circuit.gates[3,8] = circuit.gates[4,8] = CNOT_Gate(3,4)
    circuit.gates[5,8] = circuit.gates[7,8] = CNOT_Gate(5,7)
    return circuit
end


function run_genetic_search()
    params = define_parameters()
    register_lookup_array = create_lookup_array(params)
    # TODO: define the circuit here
    circuit = build_start_circuit(params)
    # TODO: block all communication qubit layers! Can be done via row check != comm_qubits,
    #TODO: Include check for no overlaps within one layer
    
    fidelity = run_simulation(params, circuit, register_lookup_array)

    print("\nFinal Steane-7 fidelity: $(fidelity.fidelity) \n")

end


end
