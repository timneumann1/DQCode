# types.jl

module Types

export SimulationParameters, SimulationFidelity

struct SimulationParameters
    register_sizes::Vector{Int}        # Number of qubits in each register
    T1_relaxation::Float64          # T1 relaxation time of all qubits   [Wang]
    T2_dephasing::Float64         # T2 dephasing time of all qubits     [Wang]
    fidelity_bell_pairs::Float64             # Fidelity of the raw Bell pairs
    simulation_time::Float64
    single_qubit_gate_exec_time::Float64   # Execution Time of a single qubit gate
    two_qubit_gate_exec_time::Float64      # Execution Time of a two qubit gate
    projective_measurement_time::Float64      # Time to peform a measurement
    init_time::Float64                         # Time to initialise the system (e.g., in the all-zero state)
    classical_communication_time::Float64   # define in relation to speed of photons in fiber cable
    bell_state_fidelity::Float64
    success_prob::Float64
    attempt_time::Float64
end

struct SimulationFidelity
    fidelity::Float64
end

end