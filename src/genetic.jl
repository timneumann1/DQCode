include("types.jl")
include("simulation.jl")

using .Types

params = Types.SimulationParameters(
        [6,3],
        12000.0,
        4200.0,
        1.0,
        10.0,  
        20e-6, # Execution Time of a single qubit gate
        200e-6,
        1e-5,
        1e-2,
        0.001,
        0.9689, # Bell state fidelity
        1.41e-4, # Bell state generation,from [Main, 2025]
        1.168e-9

#TODO: Define units and insert realistic values
#TODO: Add state preparation fidelity and single-shot readout of 99.93% [Harty], single-qubit gate fidelity of 99.99916%, two-qubit fidelity of 99.97% [Löschnauer]
# characteristic_time = 1000
# p = 1-exp(-1/characteristic_time) # define probability for Pauli Noise application (Poisson point process)
#TODO: Define single qubit error rate
)

# TODO: define the circuit here

fidelity = run_simulation(params)

print("\nFinal Steane-7 fidelity: $(fidelity.fidelity) \n")


