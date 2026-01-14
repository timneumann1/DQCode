# parameters.jl

include("types.jl")
include("simulation.jl")

using .Types

#TODO: sweep over parameters and get fidelities

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
)

#for ...
fidelity = run_simulation(params)

print("\nFinal Steane-7 fidelity: $(fidelity.fidelity) \n")

# also need to pass circuit obejct