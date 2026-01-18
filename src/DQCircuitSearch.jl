module DQCircuitSearch

include("types.jl")
include("simulation.jl")
include("genetic.jl")

using .Genetic: run_genetic_search
#export SimulationParameters, SimulationFidelity, Circuit, HadamardGate, IdentityGate, PauliXGate, PauliYGate, PauliZGate, CNOT_Gate, Gate
export run_genetic_search

end