module DQCircuitSearch

include("types.jl")
include("simulation.jl")
include("genetic.jl")

using .Genetic: run_genetic_search
export run_genetic_search

end