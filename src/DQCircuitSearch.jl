module DQCircuitSearch

include("types.jl")
include("helper.jl")
include("trivariate_bicycle_code.jl")
include("logical_encoding.jl")
include("plots.jl")
include("dtsimulation.jl")
include("circsim.jl")
include("genetic.jl")
include("parameters.jl")

using .Genetic: run_genetic_search
export run_genetic_search

using .Parameters: run_parameter_sweep
export run_parameter_sweep

using .Circuit_Plots: plot_gate_teleportation
export plot_gate_teleportation

using .LogicalEnc: run_tests
export run_tests

end