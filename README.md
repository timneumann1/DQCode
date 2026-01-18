A package for Quantum Architecture Search (QAS) on Distributed Quantum Systems. 
The module DQCircuitSearch orchestrates the execution of the search by initialising a genetic search via 
```
julia
activate .
using Revise
include("src/DQCircuitSearch.jl"); using .DQCircuitSearch; DQCircuitSearch.run_genetic_search()
```