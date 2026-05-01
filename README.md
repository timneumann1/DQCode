A package for Quantum Architecture Search (QAS) on Distributed Quantum Systems. 
The module DQCircuitSearch orchestrates the execution of the search by initialising a genetic search via 
```
julia
activate .
using Revise
include("src/DQCircuitSearch.jl"); using .DQCircuitSearch; DQCircuitSearch.run_genetic_search()
```

The module DQCircuitSearch exposes the code parameter and network specification setup, and useful functions to orchestrate circuit searches manually or at scale (with the help of the ExperimentConfig module).

DQCodePrep exposes the FT synthesis and execution

Tests provides test of functionality

Installation procedure:

install working version of MQT QECC library based on the MQT installation procedure. Make all necessary local changes in this fork. We will use the virtual environment of this library as our PyCall Python environment. Then run ENV["PYTHON"] = "/Users/tim/Tim/projects/mqt/qecc/.venv/bin/python3" Pkg.build("PyCall") from julia repo only after that using PyCall

the MQT installation already contains qiskit so this comes out of box