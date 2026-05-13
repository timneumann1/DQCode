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

Tims-MacBook-Air:thesis tim$ julia --project=envs/opt -t 2 -e 'include("src/DQCircuitSearch.jl"); using .DQCircuitSearch; DQCircuitSearch.run_genetic_search("trivariate_3_3_3_3")' 2>&1 | tee src/results/output.log

Installation procedure:

install working version of MQT QECC library based on the MQT installation procedure. Make all necessary local changes in this fork. We will use the virtual environment of this library as our PyCall Python environment. Then run ENV["PYTHON"] = "/Users/tim/Tim/projects/mqt/qecc/.venv/bin/python3" Pkg.build("PyCall") from julia repo only after that using PyCall as in the files

the MQT installation already contains qiskit so this comes out of box

To run MQT installation, follow the intallation procedure: install uv via the strange command, fork my repo/their repo, make sure Python is installed, then uv sync
uv sync --python 3.13
for visuals: uv pip install pylatexenc

Need to build qiskit with PyCall -> reference to installation