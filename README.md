# DQCode: Optimising fault-tolerant zero-state encoding on Type-II distributed quantum architectures

This repository investgates fault-tolerant logical zero state encoding for distributed Type-II quantum architectures.

It contributes two main functionalities: (i) a circuit search tool, with which the user can optimise unitary encoding circuits for arbitrary CSS stabiliser codes, and (ii) a 
DQC simulation tool, with which the respective encoding circuits can be simulated under realistic noise and networking conditions in a distributed setting.

For the optimisation/circuit search, we currently expose a Monte Carlo Tree Search [(MCTS)](src/mcts.jl) and warm-start [Genetic Algorithm](src/genetic.jl). For the [DQC simulation](src/dqc_simulator.jl),
we take some (potentially optimised) encoding circuit for the logical zero state, append a (flag-based) verification circuit from the Munich Quantum Toolkit QECC package, and then determine the logical error rate of the fault-tolerant encoding under circuit-level and Bell pair initialisation noise, mimicing realistic noise channels under telegate execution.

Available code-network configurations are defined (and can easily be appended) in the [experiment configuration](src/experiment_config.jl) file.

All optimisations and simulations are accessible via [scripts](scripts). To test differnet parameters for the optimisations, one can also change the respective entries in the [experiment configuration](src/experiment_config.jl) file. We chose to collect all hyperparameter data in this file for better reproducibility. 

For the DQC evaluation, individual noise sweeps can be specified in the [DQC simulation script](scripts/execution/dqc_sim.jl). Plitting capabilities are exposed in the respective [logical rate analysis script](scripts/analysis/logical_rate.jl).

In order to use this code, a number of setup steps should be followed.

In order to access the MQT-QECC verification circuit, one also needs to install the MQT-QECC fork tailored to work hand in hand with our Julia repository DQCode. For this, choose a project folder on your disk and 
enter it via

```
cd path/to/your/project/folder
```

Next, clone the Python work via

```
git clone https://github.com/timneumann1/qecc.git
```

To harness the capabilities of this software toolkit, we install the required dependencies in a virtual environment

```
cd qecc
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -U pip
python3 -m pip install -e .
```

Now we install our Julia repository DQCode. First, navigate back to your project folder (in which you should now find the qecc folder) using

```
cd ..
```

Then run 

```
git clone https://github.com/timneumann1/DQCode.git
```

and enter the repository via 

```
cd DQCode
```

To install the required Julia dependencies, enter the Julia REPL, activate the default project location and install the dependencies via
```
julia 
] activate .
] instantiate
```

At this point, you should have two separate local repositories in your project folder, the DQCode Julia repository, and a local copy of the MQT-QECC fork qecc.

In order to make them work together, we use the Julia library PyCall. For this, enter the path to your qecc library in DQCode, replacing the current path

```
const MQT_PATH = "/Users/tim/Tim/projects/mqt/qecc/"
```
no, just install it in the same folder

Still in the Julia REPL, we now create the required link by running

```
ENV["PYTHON"] = joinpath({MQT_PATH}, ".venv/bin/python3") # replace with your MQT_PATH set above
using Pkg
Pkg.build("PyCall")
exit()
```

PyCall acts as the brdige, making the Python repositroy readable and usable from Julia.

This completes the setup (the Python repository can also be adepted at this point, for example updated to new additions to the original library; however, we cannot guarantee for this to happen without conflicts). Now we can run scripts from the Julia REPL, for example

```
julia
] activate .
include("scripts/execution/dqc_mcts.jl")
```

In the scripts files, one can indicate the configuration (code and architecture) of the setup one wishes to investigate. Don't forget to run the [setup](scripts/execution/dqc_setup.jl) script once before
to store the required configurations. 

The examples file provides an overview over the entire pipeline to contextualise the role of each of the scripts you can call via the above procedure. While this contextualises the available functions, we recommend to use the above code snippet (replacing dqc_mcts.jl with the desired function) to run actual optimisation and evaluation experiments.

Functionality has only been tested on MacOS.

*Reproducibility Information:*

````
julia> versioninfo()
Julia Version 1.12.6
Commit 15346901f00 (2026-04-09 19:20 UTC)
Build Info:
  Official https://julialang.org release
Platform Info:
  OS: macOS (arm64-apple-darwin24.0.0)
  CPU: 8 × Apple M3
  WORD_SIZE: 64
  LLVM: libLLVM-18.1.7 (ORCJIT, apple-m3)
  GC: Built with stock GC
Threads: 8 default, 1 interactive, 4 GC (on 4 virtual cores)
```


We use
- QuantumClifford
- KaHyPar
- Qiskit, and BM alg for baseline (https://github.com/Qiskit/qiskit/blob/stable/2.4/qiskit/synthesis/clifford/clifford_decompose_bm.py#L25-L48 )
- MQT QECC
- QuantumClifford (esp. encoding circuit)
- QuantikZ arXiv:1809.03842.

Citation
author: Tim Neumann
title: DQCode -- ...
