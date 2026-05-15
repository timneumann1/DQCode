This repository investgates fault-tolerant logical state encoding for distributed Type-II quantum architectures.

It contributes two main functionalities: (i) a circuit search tool, with which the user can optimise unitary encoding circuits for arbitrary CSS stabiliser codes, and (ii) a 
DQC simulation tool, with which the respective encoding circuits can be simulated under realistic noise and networking conditions in a distributed setting.

For the optimisation/circuit search, currently a Monte Carlo Tree Search [(MCTS)](src/mcts.jl) and warm-start [Genetic Algorithm](src/genetic.jl) are available. For the [DQC simulation](src/dqc_simulator.jl),
we take a bare encoding circuit from the optimisation step, append a (flag-based) verification circuit from the Munich Quantum Toolkit QECC package, and then determine the logical error rate of fault-tolerant encoding after inserting the required
telegate operations between cores, noise channels and post-processing.

Available codes can be inspected (and easily appended) in the [experiment configuration](src/experiment_config.jl) file.

All optimisations and simulations are exposed via [scripts](scripts).

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

Still in the Julia REPL, we now create the required link by running

```
ENV["PYTHON"] = joinpath({MQT_PATH}, ".venv/bin/python3") # replace with your MQT_PATH set above
using Pkg
Pkg.build("PyCall")
exit()
```

This completes the setup. Now we can run scripts from the Julia REPL, for example

```
julia
] activate .
include("scripts/execution/dqc_mcts.jl")
```

In the scripts files, one can indicate the configuration (code and architecture) of the setup one wishes to investigate. Don't forget to run the [setup](scripts/execution/dqc_setup.jl) script once before
to store the required configurations. 


We use
- QuantumClifford
- KaHyPar
- Qiskit, and BM alg for baseline (https://github.com/Qiskit/qiskit/blob/stable/2.4/qiskit/synthesis/clifford/clifford_decompose_bm.py#L25-L48 )
- MQT QECC
- QuantumClifford (esp. encoding circuit)