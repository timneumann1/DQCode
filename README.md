# DQCode: Simulating DQC-optimised fault-tolerant logical $|0\rangle_L$ state encoding

This repository investigates fault-tolerant logical zero state encoding of quantum CSS codes for distributed Type-II quantum architectures. It contributes two main functionalities: (i) a circuit search tool, with which the user can optimise unitary encoding circuits for arbitrary CSS stabiliser codes, and (ii) a DQC simulation tool, with which the respective encoding circuits can be simulated under realistic noise and networking conditions in a distributed setting.

All source files in `src/` perform one step of the DQCode pipeline, orchestrated by `DQCode.jl`, which serves as the central entrypoint to all functionalities.

## Setup

In order to use this code, a number of setup steps should be followed. Firstly, pick a disk location and create a project folder, e.g., `ft_zero_encoding/`.

To access the verification circuit synthesis provided by MQT-QECC, one needs to install the author's MQT-QECC fork which DQCode interacts with. For this, enter your project folder 

```
cd {path/to/your/project/folder}
```

Next, clone the Python repository via

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

Now, we install our Julia repository DQCode. First, navigate back to your project root folder (in which you should now find the `qecc/` folder) using

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

At this point, you should have two separate repositories in your project folder, the DQCode Julia repository, and a local copy of the MQT-QECC fork.

Since we want those codebases to work together, we use the Julia library `PyCall`. 

Still in the Julia REPL, we create the desired link by running

```
ENV["PYTHON"] = joinpath({QECC_PATH}, ".venv/bin/python3") # replace with your MQT_PATH set above
using Pkg
Pkg.build("PyCall")
exit()
```

In the above, `QECC_PATH` is the path to your local MQT copy, i.e. 

```
{path/to/your/project/folder/qecc}
```

Now PyCall acts as the bridge, making the Python repository readable and usable from Julia.

Congrats, you've made it through; this completes the setup. (If you choose to alter the Python repository, e.g., by pulling new additions to the original library, these changes will be reflected from Julia directly). 

## Usage

Available code-network configurations are defined (and can easily be appended) in the [experiment configuration](src/experiment_config.jl) file. 

To initialise a specific code-architecture configuration, indicate the corresponding configuration string in the [setup script](scripts/execution/dqc_setup.jl) and run

```
julia
] activate .
include("scripts/execution/dqc_setup.jl")
```
which stores the `NetworkingSpecifications` and `CodeParameters` to the `data/`folder.

Using `DQCode.jl` as entrypoint to the DQCode functionalities, the `scripts/execution` folder contains the scripts needed to execute the entire DQCode pipeline, making optimisations and simulations accessible from the REPL. 

- For the optimisation/circuit search, we currently expose a Monte Carlo Tree Search [(MCTS)](src/mcts.jl) and warm-start [Genetic Algorithm](src/genetic.jl). These can be accessed via 

```
include("scripts/execution/dqc_mcts.jl")
```

and

```
include("scripts/execution/dqc_gott_ga.jl")
```



- For the [DQC simulation](src/dqc_simulator.jl), individual noise sweeps can be specified in `scripts/execution/dqc_sim.jl` and executed via

```
include("scripts/execution/dqc_sim.jl")
```

Plotting capabilities are exposed in the `scripts/analysis` folder, allowing insight into [optimiser statistics](scripts/analysis/optimiser_evolution.jl) and [logical rate analysis](scripts/analysis/logical_rate.jl).

For an overview over the entire pipeline, with sample outputs based on the Steane [[7,1,3]] code on two QPUs, please refer to the [examples](examples/example.ipynb) notebook. This notebook provides an overview over the entire pipeline and contextualises the role of each of the scripts in DQCode. 


## Reproducibility Information

```
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

Functionality has been tested on MacOS.

We chose to collect all hyperparameter data in the [configuration file](src/experiment_config.jl) for reproducibility. Upon executing the aforementioned DQCode pipeline, you can then collect all data into the `data/` folder. For completeness, we publish all data to the folder `project_data/`. 


## Acknowledgements

DQCode uses ...
- ... the QuantumClifford library (https://github.com/QuantumSavory/QuantumClifford.jl) @cite[https://arxiv.org/abs/2512.16752]
- ... the Munich Quantum Toolkit QECC library (https://github.com/munich-quantum-toolkit/qecc/tree/03a62ca1d3ccbe690265d6b5a7c59c5f72681793) @cite[https://arxiv.org/abs/2408.11894]

- ... the Karlsruhe Hypergraph Partitioning (KaHyPar) Algorithm (https://kahypar.org, https://github.com/kahypar/kahypar)
- ... Qiskit (https://github.com/Qiskit/qiskit/blob/stable/2.4/qiskit/synthesis/clifford/clifford_decompose_bm.py#L25-L48 )
- ... QuantikZ @cite[arXiv:1809.03842].

## Citation
```
{author}: Tim Neumann
{title}: Simulating DQC-optimised fault-tolerant logical $|0\rangle_L$ state encoding
{year}: 2026
```