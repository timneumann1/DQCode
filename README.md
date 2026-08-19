# DQCode: Synthesis of Fault-Tolerant $|0\rangle_L^{\otimes k}$ State Preparation Circuits for Small qLDPC-CSS Codes On Distributed Quantum Architectures

This repository investigates fault-tolerant logical zero state encoding of quantum CSS codes on distributed Type-II quantum architectures. It contributes two main functionalities: 
1. Circuit Search, allowing the user to optimise unitary encoding circuits for arbitrary CSS stabiliser codes, and
2. DQC Simulation, exeucting the respective encoding circuits under realistic noise and networking conditions in a distributed setting.

All files in `src/` perform one step of the $\texttt{DQCode}$ pipeline, orchestrated by `DQCode.jl`, which serves as the central entrypoint to all functionalities.

<p align="center">
<img src="./data/Steane/%5B4,%203%5D/2d_heatmap_ratio_FT.png" alt="2d ratio" width="250"/>
</p>

All data produced in this project can be accessed at [https://github.com/timneumann1/DQCode_Data](https://github.com/timneumann1/DQCode_Data).

## Setup

In order for you to use $\texttt{DQCode}$, a number of setup steps should be followed. Firstly, pick a disk location and create a project folder.

To access the verification circuit synthesis provided by `MQT-QECC`, one needs to install the author's `MQT-QECC` fork with which $\texttt{DQCode}$ interacts. For this, enter your project folder 

```
cd {path/to/your/project/folder}
```

Next, clone the Python repository via

```
git clone https://github.com/timneumann1/qecc.git
```

To harness the capabilities of this software toolkit, we need a working Python version. Please make sure that you have Python $3$ available on your system, or install it from https://www.python.org/downloads/. (This code
has been tested on MacOS for Python $3.13.13$, but should be functional for other Python $3$ releases as well.)
Once Python is available, install the required dependencies in a virtual environment.

```
cd qecc
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -U pip
python3 -m pip install -e .
```

Next, install the Julia repository $\texttt{DQCode}$. For this,, navigate back to your project root folder (in which you should now find the `qecc/` folder), clone the repository and enter it using

```
cd ..
git clone https://github.com/timneumann1/DQCode.git
cd DQCode
```

Before we install the required Julia dependencies, first make sure that you have Julia installed on your system. Instructions for the installation are available via the official Julia language distribution
https://julialang.org/downloads/. This code has been tested on MacOS for `julia 1.12.6` (the latest stable Julia release as of May 2026). Once available, enter the Julia REPL, activate the default project location and install the required dependencies via
```
julia 
] activate .
] instantiate
```

> [!IMPORTANT]
> At this point, you should have two separate repositories in your project folder, the $\texttt{DQCode}$ Julia repository, and a local copy of the MQT-QECC fork.

Since we want those codebases to work together, we use the Julia library `PyCall`. 

Still in the Julia REPL, we create the desired link by running

```
ENV["PYTHON"] = joinpath({QECC_PATH}, ".venv/bin/python3") # replace with your own MQT_PATH (see below)
using Pkg
Pkg.build("PyCall")
exit()
```

In the above, `QECC_PATH` is the path to your local MQT copy, i.e. 

```
{path/to/your/project/folder/qecc}
```

Now PyCall acts as the bridge, making the Python repository readable and usable from Julia.

Congrats, you've made it through; this completes the setup. 

> [!TIP]
> If you choose alter the Python repository, e.g., by pulling new additions to the original library, these changes will be reflected in your virtual environment, and thus be callalble from Julia directly. 

On the $\texttt{DQCode}$ side, you should now see have access to the following file structure:

```
├── examples
│   └── example.md
├── scripts
│   ├── analysis
│   │   ├── logical_rate.jl
│   │   └── optimiser_evolution.jl
│   └── execution
│       ├── dqc_gott_ga.jl
│       ├── dqc_mcts.jl
│       ├── dqc_mqt_baseline.jl
│       ├── dqc_qiskit_baseline.jl
│       ├── dqc_resource.jl
│       ├── dqc_setup.jl
│       ├── dqc_sim.jl
│       └── exec_exp.sh
└── src
    ├── baseline_encoding.jl
    ├── dqc_simulator.jl
    ├── DQCode.jl
    ├── encoding_gott.jl
    ├── experiment_config.jl
    ├── genetic.jl
    ├── helper.jl
    ├── km1_kKaHyPar_sea20.ini
    ├── mcts.jl
    ├── resource_estimate.jl
    ├── symplectic_double_code.jl
    ├── trivariate_bicycle_code.jl
    └── types.jl
```

## Usage

> [!NOTE]
> Here we briefly introduce how to use $\texttt{DQCode}$. For an overview over the entire pipeline, with sample outputs based on the Steane [[7,1,3]] code on two QPUs, please refer to the [examples](examples/example.md) markdown file. This provides an overview over the entire pipeline and contextualises the role of each of the scripts in $\texttt{DQCode}$. 

Available code-network configurations are defined (and can easily be appended) in the [experiment configuration](src/experiment_config.jl) file.

> [!WARNING]
> Following the steps outlined below will result in data to be written to your disk.

To initialise a specific code-architecture configuration, indicate the corresponding configuration string in the [setup script](scripts/execution/dqc_setup.jl) and run

```
julia
] activate .
include("scripts/execution/dqc_setup.jl")
```

This stores the `NetworkingSpecifications` and `CodeParameters` to the `data/`folder. Using `DQCode.jl` as entrypoint to the $\texttt{DQCode}$ functionalities, the `scripts/execution/` folder contains the scripts needed to execute the entire $\texttt{DQCode}$ pipeline, making optimisations and simulations accessible from the REPL. 

For the optimisation/circuit search, we currently expose a Monte Carlo Tree Search [(MCTS)](src/mcts.jl) and warm-start [Genetic Algorithm](src/genetic.jl). These can be accessed via 

```
include("scripts/execution/dqc_mcts.jl")
```

and

```
include("scripts/execution/dqc_gott_ga.jl")
```

For the [DQC simulation](src/dqc_simulator.jl), individual noise sweeps can be specified in `scripts/execution/dqc_sim.jl` and executed via

```
include("scripts/execution/dqc_sim.jl")
```

Plotting functions are accessible via the `scripts/analysis` folder, allowing insight into [optimiser statistics](scripts/analysis/optimiser_evolution.jl) and [logical rate analysis](scripts/analysis/logical_rate.jl).

> [!IMPORTANT]
> We refer to `small' encoding circuits in the title explicitly, since the GA and MCTS in their current form don't scale well beyond the ~20 qubit regime, and since we have capped the number of flag qubits for the verification circuit synthesis to $75$ to not allow a disproportionately high verification overhead.


## Acknowledgements

$\texttt{DQCode}$ draws upon a rich collection of (general and quantum information) software, including
- `QuantumClifford`: https://github.com/QuantumSavory/QuantumClifford.jl, https://arxiv.org/abs/2512.16752
- `Munich Quantum Toolkit QECC` : https://github.com/munich-quantum-toolkit/qecc/tree/03a62ca1d3ccbe690265d6b5a7c59c5f72681793, https://arxiv.org/abs/2408.11894
- `POMDPs`: https://github.com/JuliaPOMDP/MCTS.jl 
- `Qiskit`: https://github.com/Qiskit/qiskit/blob/main/qiskit/synthesis/clifford/clifford_decompose_full.py
- `QuantikZ`: arXiv:1809.03842
- `HiGHS`: https://doi.org/10.1007/s12532-017-0130-5
- `LsqFit`: https://github.com/JuliaNLSolvers/LsqFit.jl?tab=readme-ov-file
- `Karlsruhe Hypergraph Partitioning (KaHyPar)`: https://kahypar.org, https://github.com/kahypar/KaHyPar.jl
> [!NOTE]
> The `POMDPs` library as well as the `KaHyPar.jl` library are licensed under a MIT "Expat" license. The `KaHyPar C++` library is licensed under the GPL License: https://github.com/kahypar/KaHyPar.jl?tab=License-1-ov-file.
> The `Qiskit` library is licensed under an Apache 2.0 license.

### Reproducibility Information

We chose to collect all hyperparameter data in the [configuration file](src/experiment_config.jl) for reproducibility. Upon executing the aforementioned $\texttt{DQCode}$ pipeline, you will automatically collect all data in the `data/` folder. For reference, we published all data used in the project to the folder `project_data/`. (TODO: Publish this folder)

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

Functionality has been tested for MacOS.

### Citation
```
author: {Tim Neumann}
title: {DQCode: Fault-Tolerant Zero State Preparation for Small qLDPC Codes on Near-Term Distributed Quantum Architectures}
year: {2026}
```