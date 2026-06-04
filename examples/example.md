# DQCode Pipeline
DQCode can be used to find optimised encoding circuits for logical $|0\rangle_L$ states of Quantum Error Correction CSS codes on distributed Type-II architectures. The $\texttt{scripts/}$ folder exposes the core functionalities of the search and evaluation pipeline. In this tutorial, we contextualise these scripts to demonstrate the pipeline for the Steane [[7,1,3]] code on a 2-core architecture.

> [!IMPORTANT] 
> The purpose of this markdown scripts is to show how the DQCode scripts work together, and how to navigate the repository. However, it is not devised as a stand-alone executable. We thus ask the user to view this file as a guideline of how to use the scripts provided in the $\texttt{scripts/}$ folder in the Julia REPL. To this end, we outline the pipeline, and at every step indicate which commands to run in the REPL. For the general setup instructions, please refer to the README.md file. Having executed these instructions, one can proceed with this tutorial.

Visualisations and results shown below are extracted from the data/ folder, to which DQCode automatically writes (Note: running the listed commands will write data to your disk.)

Throughout this tutorial, we assume that the user runs the respective indicated scripts, adapting the scripts per the code in questions. For this notebook, we set the experiment label in `/scripts/dqc_setup.jl` and all other scripts to `steane_4_3`.

## Setting up CSS code and Type-II architecture

> [!WARNING]
> Following the steps outlined below will result in data to be written to your disk.

In your terminal, navigate to your `DQCode` folder. Once you entered it, activate the Julia REPL by calling 
```
julia
] activate .
```
Leave the package manager by pressing the Del key. Now set the correct label in the [setup script](/scripts/execution/dqc_setup.jl). Back in your Julia REPL, run

```
include("scripts/execution/dqc_setup.jl")
```

The setup data has been saved to `DQCode/data/Steane/[4, 3]/`
We can check that the code parameters and network architecture have been correctly initialised.

```
using Serialization

steane_dir = joinpath(@__DIR__, "data/Steane/[4, 3]")

network_specs = deserialize( joinpath(steane_dir, "network_specs.jls"))
code_params = deserialize( joinpath(steane_dir, "code_params.jls"))

network_specs, code_params
```

This yields 
```
NetworkSpecifications([4, 3], 2, [2, 3, 6, 7, 1, 4, 5], [(1, 5), (2, 5), (3, 5), (4, 6), (5, 7), (6, 7)], [5, 1, 2, 6, 7, 3, 4], [1, 1, 1, 1, 2, 2, 2], [1, 2, 3, 4, 5, 6, 7], [8, 9], 7, 2, 1, 9)
```
and 

```
CodeParameters(QECCore.Steane7(), 3, Stabilizer 1×7, Stabilizer 7×7, [0 0 … 2 0; 0 0 … 1 0; … ; 2 0 … 2 0; 1 0 … 1 0], 7, 1, 3)
```
Let's inspect the mapping determined by the hypergraph partitioning algorithm further.

```
for (idx, core_size) in enumerate(network_specs.register_sizes)
    println("Core $idx contains qubits $(network_specs.mapping[ (idx>1 ? sum(network_specs.register_sizes[1:idx-1]) +1 : 1 ) : (idx>1 ? sum(network_specs.register_sizes[1:idx-1]) +1 : 1 ) + core_size-1]  )")
end
```
```
Core 1 contains qubits [2, 3, 6, 7]
Core 2 contains qubits [1, 4, 5]
```

## Qiskit baseline encoding

Next, we use the qiskit library to obtain a baseline encoding. The method exposed in qiskit is meant to synthesised general stabiliser tableaus on a monolithic all-to-all architecture, so we expect it to perform poorly when evaluated with respect to telegates.

Setting the correct label in the [Qiskit baseline script](/scripts/execution/dqc_qiskit_baseline.jl), we run

```
include("scripts/execution/dqc_qiskit_baseline.jl")
```

This saves the qiskit baseline data to `data/Steane/[4, 3]/qiskit_encoding`.
Here, you should fine the encoding circuit produced by Qiskit

<img src="../data/Steane/%5B4,%203%5D/qiskit_encoding/qiskit_encoding_circuit.png" alt="Qiskit encoding circuit" width="450"/>

You can also take a look at the gate counts here:
[Qiskit encoding circuit gate counts](/data/Steane/%5B4,%20%33%5D/qiskit_encoding/qiskit_encoding_stats.csv).

Similarly, you can perform the MQT-QECC baseline encoding by changing the label in the [MQT baseline script](/scripts/execution/dqc_mqt_baseline.jl) and running

```
include("scripts/execution/dqc_mqt_baseline.jl")
```

Here, we chose the `heuristic` method for the verification circuit, yielding the circuit 

<img src="../data/Steane/%5B4,%203%5D/mqt_encoding/mqt_encoding_circuit.png" alt="MQT encoding circuit" width="250"/>

## Monte Carlo Tree Search

Having explored the baselines, we now perform a Monte Carlo Tree Search for an efficient $|0\rangle_L$ encoding circuit. For this, set the correct label in the [GA script](/scripts/execution/dqc_mcts.jl) and set your the MCTS hyperparameters in the [experiment configuration file](/src/experiment_config.jl). In the current configuration, we use the hyperparameters

```
max_steps=15
fitness_weights=[1e6,1, 5, 2.5e4]
discount_factor=0.999
tableau_metric="jaccard"
reuse_tree=true
depth=3
n_iterations=5e4
exploration_constant=10.0
```

Then run 
```
include("scripts/execution/dqc_mcts.jl")
```

Let's take a look at the optimised circuit the genetic search produced:

<img src="../data/Steane/%5B4,%203%5D/mcts/MCTS_circuit.png" alt="MCTS encoding circuit" width="320"/>

How did the optimisation perform in a DQC setting? Find out in the [MCTS statistics data .csv file](/data/Steane/%5B4,%20%33%5D/mcts/mcts_stats.csv).

As the file reveals, we need four telegates for the distributed implementation!

## Genetic Search for efficient encoding circuits

As a second optimiser, it is time initialise our genetic algorithm, which is seeded with optimised Gottesman encoding circuits as warm-start circuit. 

For this, set the correct label in the [GA script](/scripts/execution/dqc_gott_ga.jl) and set your the genetic search hyperparameters in the [experiment configuration file](/src/experiment_config.jl) . In the current configuration, we use the hyperparameters

```
num_individuals=7500
num_generations=1500
max_len=100
mutation_rate=0.85
tournament_size=5
num_elite=1
fitness_weights=[1e4, 1, 10, 1e2]
tableau_metric='jaccard'
```

Then run 
```
include("scripts/execution/dqc_gott_ga.jl")
```

The genetic search can take some time to complete. Once it completes, it has saved all relevant data from the genetic search to the `/data/Steane/[4, 3]/warmstart_ga/` folder. We first inspect the Gottesman encoding circuit:

<img src="../data/Steane/%5B4,%203%5D/warmstart_ga/gott_encoding_circuit.png" alt="Gottesman encoding circuit" width="350"/>

Then the DQC-compiled version of it:

<img src="../data/Steane/%5B4,%203%5D/warmstart_ga/gott_circuit_dqc_compiled.png" alt="DQC-compiled Gottesman encoding circuit" width="300"/>

Finally, let's take a look at the optimised circuit the genetic search produced:

<img src="../data/Steane/%5B4,%203%5D/warmstart_ga/GA_circuit.png" alt="Warm-start GA encoding circuit" width="320"/>

How did the optimisation perform in a DQC setting? Find out in the [GA statistics data .csv file](/data/Steane/%5B4,%20%33%5D/warmstart_ga/warm_start_ga_stats.csv).

As the file reveals, we need two telegates for the distributed implementation!

We can also inspect how the optimiser reached its goal.

With the `/scripts/analysis/optimiser_evolution.jl` script, we can visualise the evolution of fitness and fidelities. For this, enter the code (here: "Steane") and the qpu_size (here: "[4, 3]") in the [optimiser_evolution script](/scripts/analysis/optimiser_evolution.jl) and run

```
include("scripts/analysis/optimiser_evolution.jl")
```




This produces 

<img src="../data/Steane/%5B4,%203%5D/warmstart_ga/optimisation_evolution.png" alt="Optimiser Evolution" width="400"/>


## DQC Simulation

Using the optimised GA circuit, we now run the DQC simulation. In 
the [DQC simulation script](/scripts/execution/dqc_sim.jl)
`scripts/dqc_sim.jl`, we can set the corresponding noise range. For the below results, we swept noise $p \in [5e-5,1e-3]$, with $p_{Bell} \in [1e-3, 5e-2]$, with $15$ values each. We have further set the telegate idle depth to $10$, the single-qubit error rate to $p/100$ and the idling rate to $p/10$.

For the simulation, we take the optimmised encoding circuit for the logical zero state, append a (flag-based) verification circuit from the Munich Quantum Toolkit QECC package, and then determine the logical error rate of the fault-tolerant encoding under circuit-level and Bell pair initialisation noise, mimicing realistic noise channels under telegate execution. In this example, we choose the `optimal` verification method from MQT-QECC.

Then we run 
```
include("scripts/execution/dqc_sim.jl")
```

This saves the simulation data to `/data/Steane/[4, 3]/simulation_FT'. Let us inspect the verification circuit that was appended to our optimised circuit in order to make it fault-tolerant.

<img src="../data/Steane/%5B4,%203%5D/simulation_FT/verification_circuit.png" alt="Verification circuit" width="150"/>

We can see that the verification circuit for the Steane code is fairly simple: we check for logical X errors via Z-type ancilla measurement. The canonical logical X operator is $X_3X_5X_6$, but recalling the quotient group equivalence, it is equivalent to $(X_3X_5X_6)(X_2X_3X_6X_7) = X_2X_5X_7$, where $X_2X_3X_6X_7$ is a stabiliser of the code. (The above circuit measures $X_2X_5X_7$).

This verification circuit guarantees that no low-weight errors propagate to harmful high-weight errors without being flagged by a -1 eigenstate measurement in either of the ancillas or the flag qubit.
as or the flag qubit.

Since in the DQC simulation, we simulate the execution of the raw encoding circuit and the verification circuit on the distributed architecture, let us take a look at the DQC executable circuit. This includes the circuit-level noise and Bell pair initialisation noise we wish to examine.

<img src="../data/Steane/%5B4,%203%5D/simulation_FT/DQC_circuit.png" alt="DQC circuit" width="750"/>

There's a lot going on in this circuit. Firstly, the first seven qubits are our data qubits. The (noiseless) SWAP operations at the very beginning of the circuit serve to visualise that we have mapped our qubits in a specific way (based on the hypergraph partitioning determined by distributed QEC cycles). (The SWAP operations are noiseless, since in a real experiment, one would simply relabel the qubits at the beginning of the circuit; this is different from SWAP operations that occur in the middle of the circuit, which might necessitate tele-operations). Additionally, we find two communication qubits (qubits $8$ and $9$), which enable our telegates between the two cores. Earlier, we identified that two CNOT gates are telegates: CNOT(1,3) and CNOT(4,7). These telegates have been replaced with the corresponding the EJPP protocol, making use of the respective communication qubits. The $10^{\text{th}}$ qubit stems from the verification circuit, and is used to discard runs that violated the FT criterion (and resulted in harmful errors from low-weight errors). Also, there are classical lines storing the measurement information from the circuit.

The full circuit, including noise-free stabiliser measurements, is stored as a [.tex file](/data/Steane/[4,%203]/simulation_FT/full_circuit.tex).


We can inspect the data stored by the above simulation in the 
[simulation data file](/data/Steane/%5B4,%20%33%5D/simulation_FT/dqc_sim_data.csv).

With the `/scripts/analysis/logical_rate.jl` script, we can now visualise this data for further insight. For this, enter the code (here: "Steane") and the qpu_size (here: "[4, 3]") in the [logical_rate script](/scripts/analysis/logical_rate.jl) and run

```
include("scripts/analysis/logical_rate.jl")
```

This saves a log-log plot containing the logical initialisation error scaling per physical initialisation error rate `p`

<img src="../data/Steane/%5B4,%203%5D/simulation_FT/qec_threshold.png" alt="log-log plot" width="350"/>

as well as 2d heatmaps including the Bell state initialisation error probability `p_Bell` to `/data/Steane/[4, 3]/simulation_FT`.

<img src="../data/Steane/%5B4,%203%5D/simulation_FT/2d_heatmap.png" alt="2d" width="450"/>

<img src="../data/Steane/%5B4,%203%5D/simulation_FT/2d_heatmap_ratio.png" alt="2d ratio" width="450"/>


We can also create other visualisations based on the data. 

(TODO: Add latest results and other visualisation)

## Resource estimation
Finally, let's analyse the resources that the FT encoding circuit consumes. For this, set the correct label in the [GA script](/scripts/execution/dqc_resource.jl) and run

```
include("scripts/execution/dqc_resource.jl")
```

This saves data into our `/data/Steane/[4, 3]/simulation_FT` folder, containing the number of ancillas, depth, gate counts and effective required QPU sizes (including ancilla qubits) for both the [FT encoding circuit](/data/Steane/%5B4,%20%33%5D/simulation_FT/resources_info_circ.csv) and a [1-round distributed stabiliser measurement approach](/data/Steane/%5B4,%20%33%5D/simulation_FT/resources_info_circ.csv).

For the Steane-code, we find that the FT encoding circuit needs less telegates and also less telegate layers, while only using one additional ancilla, compared to one round of distributed stabiliser measurements (which would not be FT per se). Also, we find that the ancilla from verification can be reused as an ancilla for distributed stabiliser measurement!

