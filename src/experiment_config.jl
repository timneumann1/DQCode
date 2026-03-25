module ExperimentConfig

using ..Types
using QECCore
using QECCore: Steane7, QuantumTannerGraphProduct, CyclicQuantumTannerGraphProduct, Triangular488, distance
using ..TrivariateBicycleCode

#distributed_qec_code = Steane7()
#distributed_qec_code = Shor9()
distributed_qec_code = TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)])
# distributed_qec_code = BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)])
# distributed_qec_code = Triangular488(5)
# distributed_qec_code = BivariateBicycleViaCirculantMat(12, 6, [(:x, 3), (:y, 1), (:y, 2)], [(:y, 3), (:x, 1), (:x, 2)])

#type_two_register_sizes = [4,3]
#type_two_register_sizes = [3,3,3]
type_two_register_sizes = [3,3,3,3]
#type_two_register_sizes = [3,3,3,3,3,3]


##### Optimisation Parameters #####

opt_params = OptimisationParameters(
    "jaccard" # tableau distance metric
) 

genetic_params = GeneticParameters(
    300, # individuals
    100, # generations
    100, # max length of (raw) circuit individual
    0.8,  # mutation rate
    5, # tournament size
    0.5, # selection_ratio
    1, # num_elite
    false, # warm_start
    [1,3,10] # fitness weights
    )

mcts_params = MCTSParameters(
    20, # steps before termination
    2000, #iterations
    1.414265 # exploration constant
)

end