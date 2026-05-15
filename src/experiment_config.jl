module ExperimentConfig

using ..Types: GeneticParameters, MCTSParameters
using QECCore: Steane7, Shor9, AbstractCSSCode, BivariateBicycleViaCirculantMat#, distance
using ..TrivariateBicycleCode

# This file indicates the code-network configurations we are investigating, with code and gpu_sizes indicating the differnet setups,
# and the parameters being used as defaults for the circuit search (can be changed in scripts).

export experiment_configurations

function experiment_configurations()
    ConfigKeyType = String
    ConfigValueType = @NamedTuple{
        code::AbstractCSSCode,
        qpu_sizes::Vector{Int},
        genetic_params::GeneticParameters, 
        mcts_params::MCTSParameters, 
    }
    configs = Dict{ConfigKeyType, ConfigValueType}()

    # ------- Steane7 - [4,3] ---------

    configs["steane_4_3"] = (
        Steane7(),
        [4,3],
        GeneticParameters(7500, 1500, 100, 0.85, 5, 0.5, 1, [1e4, 1, 10, 1e2], "jaccard"),
        MCTSParameters(15, [1e6,1, 5, 2.5e4], 0.999, "jaccard", true, 3, 5e4, 10.0)
        )
    
    # ------- Shor9 - [3,3,3] ---------

    configs["shor_3_3_3"] = (
        Shor9(),
        [3,3,3],
        GeneticParameters(7500, 1500, 100, 0.85, 5, 0.5, 1, [1e4, 1, 10, 1e2], "jaccard"),
        MCTSParameters(15, [1e6,1,5,1e2], 0.9, "jaccard", true, 5, 1e5, 3.5)
        )

    # ------- Trivariate Bicycle - [3,3,3,3] ---------

    configs["trivariate_3_3_3_3"] = (
        TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]),
        [3,3,3,3],
        GeneticParameters(25000, 1000, 100, 0.85, 5, 0.5, 1, [1e4, 1, 10, 5e2], "jaccard"), 
        MCTSParameters(30, [1e6,1,5,1e4], 0.99, "jaccard", false, 3, 5e4, 10.0),#, 5, 1e5, 10)
        )
    
    # ------- Trivariate Bicycle - [4,4,4] ---------

    configs["trivariate_4_4_4"] = (
        TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]),
        [4,4,4],
        GeneticParameters(25000, 2500, 100, 0.85, 5, 0.5, 1, [1e4, 1, 10, 5e2], "jaccard"),
        MCTSParameters(30, [1e6,1,5,2.5e4], 0.999,"jaccard", false, 3, 5e4, 10.0)
        )

    # ------- Trivariate Bicycle - [6,6] ---------

    configs["trivariate_6_6"] = (
        TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]),
        [6,6],
        GeneticParameters(25000, 2500, 100, 0.85, 5, 0.5, 1, [1e4, 1, 10, 5e2], "jaccard"),
        MCTSParameters(30, [1e6,1,5,4e4], 0.99,"jaccard", true, 2, 5e4, 5.0)
        )

    # ------- Bivariate Bicycle [[18,4,4]] - [3,3,3,3,3,3] ---------

    configs["bivariate_3_3_3_3_3_3"] = (
        BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]),
        [3,3,3,3,3,3],
        GeneticParameters(15000, 5000, 100, 0.85, 5, 0.5, 1, [3e4, 1, 10, 1e3], "jaccard"),
        MCTSParameters(70, [1e6,1,5,5e3], 0.999, "jaccard", false, 3, 5e4, 10.0)
        )
    
    # # # ------- Bivariate Bicycle [[18,4,4]] - [6,6,6] ---------

    configs["bivariate_6_6_6"] = (
        BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]),
        [6,6,6],
        GeneticParameters(15000, 5000, 100, 0.85, 5, 0.5, 1, [4e4, 1, 10, 1e3], "jaccard"),
        MCTSParameters(70, [1e6,1,5,2e3], 0.999, "jaccard", false, 3, 5e4, 10.0)
        )

    # # # ------- Bivariate Bicycle [[18,4,4]] - [9,9] ---------

    configs["bivariate_9_9"] = (
        BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]),
        [9,9],
        GeneticParameters(15000, 5000, 100, 0.85, 5, 0.5, 1, [5e4, 1, 10, 1e3], "jaccard"),
        MCTSParameters(85, [1e6,1,5,1e3], 0.999, "jaccard", false, 3, 5e4, 10.0)
        )
    
    # # # ------- Bivariate Bicycle [[36,4,6]] - [6,6,6,6,6,6] ---------

    # configs["bivariate_6_6_6_6_6_6"] = (
    #     BivariateBicycleViaCirculantMat(3, 6, [(:x, 1), (:y, 2), (:y, 3)], [(:x, 0), (:y, 1), (:x, 2)]),
    #     [6,6,6,6,6,6],
    #     GeneticParameters(5000, 5000, 100, 0.85, 5, 0.5, 1, [5e4, 1, 10, 1e3], "jaccard"),
    #     MCTSParameters(70, [1e5, 1, 5, 1e2], 0.999, "jaccard", false, 3, 1e6, 10.0),
    #     )

    # # ------- Bivariate Bicycle - [[144,12,12]] ---------

    # configs["bivariate_12"] = (
    #     BivariateBicycleViaCirculantMat(12, 6, [(:x, 3), (:y, 1), (:y, 2)], [(:y, 3), (:x, 1), (:x, 2)]),
    #     [12,12,12,12,12,12,12,12,12,12,12,12],
    #     GeneticParameters(15000, 5000, 100, 0.85, 5, 0.5, 1, [5e4, 1, 10, 1e3], "jaccard"),
    #     MCTSParameters(2000, [1e6, 1, 5, 1e2], 0.999, "jaccard", 2, 1e4, 10.0),
    #     joinpath(@__DIR__, "results", string(code_dirname(BivariateBicycleViaCirculantMat(12, 6, [(:x, 3), (:y, 1), (:y, 2)], [(:y, 3), (:x, 1), (:x, 2)]))), string([12,12,12,12,12,12,12,12,12,12,12,12]))
    #     )

    return configs

end


end

