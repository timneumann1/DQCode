# experiment_config.jl

"""
Stores the code-network configurations to be investigated, with `code_architecture_setup` containing
the CSS code and GPU sizes of a setup. 
Also allows to define hyperparameters for genetic search and MCTS for each code-network configuration.
"""
module ExperimentConfig

export experiment_configurations

using ..Types: GeneticParameters, MCTSParameters
using ..TrivariateBicycleCode
using ..SymplecticDoubleCode

using QECCore: Steane7, Shor9, AbstractCSSCode, BivariateBicycleViaCirculantMat, Triangular488

function experiment_configurations()

    # ----------------------------- Type Setup for configurations -----------------------------
    SetupCombo = @NamedTuple{code::AbstractCSSCode, qpu_sizes::Vector{Int}}
    code_architecture_setup = Dict{String, SetupCombo}()
    genetic_params = Dict{String, GeneticParameters}()
    mcts_params = Dict{String, MCTSParameters}()

    # -----------------------------------------------------------------------------------------
    # ----------------- Codes, Architectures and optimisation parameters ---------------------- 
    # -----------------------------------------------------------------------------------------

    # ------- Steane7 [[7,1,3]]  - [4,3] ---------

    code_architecture_setup["steane_4_3"] = ( Steane7(), [4,3] )
    genetic_params["steane_4_3"] = GeneticParameters(7500, 1500, 100, 0.85, 5, 1, [1e4, 1, 10, 1e2], "jaccard")
    mcts_params["steane_4_3"] =  MCTSParameters(15, [1e6,1, 5, 2.5e4], 0.999, "jaccard", true, 3, 5e4, 10.0)

    # ------- Shor9 [[9,1,3]] - [3,3,3] ---------

    code_architecture_setup["shor_3_3_3"] = ( Shor9(), [3,3,3] )
    genetic_params["shor_3_3_3"] = GeneticParameters(7500, 1500, 100, 0.85, 5, 1, [1e4, 1, 10, 1e2], "jaccard")
    mcts_params["shor_3_3_3"] =  MCTSParameters(15, [1e6,1,5,1e2], 0.9, "jaccard", true, 5, 1e5, 3.5)

    # ------- Trivariate Bicycle [[12,2,3]]- [3,3,3,3] ---------

    code_architecture_setup["trivariate_3_3_3_3"] = ( TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]), [3,3,3,3] )
    genetic_params["trivariate_3_3_3_3"] = GeneticParameters(25000, 1000, 100, 0.85, 5, 1, [1e4, 1, 10, 5e2], "jaccard")
    mcts_params["trivariate_3_3_3_3"] = MCTSParameters(30, [1e6,1,5,1e4], 0.99, "jaccard", false, 3, 5e4, 10.0)

    # ------- Trivariate Bicycle [[12,2,3]] - [4,4,4] ---------

    code_architecture_setup["trivariate_4_4_4"] = ( TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]), [4,4,4] )
    genetic_params["trivariate_4_4_4"] = GeneticParameters(25000, 2500, 100, 0.85, 5, 1, [1e4, 1, 10, 5e2], "jaccard")
    mcts_params["trivariate_4_4_4"] = MCTSParameters(30, [1e6,1,5,2.5e4], 0.999,"jaccard", false, 3, 5e4, 10.0)
    
    # ------- Trivariate Bicycle [[12,2,3]] - [6,6] ---------

    code_architecture_setup["trivariate_6_6"] = ( TrivariateBicycleViaCirculantMat(2, 3, [(:x, 1), (:y, 2)],[(:x, 0), (:z, 4)]), [6,6] )
    genetic_params["trivariate_6_6"] = GeneticParameters(25000, 2500, 100, 0.85, 5, 1, [1e4, 1, 10, 5e2], "jaccard")
    mcts_params["trivariate_6_6"] = MCTSParameters(50, [1e6,1,5,4e4], 0.99,"jaccard", true, 2, 5e4, 5.0)
    
    # ------- Color [[17,1,5]] - [8,9] ---------

    code_architecture_setup["color_8_9"] = ( Triangular488(5), [8,9] )
    genetic_params["color_8_9"] = GeneticParameters(15000, 2500, 100, 0.85, 5, 1, [2e4, 1, 10, 2e3], "jaccard")
    #mcts_params["color_8_9"] = MCTSParameters(50, [1e6,1,5,2e3], 0.99, "jaccard", false, 3, 5e4, 5.0)
    
    # ------- Bivariate Bicycle [[18,4,4]] - [3,3,3,3,3,3] ---------

    code_architecture_setup["bivariate_3_3_3_3_3_3"] = ( BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]), [3,3,3,3,3,3] )
    genetic_params["bivariate_3_3_3_3_3_3"] = GeneticParameters(15000, 5000, 100, 0.85, 5, 1, [3e4, 1, 10, 1e3], "jaccard")
    mcts_params["bivariate_3_3_3_3_3_3"] = MCTSParameters(70, [1e6,1,5,5e3], 0.999, "jaccard", false, 3, 5e4, 10.0)
    
    # ------- Bivariate Bicycle [[18,4,4]] - [6,6,6] ---------

    code_architecture_setup["bivariate_6_6_6"] = ( BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]), [6,6,6] )
    genetic_params["bivariate_6_6_6"] = GeneticParameters(15000, 5000, 100, 0.85, 5, 1, [4e4, 1, 10, 1e3], "jaccard")
    mcts_params["bivariate_6_6_6"] = MCTSParameters(70, [1e6,1,5,2e3], 0.999, "jaccard", false, 3, 5e4, 10.0)

    # ------- Bivariate Bicycle [[18,4,4]] - [9,9] ---------

    code_architecture_setup["bivariate_9_9"] = ( BivariateBicycleViaCirculantMat(3, 3, [(:x, 0), (:x, 1), (:y, 1)], [(:y, 0), (:x, 2), (:y, 2)]), [9,9] )
    genetic_params["bivariate_9_9"] = GeneticParameters(15000, 5000, 100, 0.85, 5, 1, [5e4, 1, 10, 1e3], "jaccard")
    mcts_params["bivariate_9_9"] = MCTSParameters(85, [1e6,1,5,1e3], 0.999, "jaccard", false, 3, 5e4, 10.0)


    # ------- Symplectic Double Code [[30,6,5]] - [15, 15] ---------

    code_architecture_setup["sd_15_15"] = (  SymplecticDouble(), [15,15] )

    # ------- Trivariate Bicycle [[30,4,5]] - [15,15] ---------------

    code_architecture_setup["trivariate_15_15"] = ( TrivariateBicycleViaCirculantMat(5, 3, [(:x, 4), (:x, 2)], [(:x, 1), (:x, 2), (:y, 1), (:z, 2), (:z, 3)]), [15,15] )
  
    # ------- Color [[31,1,7]] - [15,16] ---------

    code_architecture_setup["color_15_16"] = ( Triangular488(7), [15,16] )

    # ------- Bivariate Bicycle [[36,4,6]] - [12,12,12] ---------

    code_architecture_setup["bivariate_12_12_12"] = ( BivariateBicycleViaCirculantMat(3, 6, [(:x, 1), (:y, 2), (:y, 3)], [(:x, 0), (:y, 1), (:x, 2)]), [12,12,12] )

    #  ------- Bivariate Bicycle - [[144,12,12]] ---------

    code_architecture_setup["bivariate_12_12_12_12_12_12_12_12_12_12_12_12"] = ( 
        BivariateBicycleViaCirculantMat(12, 6, [(:x, 3), (:y, 1), (:y, 2)], [(:y, 3), (:x, 1), (:x, 2)]),
        [12,12,12,12,12,12,12,12,12,12,12,12]
    )
    
    return code_architecture_setup, genetic_params, mcts_params
end


end

