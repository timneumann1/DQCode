module DQCodePrep

include("DQCircuitSearch.jl")
include("types.jl")
include("trivariate_bicycle_code.jl")
include("helper.jl")
include("experiment_config.jl")

#include("circsim.jl")
#include("logical_encoding.jl")
#include("plots.jl")
#include("dtsimulation.jl")
#include("genetic.jl")
#include("mcts.jl")
include("qec_tools.jl")
include("dqc_state_prep_sim.jl")
#include("deep_q.jl")
#include("parameters.jl")



using QECCore
using QECCore: distance
using .TrivariateBicycleCode
using .DQCircuitSearch

using Random
using QuantumClifford
using QuantumClifford: MixedDestabilizer, Stabilizer, Tableau, stabilizerview, logicalxview, logicalzview, canonicalize_rref!, tab, AbstractOperation
using QuantumClifford.ECC: DistanceMIPAlgorithm
using Serialization
using StatsBase


using .Helper: tableau_to_bitmatrix, data_qubit_partitioning, perm_to_transpositions, create_lookup_array, verify_success, execute_circuit, qc_circuit_to_qasm
using .Types

using .QECTools
using .DQCLogicalStatePrepSimulator

using .ExperimentConfig: experiment_configurations

using PyCall
#np = pyimport("numpy")

# py"""
# import random, numpy as np, z3
# random.seed(0)
# np.random.seed(0)
# z3.set_param("smt.random_seed", 0)
# z3.set_param("sat.random_seed", 0)
# """



# logging = pyimport("logging")
# logging.basicConfig(
#     level=logging.INFO,
#     format="%(asctime)s %(name)s %(levelname)s: %(message)s"
# )
# logging.getLogger("mqt.qecc").setLevel(logging.INFO)

# #CSSCode = pyimport("mqt.qecc").CSSCode
# mqt_synthesis = pyimport("mqt.qecc.circuit_synthesis")

# # sim = pyimport("mqt.qecc.simulation")
# # noise = pyimport("mqt.qecc.noise")

qiskit = pyimport("qiskit")
qi = pyimport("qiskit.quantum_info")
qasm2 = pyimport("qiskit.qasm2")
#synth = pyimport("qiskit.synthesis")
plt = pyimport("matplotlib.pyplot")

# #heuristic_prep_circuit         = cs.heuristic_prep_circuit
# gate_optimal_verification_circuit = mqt_synthesis.gate_optimal_verification_circuit
# mqt_circuits = mqt_synthesis.circuits# pyimport("mqt.qecc.circuits")
# mqt_state_prep = mqt_synthesis.state_prep
#VerificationNDFTStatePrepSimulator = cs.VerificationNDFTStatePrepSimulator
#CircuitLevelNoiseIdlingParallel    = cs.CircuitLevelNoiseIdlingParallel


function run_dqc_state_prep(exp_label::String)
    Random.seed!(42) 
    configs = experiment_configurations()

    # data_circuit = Vector{AbstractOperation}()
    # verification_qasm = ""
    cfg = ""

    if haskey(configs, exp_label)
        cfg = configs[exp_label]
        if !isfile(joinpath(cfg.folder, "network_specs.jls")) || !isfile(joinpath(cfg.folder, "code_params.jls"))
            error("The serialized specification and parameter files for this experiment are missing. Please run create_code_network_data($exp_label).")
        end
    else
        error("The configuration label $exp_label was not found. Please add the respective data to the configuration file first.")
    end

    print(cfg.folder)
    network_specs = deserialize( joinpath(cfg.folder, "network_specs.jls"))
    code_params = deserialize( joinpath(cfg.folder, "code_params.jls"))
    filepath = "$(cfg.folder)/dqc_comp_gott>GA/GA_circuit.jls" 


    circuit = deserialize(filepath)
    data_circuit = copy(circuit)

    println("Initial circuit: $circuit")

    qasm = qc_circuit_to_qasm(circuit)
    
    println("QASM version: $qasm")

    # Procedure: We pass a qasm string with the optimised encoding circuit, and get back a qasm string with the verification (in MQT QECC, the result was a qiskit circuit, which was then converted to qasm in order to make the output stream usable, and then 
    # converted to qiskit here again)
    
    
    
    
    # qiskit_circ = qasm2.loads(qasm) # loads the string
    # print(qiskit_circ)
    # #qiskit_circ.draw(output="mpl", initial_state=true, fold=-1, scale=0.4)
    # #plt.show()

    # cnot_circ = mqt_circuits.CNOTCircuit.from_qiskit_circuit(qiskit_circ, init_all = true)
    
    
    # print("CNOTCIRCUIT: $(cnot_circ.cnots)")


    # # get code from cnot circ
    # css_code = cnot_circ.get_code()
    # # and check whether the stabilisers and logicals match
    # println("X checks: $(css_code.Hx)")
    # println("Z checks: $(css_code.Hz)")
    # println("Logical X: $(css_code.Lx)")
    # println("Logical Z: $(css_code.Lz)")

    # # everthing is now collapsed in the Z check matrix (Z stabilisers and logical Zs), since they are equaivlent as stabilisers of the desired state
    # # and we can measrue both to obtain the correct one (of course cofrectinng for hook errors); our initial definition of logical still applies and
    # # can be used for the noiseless procedure

    # println(code_params.distance//2)
    # t = Int(floor(code_params.distance/2))
    # faulty_prep_circuit = mqt_state_prep.FaultyStatePrepCircuit(cnot_circ,t,t)
    # println("\n\n\n\n\n FAULT circuit $(faulty_prep_circuit.circ.cnots)\n\n\n")

    verification_circ_qasm = readchomp(`/Users/tim/Tim/projects/mqt/qecc/.venv/bin/python3 /Users/tim/Tim/projects/mqt/qecc/scripts/verification_circuit.py $qasm $(code_params.distance)`)
    #verification_circ_qasm = gate_optimal_verification_circuit(faulty_prep_circuit)#, min_timeout=600,max_timeout=3600)
    print("Verification QASM: \n$verification_circ_qasm")
    verification_circ = qasm2.loads(verification_circ_qasm) 

    println("Verification Cirucit: $verification_circ")
    #println("DAAATA: $(verification_circ.data)")


    #verification_circ.draw(output="mpl", initial_state=true, fold=-1, scale=0.4)
    #plt.show()


    #qasm_prog = qasm2.dumps(verification_circ) # dumps writes to a string
    # qasm program could be reordered but should still be intact:
    #println("Final QASM Program orig + verification: $qasm_prog")

    #verification_qasm = copy(qasm_prog)

    quantum_clifford_verification_circ = Vector{AbstractOperation}()
    # We actually parse throught the qiskit circuit now, which is more convenient than the QASM string

    #um_z_anc = filter verification_qasm to find line with qreg z_anc and then extract the number
    # num_x_anc = filter verification_qasm to find line with qreg x_anc and then extract the number
    # num_flag = filter verification_qasm to find line with qreg flag and then extract the number

    ancilla_data_interactions = Dict{Tuple{String, Int}, Vector{Int}}() # will contain ancilla qubits, e.g., Z_anc, 1 as keys, and interacting data qubits as values
    println(typeof(verification_circ))
    println(verification_circ.data)
    #registers = [reg for reg in verification_circ.qregs]
    num_q = 0
    num_z_anc = 0
    num_x_anc = 0
    num_flags = 0
    
    for reg in verification_circ.qregs
        if reg.name == "q"
            num_q = reg.size
        elseif reg.name == "z_anc"
            num_z_anc = reg.size
        elseif reg.name == "x_anc"
            num_x_anc = reg.size
        elseif reg.name == "flag"
            num_flags = reg.size # this is an upper bound, as set in MQT QECC
        end
    end
    @assert num_q == network_specs.num_data_qubits
    #num_ancillas = num_z_anc + num_x_anc + num_flags
    # Assume that we always add z_anc first, then x_anc, then flags in register ordering
    # num_z_anc = verification_circ.qregs[2]
    # num_x_anc = verification_circ.qregs[3]
    # num_flags = verification_circ.qregs[4]
    #print(num_q, num_z_anc, num_x_anc, num_flags)
    
    for instruction in verification_circ.data
        #println(verification_circ.data)
        gate = instruction.operation.name
        
        if gate == "h"
            
            qubits = instruction.qubits
            bit_info = verification_circ.find_bit(qubits[1])
            reg_name = String(bit_info.registers[1][1].name)
            if reg_name == "q"
                continue
            end
            index = Int(bit_info.index)
            
            println(gate, index, reg_name)
            # qc_index = 0 
            # if reg_name == "z_anc"
            #     qc_index = num_q + index
            # elseif reg_name == "x_anc"
            #     qc_index = num_q + num_z_anc + index
            # elseif reg_name == "flag"
            #     qc_index = num_q + num_z_anc + num_x_anc + index
            # end

            push!(quantum_clifford_verification_circ, sHadamard(index+network_specs.num_comm_qubits+1))

        elseif gate == "cx"
            qubits = instruction.qubits
            control_info = verification_circ.find_bit(qubits[1])
            control_reg_name = String(control_info.registers[1][1].name)
            control = Int(control_info.index)
            target_info = verification_circ.find_bit(qubits[2])
            target_reg_name = String(target_info.registers[1][1].name)
            target = Int(target_info.index)

            # target_reg_name holds the name of the register
            # target_index_reg holds the index WITHIN this reigster

            # the bitinfo index refers to num_data qubits (not comm qubut which are indexed in between in our architecture, but excluded in qskit) + 
            # whatever index in the 3 ancilla registers we have
            if control_reg_name == "q" && target_reg_name == "q"
                continue
            elseif control_reg_name == "q" 
                #println(target)
                if target_reg_name == "z_anc"
                    target_index_reg = target-network_specs.num_data_qubits#num_q + target
                elseif target_reg_name == "x_anc"
                    target_index_reg = target-network_specs.num_data_qubits-num_z_anc #target_index = num_q + num_z_anc + index
                elseif target_reg_name == "flag"
                    target_index_reg = target-network_specs.num_data_qubits-num_z_anc-num_x_anc#target_index = num_q + num_z_anc + num_x_anc + index
                end
                key = (target_reg_name,target_index_reg+1)
                push!(get!(ancilla_data_interactions, key, Int[]), network_specs.register_lookup_array[network_specs.inv_map[control+1]])
                #ancilla_data_interactions[] = network_specs.register_lookup_array[network_specs.inv_map[control+1]]
                println("Target: $target, network_specs.num_comm_qubits: $(network_specs.num_comm_qubits)")
                push!(quantum_clifford_verification_circ, sCNOT(control+1,target+network_specs.num_comm_qubits+1))
            elseif target_reg_name == "q" 
                
                if control_reg_name == "z_anc"
                    control_index_reg = control-network_specs.num_data_qubits
                elseif control_reg_name == "x_anc"
                    control_index_reg = control-network_specs.num_data_qubits-num_z_anc#num_q + num_z_anc + control
                elseif control_reg_name == "flag"
                    control_index_reg =  control-network_specs.num_data_qubits-num_z_anc-num_x_anc# num_q + num_z_anc + num_x_anc + control
                end
                key = (control_reg_name,control_index_reg+1)
                push!(get!(ancilla_data_interactions, key, Int[]), network_specs.register_lookup_array[network_specs.inv_map[target+1]])
                #ancilla_data_interactions[(control_reg_name,control_index_reg+1)] = network_specs.register_lookup_array[network_specs.inv_map[target+1]]
                push!(quantum_clifford_verification_circ, sCNOT(control+network_specs.num_comm_qubits+1,target+1))
            else # both are ancilla qubits
                if control_reg_name == "z_anc"
                    control_index_reg = control-network_specs.num_data_qubits
                elseif control_reg_name == "x_anc"
                    control_index_reg = control-network_specs.num_data_qubits-num_z_anc
                elseif control_reg_name == "flag"
                    control_index_reg = control-network_specs.num_data_qubits-num_z_anc-num_x_anc
                end
                if target_reg_name == "z_anc"
                    target_index_reg = target-network_specs.num_data_qubits
                elseif target_reg_name == "x_anc"
                    target_index_reg = target-network_specs.num_data_qubits-num_z_anc
                elseif target_reg_name == "flag"
                    target_index_reg = target-network_specs.num_data_qubits-num_z_anc-num_x_anc
                end
                control_key = (control_reg_name,control_index_reg+1)
                target_key = (target_reg_name,target_index_reg+1)
                # If both qubits are ancilla qubits, we would still like to record information about their mapping, so we simply map the current mode of the counterpart ancilla
                control_interactions = get!(ancilla_data_interactions, control_key, Int[])
                control_mode = isempty(control_interactions) ? 1 : mode(control_interactions)

                target_interactions = get!(ancilla_data_interactions, target_key, Int[])
                target_mode = isempty(target_interactions) ? 1 : mode(control_interactions)

                push!(get!(ancilla_data_interactions, control_key, Int[]), target_mode )
                push!(get!(ancilla_data_interactions, target_key, Int[]), control_mode )

                push!(quantum_clifford_verification_circ, sCNOT(control+network_specs.num_comm_qubits+1,target+network_specs.num_comm_qubits+1))
                    
                    #                push!(get!(ancilla_data_interactions, key, Int[]), network_specs.register_lookup_array[network_specs.inv_map[target+1]])

            end

        elseif gate=="measure" 

            qubits = instruction.qubits
            println("\n\n $qubits")
            bit_info = verification_circ.find_bit(qubits[1])
            reg_name = String(bit_info.registers[1][1].name)
            @assert reg_name != "q" # only ancillas will be measured
    
            index = Int(bit_info.index)
            #println("INDEX $index")
            #println(gate, control, target)
            # q_index = 0 
            # b_index = 0
            # if reg_name == "z_anc"
            #     q_index = num_q + index
            #     b_index = network_specs.num_comm_qubits + index
            #     println("\n\n\n\n\n The b index is $b_index = $index ")
            # elseif reg_name == "x_anc"
            #     q_index = num_q + num_z_anc + index
            #     b_index = network_specs.num_comm_qubits + num_z_anc + index
            # elseif reg_name == "flag"
            #     q_index = num_q + num_z_anc + num_x_anc +  index
            #     b_index = network_specs.num_comm_qubits + num_z_anc + num_x_anc + index
            # end
            push!(quantum_clifford_verification_circ, sMRZ(index+network_specs.num_comm_qubits+1, index-network_specs.num_data_qubits+network_specs.num_comm_qubits+1))
            
        end
        
        
    end
    #     see whether z, x anc or flag is part of it and if yes which
    #         determine core of qubit that maps to iff cx is the gate: n.register_lookup_array[n.inv_map[op.q1]]
    #         add entry to dictionary
    #     add QC operation sCNOT between unmpapped q1 and num_data_and_comm_qubits + i
    #     elseif operation is h
    #                 add QC operation H to num_data_and_comm_qubits + i
    #             elseif operation is measure
    #                 add QC operation sMZ from num_data_and_comm_qubits + i to num_comm_qubits + i
    #             end


    

        # core_mapping of size num_z_anc+num_x_anc + num_flag # ancillas will be listed lexicographical after data and comm qubits, but we store the core that the ancilla is actually mapped to (for telegates)

        #     choose the core with the max occurences fr the ancilla to be placed and store it in core_mapping

        

    
    # for i in 1:num_z 
    #     cores_addressed
    #     for line in verification_qasm
    #         if anc in line with index in line
    #             if operation is cx
    #                 determine core of qubit that maps to iff cx is the gate: n.register_lookup_array[n.inv_map[op.q1]] 
    #                 add core to list
    #                 add QC operation sCNOT between unmpapped q1 and num_data_and_comm_qubits + i
    #             elseif operation is h
    #                 add QC operation H to num_data_and_comm_qubits + i
    #             elseif operation is measure
    #                 add QC operation sMZ from num_data_and_comm_qubits + i to num_comm_qubits + i
    #             end
    #         end
    #     end
    #     choose the core with the max occurences fr the ancilla to be placed and store it in core_mapping
    # end
            
    # for i in 1:num_x
    #     cores_addressed
    #     for line in verification_qasm
    #         if anc in line with index in line
    #             if operation is cx
    #                 determine core of qubit that maps to iff cx is the gate: n.register_lookup_array[n.inv_map[op.q1]] 
    #                 add core to list
    #                 add QC operation sCNOT between unmpapped q1 and num_data_and_comm_qubits + num_z_anc + i
    #             elseif operation is h
    #                 add QC operation H to num_data_and_comm_qubits +  num_z_anc +_i
    #             elseif operation is measure
    #                 add QC operation sMZ from num_data_and_comm_qubits + i to num_comm_qubits + num_z_anc + i
    #             end
    #             end
    #         end
    #     end
    #     choose the core with the max occurences fr the ancilla to be placed and store it in core_mapping
    # end
           

    # for i in 1:num_flag
    #     for line in verification_qasm
    #         if anc in line with index in line
    #             if operation is cx
    #                 note down core of the ancilla that flag interacts with
    #                 if z_anc
    #                     add sCNOT from num_data_and_comm_qubits + index of z_anc to num_data_and_comm_qubits + num_z_anc + x_anc  + i
    #                 elseif x_anc
    #                     add sCNOT from num_data_and_comm_qubits + num_z_anc index of x_anc to num_data_and_comm_qubits + num_z_anc + x_anc  + i
    #                 end
    #             elseif if operation is h
    #                 add H to num_data_and_comm_qubits + num_z_anc + num_x_anc + i 
    #             elseif operation is measure
    #                 add QC operation sMZ from num_data_and_comm_qubits + num_z_anc + i to num_comm_qubits + num_z_anc + num_x_anc + i
    #             end
    #         end
    #     end
    #     determine ancilla that interacts with it (should onlt be one) and add it to the same core in core_mapping
    #     if lenght(cores) >1 
    #         throw error
    #     else
    #         core_mapping - core[1]
    # end


        
    # provide data_circuit and verifation_circuit to dqc_state_prep (separetely is ok)
    # Data qubits should continue to experience dephasing, but anciall qubits should only be initialised right before they are used
    # dqc_circ =  DQCLogicalStatePrepSimulator.dqc_state_prep(data_circuit, quantum_clifford_verification_circ, code_params, network_specs, noise_model)


    # determine the best placement for the ancilla qubits, noting which cnots have to be applied to them, and noting the index of the measueretn (on which we
    #later base the discarding)
    # Then, we add the verfication part in the DQC setting (on the mapped data qubits, which are afain in normal ordering, such that we don't
    # need to consider inv_map for indices. However, we do need )
    # determine mapping of ancillas

    # add idling during verification as well

    # this all before uncomputing, then uncompute and do the noisefree stabiliser gadget


    # pass data_circuit,quantum_clifford_verification_circ,ancilla_data_interactions to exectuable, where data qubits will be executed, then 
    # verificaion will be appended baed on the ancilla_data_interactions dict, then unmapping

    #println("\n\n\nqasm: $qasm_prog\n")
    println("\nData Circuit: $data_circuit \n")
    println("Ver. Circuit:$quantum_clifford_verification_circ \n")
    println("Ancilla interactis: $ancilla_data_interactions")
    #ancilla_cores = 

    # create ancilla mapping list

    # TODO: For the flags that were not part of any interaction, discard them and reduce the number of flags accordinglt!
    # NOOOOT NECESSART: Since the flag qubits are appended in the end anyways
    # Also reduce the number of num_ancillas correspodnllgy

    ancilla_map = Vector{Int}()# DONT initialie with fixed number since num_flags is an upper boundundef, num_z_anc + num_x_anc + num_flags)
    #println(num_z_anc, num_x_anc , num_flags)
    anc_order = Dict("z_anc" => 1, "x_anc" => 2, "flag" => 3)
    sorted_dict = sort(collect(ancilla_data_interactions); by = x -> (anc_order[x[1][1]], x[1][2]))
    
    for (index, (ancilla, interactions)) in enumerate(sorted_dict)
        println("Index: $index, key $ancilla, value $interactions")
        #ancilla_map[index] = mode(interactions)
        push!(ancilla_map, mode(interactions))
    end

    num_ancillas = length(ancilla_map) # collects all ancillas that have been used in some interaction
    print("num_ancillas is $num_ancillas, where there are $num_z_anc z ancillas and $num_x_anc x ancillas")
    # for anc_type in ("z_anc", "x_anc", "flag")
    #     for index in sort(collect(keys(ancilla_data_interactions)))
    # for (index, dict_entry) in enumerate(sorted(ancilla_data_interactions)) 
    #     interactions = dict_entry[2](anc_type, index)  
    #     println((anc_type, index), interactions)
    #     ancilla_map ancilla_data_interactions[(anc_type, index)] = mode( interactions )
    # end

    println("Ancilla map: $ancilla_map")
    p = 1e-4
    num_samples = 1e6
    #noise_model = NoiseSpecs(num_samples,p,p,p,p,p,p,p,p,p,p,p,p)
    #noise_model = NoiseSpecs(num_samples,1e-4,1e-4,1e-3,1e-4,1e-3,1e-3,1e-3,1e-1,1e-1,1e-4,1e-2,0)
    tele_p = 1e-2
    noise_model = NoiseSpecs(num_samples,p,p,tele_p,p,p,p,p,tele_p,tele_p,tele_p,tele_p,0)
    # init_noise::Float64               # Initialisation noise
    # idle_depolarising_noise::Float64  # idling depolarising probability
    # idle_depolarising_noise_tele::Float64 # idle depolarising probability under telegate
    # single_q_gate_noise::Float64      # single qubit gate noise probability
    # two_q_gate_noise::Float64         # two-qubit gate noise probability
    # measurement_noise::Float64        # Measurement noise
    # two_q_gate_noise_diff_species::Float64     # two-qubit gate noise probability between communication and memory qubit
    # comm_qubit_init_noise::Float64             # Communication qubit init noise, de facto two qubit depolarising noise to mimic the imperfect creation of Bell pairs
    # comm_idle_depolarising_noise::Float64      # Communication qubit idling depolarising probability
    # single_comm_q_gate_noise::Float64          # Communication qubit single gate depolarising probability
    # comm_qubit_measurement_noise::Float64      # Communication qubit measurement depolarising probability
    # classical_comm_noise::Float64              # Classical communication error

    dqc_state_prep(data_circuit, quantum_clifford_verification_circ, num_ancillas, ancilla_map, code_params, network_specs, noise_model)
    return 42
end



export run_dqc_state_prep

function noise_run()
    p = 1e-3
    noise_model = NoiseSpecs(p,p,p,p,p,p,p,p,p,p,p,p)
    num_samples = 1e5


    mem_errors = 0.001:0.0005:0.01
    #need proper conversion between real T1 and T2 times to noise probs in order to gauge where we are currently at
    #need more elaborate sweep of noise parameter array
    codes = [Shor9()]
    results = zeros(length(mem_errors), 2)


    for (ic, c) in pairs(codes)
        for (i,m) in pairs(mem_errors)
            setup = NaiveSyndromeECCSetup(m)
            decoder = TableDecoder(c)
            #call the dqc_state_prep_sim file
            r = evaluate_decoder(decoder, setup, 10000)
            results[ic,i,:] .= r
        end
    end

    make_decoder_figure(mem_errors, results, "Shor's code with a lookup table decoder")
end

"""
Plot logical vs. physical error rate
-> Level of dependency on differnet kind of noise (e.g., all telegate noises p* and all regular p, see how it behaves)
Threshold would be where DQC logical state fidelity is larger than measurement of two qubits (prep of physical two qubit zero state under measuremtn noise)?
"""

using CairoMakie

function make_decoder_figure(phys_errors, results, title="")
    minlim = min(minimum(phys_errors),minimum(results[results.!=0]))
    maxlim = min(1, max(maximum(phys_errors),maximum(results[results.!=0])))

    fresults = copy(results)
    fresults[results.==0] .= NaN

    f = Figure()
    a = Axis(f[1,1],
        xscale=log10, yscale=log10,
        limits=(minlim,maxlim,minlim,maxlim),
        aspect=DataAspect(),
        xlabel="physical error rate",
        ylabel="logical error rate",
        title=title)
    lines!(a, [minlim,maxlim],[minlim,maxlim], color=:black)
    for (i,sresults) in enumerate(eachslice(fresults, dims=1))
        scatter!(a, phys_errors, sresults[:,1], marker=:+, color=Cycled(i))
        scatter!(a, phys_errors, sresults[:,2], marker=:x, color=Cycled(i))
    end
    f
end


function plot_sweep(parameter_values, state_fidelities, params)
    fig = Figure(resolution = (700, 500))
    ax = Axis(
        fig[1, 1],
        xlabel = "Time (in seconds)",
        ylabel = "Final Steane-7 fidelity",
        title = "Simulation fidelity vs depolarising char. times")

    lines!(
        ax,
        parameter_values,
        state_fidelities,
        linewidth = 3)
    
    vlines!(
    ax,
    0.15, # currently usual execution time of entire circuit TODO: should be some meaningful comparison params.classical_communication_time;            
    linestyle = :dash,
    linewidth = 2,
    label = "")

    scatter!(
        ax,
        parameter_values,
        state_fidelities,
        markersize = 8)
    
    save("src/plots/fidelity_vs_depolarising_noise_time.png", fig)

end

function run_parameter_sweep()
    
    params = retrieve_parameters(1) #  TODO: Refactor this; dummy parms vector to enable the creation of the register lookup and circuit (both depend only on register_size!)
    register_lookup_array, register_start_indices = create_lookup_array(params.register_sizes)      # create lookup array
    circuit = steane_encoding_circuit(params)                # build initial circuit
    #TODO: block all communication qubit layers! Can be done via row check != comm_qubits,
    #TODO: Include check for no overlaps within one layer
    
    depolarising_times = collect(1e2:1e4:1e3) # start:size:stop
    state_fidelities = Float64[]

    #TODO: Add standard deviation 
    #TODO: sweep over other parameters as well
    for depolarising_time in depolarising_times
        params = retrieve_parameters(depolarising_time)  # retrieve parameters, changing bell state fidelity in every sweep
        fid_depol = 0
        num_runs = 1
        for _ in 1:num_runs
            sim_fid = run_dtsimulation(params, circuit, register_lookup_array, register_start_indices)
            fid_depol += sim_fid.fidelity
        end
        
        fid_depol = fid_depol/num_runs
        push!(state_fidelities, fid_depol)
        @info "For depolarising tau $depolarising_time, we obtain final state fidelity $(fid_depol)"
    end

    # Profiling
    #@btime run_dtsimulation($params, $circuit, $register_lookup_array, $register_start_indices) # running profiling with last set of params
    #@timed run_dtsimulation(params, circuit, register_lookup_array, register_start_indices)
    #@profile run_dtsimulation(params, circuit, register_lookup_array, register_start_indices)
    #@profile for _ in 1:100
    #    run_dtsimulation(params, circuit, register_lookup_array, register_start_indices)
    #end
    #Profile.print(format=:flat, sortedby=:count)
    #plot_sweep(depolarising_times, state_fidelities, params) # !! takes the last params iteration
end

end