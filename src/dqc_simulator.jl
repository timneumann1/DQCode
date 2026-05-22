"""
Adapted from https://github.com/QuantumSavory/QuantumClifford.jl/blob/master/src/ecc/circuits.jl
and https://github.com/QuantumSavory/QuantumClifford.jl/blob/2664eba07a461441ea051a17cad7c725f9576176/src/ecc/decoder_pipeline.jl
build layers adapted from MQT QECC
"""

module DQCodeSimulator

using ..Types
using ..Helper

using QuantumClifford
using QuantumClifford: MixedDestabilizer, sHadamard, sCNOT, @S_str, Register, continue_stat, AbstractSingleQubitOperator, sMZ, AbstractTwoQubitOperator, AbstractOperation, AbstractStabilizer, AbstractNoise, apply_single_x!, apply_single_y!, apply_single_z! #, sX, SY, SZ
using QECCore
using QuantumClifford.ECC: DecoderCorrectionGate, CSSTableDecoder, decode,  AbstractSyndromeDecoder, faults_matrix, ClassicalTableDecoder, create_lookup_table, parity_checks, decode
using Combinatorics: combinations
using ProgressMeter
using PyCall
using StatsBase


import QuantumClifford: apply!, affectedqubits, applynoise! # we want to extend this with ConditionalGate

qiskit = pyimport("qiskit")
qasm2 = pyimport("qiskit.qasm2")

export dqc_ft_encoding_simulation, dqc_non_ft_encoding_simulation
export construct_DQC_executable_circuit
export dqc_logical_evaluation


function dqc_non_ft_encoding_simulation(code_params::CodeParameters, network_specs::NetworkSpecifications, circuit::Vector{AbstractOperation})

    data_circuit = deepcopy(circuit)
    
    noise_verif = NoiseSpecs(1e5,0,0)

    quantum_clifford_verification_circ = Vector{AbstractOperation}() # no verification circuit ⇔ non-FT
    DQC_circuit = construct_DQC_executable_circuit(data_circuit, quantum_clifford_verification_circ, 0, [], network_specs, noise_verif, 0.0)

    initial_state = Register(one(MixedDestabilizer, network_specs.num_data_and_comm_qubits+num_ancillas),network_specs.num_comm_qubits + num_ancillas)
    state, stat = mctrajectory!(initial_state, vcat(DQC_circuit, VerifyOp(code_params.target_state, network_specs.data_qubits))) #, trajectories=noise_verif.n_samples)
    @assert stat == true_success_stat "Adding the verifcation circuit compromised the data circuit"

    # Also verify the dqc_evaluation, which should give an error rate of zero in the noiseless setting
    logical_error_rate, acceptance_ratio = dqc_logical_evaluation(data_circuit, quantum_clifford_verification_circ, 0, [], code_params, network_specs, noise_verif, 0.0)
    @assert logical_error_rate == 0.0 "Logical error rate is non-zero in noiseless setting; please check your circuit"  
    @assert acceptance_ratio == 1.0 "Not all runs accepted in noiseless setting; please check your verification circuit"
    @info "Noiseless testing of FT encoding circuit successful."

    num_samples = 5e5
    ps = 10 .^ range(log10(5e-5),-3,length=30) # 5e-5:2e-5:1e-3  # 50
    p_bells = 10 .^ range(-3,log10(5e-2), length = 30) #1e-3:2e-3:5e-2  #25
    telegate_idle_depth = 8

    data = NamedTuple[]
    progress = Progress(length(ps) * length(p_bells); desc="(p, p_bell) sweep", dt=1)

    for (p, p_bell) in Iterators.product(ps, p_bells)
        #  for each combination, initialise NoiseSpecs(...) according to above considerations
        noise_model = NoiseSpecs(num_samples, p,p_bell)
        p_idle_telegate_layer = 1-(1-noise_model.p)^telegate_idle_depth # compute the telegate_layer idle error prob from the current p
        logical_error_rate, acceptance_ratio = dqc_logical_evaluation(data_circuit, quantum_clifford_verification_circ, num_ancillas, ancilla_map, code_params, network_specs, noise_model, p_idle_telegate_layer)
        push!(data, (p=p, p_bell=p_bell, logical_error_rate=logical_error_rate, acceptance_ratio=acceptance_ratio))
        next!(progress; showvalues=[(:p, p), (:p_bell, p_bell), (:logical_error_rate, logical_error_rate), (:acceptance_ratio,acceptance_ratio)])
    end
   
    return data, data_circuit, DQC_circuit
end


function dqc_ft_encoding_simulation(num_samples, ps, p_bells, telegate_idle_depth, code_params::CodeParameters, network_specs::NetworkSpecifications, mqt_path::String, circuit::Vector{AbstractOperation}, method::String)
    
    data_circuit = deepcopy(circuit)
    #println("Initial circuit: $circuit")

    qasm = qc_circuit_to_qasm(circuit)
    
    #println("QASM version: $qasm")

    # Procedure: We pass a qasm string with the optimised encoding circuit, and get back a qasm string with the verification 
    #(in MQT QECC, the result was a qiskit circuit, which was then converted to qasm in order to make the output stream usable, and then converted to qiskit here again)
    
    # everthing is now collapsed in the Z check matrix (Z stabilisers and logical Zs), since they are equaivlent as stabilisers of the desired state
    # and we can measrue both to obtain the correct one (of course cofrectinng for hook errors); our initial definition of logical still applies and
    # can be used for the noiseless procedure

    python_bin = joinpath(mqt_path, ".venv/bin/python3")
    script_path =  joinpath(mqt_path, "scripts/verification_circuit.py")

    @info "Retrieving verification circuit"

    verification_circ_qasm = readchomp(`$(python_bin) $(script_path) $qasm $(code_params.distance) $(method)`)
    verification_circ = qasm2.loads(verification_circ_qasm) 
    @info "Retrieved verification circuit"
    println("Verification Cirucit: $verification_circ_qasm")
    #verification_circ.draw(output="mpl", initial_state=true, fold=-1, scale=0.4)
    #plt.show()

    # qasm program could be reordered but should still be intact:
    #println("Final QASM Program orig + verification: $qasm_prog")

    quantum_clifford_verification_circ = Vector{AbstractOperation}()
    # We actually parse throught the qiskit circuit now, which is more convenient than the QASM string

    ancilla_data_interactions = Dict{Tuple{String, Int}, Vector{Int}}() # will contain ancilla qubits, e.g., Z_anc, 1 as keys, and interacting data qubits as values
   
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
    
    for instruction in verification_circ.data
        gate = instruction.operation.name
        
        if gate == "h"
            qubits = instruction.qubits
            bit_info = verification_circ.find_bit(qubits[1])
            reg_name = String(bit_info.registers[1][1].name)
            if reg_name == "q"
                continue
            end
            index = Int(bit_info.index)
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
                #println("Target: $target, network_specs.num_comm_qubits: $(network_specs.num_comm_qubits)")
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
                    
            end

        elseif gate=="measure" 
            qubits = instruction.qubits
            bit_info = verification_circ.find_bit(qubits[1])
            reg_name = String(bit_info.registers[1][1].name)
            @assert reg_name != "q" # only ancillas will be measured
            index = Int(bit_info.index)
            push!(quantum_clifford_verification_circ, sMZ(index+network_specs.num_comm_qubits+1, index-network_specs.num_data_qubits+network_specs.num_comm_qubits+1)) 
        else
            throw("Verification circuit contains gates that are currently not supported.")
        end
    end
    
    # determine the best placement for the ancilla qubits, noting which cnots have to be applied to them, and noting the index of the measueretn (on which we later base the discarding)
    # choose the core with the max occurences fr the ancilla to be placed and store it in core_mapping

    # println("\nData Circuit: $data_circuit \n")
    # println("Ver. Circuit:$quantum_clifford_verification_circ \n")
    # println("Ancilla interactis: $ancilla_data_interactions")
   
    ancilla_map = Vector{Int}() # don't initialise with fixed number, since num_flags is an upper bound
    anc_order = Dict("z_anc" => 1, "x_anc" => 2, "flag" => 3)
    sorted_dict = sort(collect(ancilla_data_interactions); by = x -> (anc_order[x[1][1]], x[1][2]))
    
    for (index, (ancilla, interactions)) in enumerate(sorted_dict)
        #println("Index: $index, key $ancilla, value $interactions")
        push!(ancilla_map, mode(interactions))
    end

    num_ancillas = length(ancilla_map) # collects all ancillas that have been used in some interaction
    print("Verification circuit uses $num_ancillas ancillas, with $num_z_anc z-ancillas, $num_x_anc x-ancillas, and $(num_ancillas-num_z_anc-num_x_anc) flag qubits")
    
    #println("Ancilla map: $ancilla_map")
    @info "Data Circuit: $data_circuit"
    
    @info "Verification Circuit: $quantum_clifford_verification_circ"
    # Now, we have retrieved the verification circuit. Next, we want to execute it for multiple noises, which is handled by dqc_state_prep, which calls make_executable and then executes

    # TODO: save data circuit and verification circuit in noise folder
    # TODO: Add noise sweep here

    # this is the most important object

  

    #evaluate encoding of Z state by identifying decoding correction and checking via measured fault in Z logical classically (we check whether a logical X error has been applied, leading to a logical Z error in the stabiliser)

    # what would the logical Z measuremten show if the state was + and thus stabilised by X (what else than 0 for 0 state and 1 for 1 state)?
    # It would show a probabilistic 0 and 1. Hence, we can test whether our circuit really works by setting all noise to zero and measuring success.
    # If we actually encoded the plus state, the result would only be true in 50% of the cases; we run this verification ccheck ones before initialising the noisy simulaton

    @info "Noiseless testing of FT encoding circuit ..."

    noise_verif = NoiseSpecs(1e5,0,0)
    DQC_circuit = construct_DQC_executable_circuit(data_circuit, quantum_clifford_verification_circ, num_ancillas, ancilla_map, network_specs, noise_verif, 0.0)

    initial_state = Register(one(MixedDestabilizer, network_specs.num_data_and_comm_qubits+num_ancillas),network_specs.num_comm_qubits + num_ancillas)
    state, stat = mctrajectory!(initial_state, vcat(DQC_circuit, VerifyOp(code_params.target_state, network_specs.data_qubits))) #, trajectories=noise_verif.n_samples)
    @assert stat == true_success_stat "Adding the verifcation circuit compromised the data circuit"

    # Also verify the dqc_evaluation, which should give an error rate of zero in the noiseless setting
    logical_error_rate, acceptance_ratio = dqc_logical_evaluation(data_circuit, quantum_clifford_verification_circ, num_ancillas, ancilla_map, code_params, network_specs, noise_verif, 0.0)
    @assert logical_error_rate == 0.0 "Logical error rate is non-zero in noiseless setting; please check your circuit"  
    @assert acceptance_ratio == 1.0 "Not all runs accepted in noiseless setting; please check your verification circuit"
    @info "Noiseless testing of FT encoding circuit successful."

    # We now define noise regimes, motivated by literature on ion traps with photonic interconnects.
    # Here, we assume data/memory/ancilla qubits, which make up the code (and the ancillas for verification), and networking/communication qubits, which serve to connect modules via photonic interconnects.
    # The precise process of how photons, for example, are measured in Beam Splitters to entangle the communication qubits and thereafter perform the EJPP protocol,
    # are not of relevance here. We note however, that this process creates unproportionally much noise, and takes a long time, which is why we introduce the comm_init_noise noise,
    # defining the error probability of Bell pair entanglement of the communication qubits via photonic interconnects. This fidelity is usually well above 90% (Main: 97%, https://arxiv.org/pdf/1911.10841: 94%).
    # Besides from this, there are local operations, such as single- and two-qubit gate noise, measurement noise and initalisation error of data qubits. 
    #Quantinuum Helios noise specs (typical): 1q fidelity ~ 3e-5, 2q fidelity 8e-4, state prep + msm 5e-4, memory error depth 1 6e-4, IonQ 2q with fidelity 99.99%, i.e. 1e-4 fidelity
    #For simplicity, we assume all these to be comparable in magnitude, where it is reasonable to sweep a range from 5e-5 to 1e-3 in increments of 2.5e-5 (50 evals). (By the way that we define the 2-qubit depolarising
    # channel, the expected number of noise applications is actually 1.6x larger for 2-qubit gates than for 1-qubit gates, which benefits the distinction of the two). The communication qubits do not experience much decoherence on top of the initialisation error in our setup, which is why we exclude it from the simulation. Furthermore, we
    # assume that single-qubit gate and measurement noise are comparable between memory and communicatioin qubits. Lastly, we have two-qubit gate noise between ions of different species, which 
    # can be assumed to be roughly 98% (Main), but which we will subsume under p for simplicity (in the end, this is a local operation), and due to the probabilistic Bell state creation delay, we need to account for additional decoherence of the memory qubits. Assuming a two-qubit gate takes on the order of hundredes of μs incl. cooling, whereas Bell state creation takes low ms regime ( https://arxiv.org/pdf/1911.10841)
    # we are probably dealing with a factor of d ~ 1e1 to 1e3, which we account for by increasing the noise probability of the idling depolarising channel according to 1-(1-p)^d, where p is memory_idle_depolarising_noise, the usual 1 layer
    # memory noise, which is on the order of 1e-4 to 1e-3 (according to Quantinuum) (and will be swept in roughly this regime). However, as mentioned in Main, for example, we can use DD to mitigate this dephasing effect,
    # which justigies choosing d on the lower end of the spectrum. In particular, we adapt Floquette's claim of 5 gate cycles for Bell state creation, and add three more for the EJPP protocol two-qubit gate, measurement and correction gates, giving an added noise corresponding to d=8 layers (this might be realistic for near-term experiments).
    # If a layer is a telegate layer, we will thus apply the corresponding p_idle_telegate_layer channel to all memory/ancilla qubits, whereas for other layers, we only assume decoherence for passive qubits.

    # In summary, we have local noise (memory init, single and two qubit noise on data/comm, measurement noise on data/comm), which we denote with p, and additionally
    # idling noise during telegate operations, which is a function of p (see above), two-qubit gate noise between species, which we denote p_mixed ≈1e-2, and most importantly,
    # the communication qubit Bell state initialisation error probabiliy p_bell.
    # As seen above, it is reasonable to sweep 5e-5:2e-4:1e-3 for p, and to sweep 1e-3 to 5e-2 for p_bell. (this way, p_bell is always larger than p, and comparable with p_mixed)
    # It is reasonable to assume that classical communication noise is negligible.

    # We then hope to extract relations between p and p_bell at various pseudothresholds (similar to Floquet), in order to gauge how much noise a QEC system can tolerate
    # The pseudothreshold lies at p = p_logical, since initialisation \textit{any} bare qubit has error probability p, whereas initialising \textit{any} logical qubit has max error rate logical_rate



    # num_samples = 1e3
    # # Evaluating 50x25 = 1250 logical error rates, with 1e6 samples each (but early stopping condition)
    # ps = 10 .^ range(-2,-1,length=2) # 5e-5:2e-5:1e-3  # 50
    # #p_idle_telegate_layer = 1-(1-p)^d
    # p_bells = 10 .^ range(-2,-1, length = 2) #1e-3:2e-3:5e-2  #25
    # telegate_idle_depth = 1


    data = NamedTuple[]
    progress = Progress(length(ps) * length(p_bells); desc="(p, p_bell) sweep", dt=1)

    for (p, p_bell) in Iterators.product(ps, p_bells)
        #  for each combination, initialise NoiseSpecs(...) according to above considerations
        noise_model = NoiseSpecs(num_samples, p,p_bell)
        p_idle_telegate_layer = 1-(1-noise_model.p)^telegate_idle_depth # compute the telegate_layer idle error prob from the current p
        logical_error_rate, acceptance_ratio = dqc_logical_evaluation(data_circuit, quantum_clifford_verification_circ, num_ancillas, ancilla_map, code_params, network_specs, noise_model, p_idle_telegate_layer)
        push!(data, (p=p, p_bell=p_bell, logical_error_rate=logical_error_rate, acceptance_ratio=acceptance_ratio))
        next!(progress; showvalues=[(:p, p), (:p_bell, p_bell), (:logical_error_rate, logical_error_rate), (:acceptance_ratio,acceptance_ratio)])
    end
   
    return data, data_circuit, quantum_clifford_verification_circ, DQC_circuit, num_ancillas, num_x_anc, num_z_anc, ancilla_map
end


# ---------------------------------------------------------------
# ------------ DQC Logical Evaluation ---------------------------
# ---------------------------------------------------------------


function dqc_logical_evaluation(data_circuit, verification_circuit, num_ancillas, ancilla_map, code_params, network_specs, noise, p_idle_telegate_layer::Float64)

    DQC_circuit = construct_DQC_executable_circuit(data_circuit, verification_circuit, num_ancillas, ancilla_map, network_specs, noise, p_idle_telegate_layer)

    # ----- For the noiseless decoding ------
    css_lut_decoder = CSSTableDecoder(code_params.qec_code, error_weight =  Int(floor((code_params.distance-1)/2)) )
    #println("X Lookup Table: $(css_lut_decoder.tabledecoderx.lookup_table)")
    #println("Z Lookup Table: $(css_lut_decoder.tabledecoderz.lookup_table)")
    H = parity_checks(css_lut_decoder)
    # fault matrix is a (2k)x(2n) dimensional matrix, and to determine the logical Z part, we need the last k rows: O[end÷2+1:end,:]
    faults_matrix_z = css_lut_decoder.faults_matrix[end÷2+1:end,:] #  decoder has the faults_matrix as attribute

    #println("PARITY CHECKS: $H")
    # Build perfect syndrome circuit, with ancillas being appended to DQC cores (consisting of data qubits) and comm_qubits register
    # and bit string measurements being written into classical register n.num_registers*(n.num_registers-1) + #verification classical registers + 1
    noisefree_syndrome_circ, num_noisefree_syndrome_ancillas, noisefree_syndrome_bits = syndrome_circuit(H, network_specs.num_data_and_comm_qubits + num_ancillas + 1, network_specs.num_comm_qubits + num_ancillas +1, network_specs)

    noisefree_logical_Z_circ, num_noisefree_logical_Z_ancillas, noisefree_logical_Z_bits = syndrome_circuit(code_params.logical_Zs, network_specs.num_data_and_comm_qubits + num_ancillas + num_noisefree_syndrome_ancillas + 1, last(noisefree_syndrome_bits)+1, network_specs )
    # prior: QECTools.
    total_number_qubits = network_specs.num_data_and_comm_qubits + num_ancillas +  num_noisefree_syndrome_ancillas + num_noisefree_logical_Z_ancillas
    total_number_classical_regs = network_specs.num_comm_qubits + num_ancillas + num_noisefree_syndrome_ancillas + num_noisefree_logical_Z_ancillas

    # HERE, the verificaion circuit is already part of the DQC circuit
    circuit = vcat(DQC_circuit, noisefree_syndrome_circ, noisefree_logical_Z_circ ) #sX(9),sX(13),sX(14),sX(19), sX(21) gives a logical X error, [sHadamard(1),sHadamard(2), sHadamard(3), sHadamard(7), sHadamard(8), sHadamard(9), sHadamard(13),sHadamard(14), sHadamard(15), sHadamard(19), sHadamard(20), sHadamard(21)] a logical H
    #logical_failures_pre_decoding = 0

    logical_failures = zeros(code_params.k) # we determine  the logical Z failures per qubit, then take the max over all > later compare to initialisation failure prop of a single physical qubit
    discarded_runs = 0
    #apply_correction = false # for code testing
    #println("\n\n\n $circuit \n\n\n")

    n_samples = noise.n_samples
    initial_state = Register(one(MixedDestabilizer,total_number_qubits),total_number_classical_regs)
    for sample in 1:n_samples

        state, stats = mctrajectory!(copy(initial_state), circuit)
        
        # HERE, determine whether or not to discard the state based on the ancilla bits being triggered or not

        verification_bits = Vector{Bool}(state.bits[network_specs.num_comm_qubits+1:network_specs.num_comm_qubits+num_ancillas]) # num_ancillas contains z_anc, x_anc and flags

        if any(verification_bits)==1
            #@info "Discarded"
            discarded_runs +=1
            continue
        end

        syndrome = Vector{Bool}(state.bits[noisefree_syndrome_bits])
        
        # if any(measured_logical_Z_bits)==1 || any(syndrome)==1 
        #     logical_failures_pre_decoding +=1
        # end

       
        error_guess = decode(css_lut_decoder, syndrome)
        # error_guess collects n guesses for x errors, then n guess for z errors, where we are interested in the first n guesses (whether the decoder predicts that a certain X error has happened)
        #@info "ERROR GUESS: $error_guess"

        #@info "Syndrome is: $syndrome"
        #@info "Error guess: $error_guess"
        

        if isnothing(error_guess)
            # Table decoder can't find a matching syndrome bc error weight is too high
            # this edge case only occurs if the syndrome is non-trivial, hence it is safe to say that if it occurs, there will be a logical error (since the state is incorrect yet nothing got corrected)
            # Since we are counting logical errors per logical qubit, we thus manually set the error guess to be all-zero, such that the faulty logical qubits
            # accumulate the error from non-decoding via detection of a logical error in the next loop
            #error_guess = zeros(code_params.n) # give n zero-guesses for X errors
            logical_failures .+= 1
            continue
        end

        measured_logical_Z_bits = state.bits[noisefree_logical_Z_bits] # since we encode the logical zero state, the bit value is effectively a fault value: 0 means logical Z, 1 means logical -Z, which is a fault in the respective qubit
        
        k = size(faults_matrix_z, 1) 
        n = size(faults_matrix_z, 2)
        @assert k == code_params.k
        @assert n == 2*code_params.n


        #no_Z_error = true

        for j in 1:k # iterate over the k logical Z operators
            sum_mod = 0
            @inbounds @simd for q in 1:n # iterate over all the physical qubits (/error locations)
                sum_mod += faults_matrix_z[j, q] * error_guess[q] 
            end
            sum_mod %= 2  # will be 1 if there is an odd number of agreed indiced between fault matrix and error_guess for the given jth logical Z operator
            # the logic behind this is error degeneracy of the code (since we are able to decode, we will be back in the codespace again, so now we check if we are in the correct state within the codespace):
            # If there are is an even number of corrections in error guess coinciding with locations that actually lead to
            # a logical flip of this operator, then applying the even number cancels out the flip. Then if there was a logical error for this logical operator, we will not correct it (if there was not)
            # i.e., sum_mod == measured_fault, we are fine since the even number of corrections fixes the state.
            # If there is an odd number of corrections matching, we perform a logical correction. Now if there was an error (agian sum_mod == mneasured fault), this is exactly what we want. If not, we INTRODUCE an error, 
            # again leading to a logical failure. 
            if sum_mod != measured_logical_Z_bits[j] # measured_logical_Z_bits = measured_Z_faults
                logical_failures[j] += 1
                #no_Z_error = false
            end
        end


        # Alternatively, apply correction and see if we ended up in the correct state:
        # When checking analyitcally for logical errors with the error_guess, we can also use the error guess to apply a correction which maps back to codespace, 
        # then the verifyOp should give exactly the same outcome as the Z logical analytical computation

        # Here, we apply the correction gatehttps://github.com/QuantumSavory/QuantumClifford.jl/blob/444f341a50d2926b16b63d98586b8b06a7b6ac10/src/ecc/decoder_correction_gate.jl,
        # whihc maps back to the codespace  (since there is a correctable syndrome, this measn that the correction will bring us back to the codespace) so stabilisers are satisfied. 
        
        # Then we check whether the state overlaps fully with the target zero state, and whether this info matches with the analytical computation above.


        # correction_gate = DecoderCorrectionGate(css_lut_decoder, network_specs.data_qubits, noisefree_syndrome_bits )
        # final_state,_ = mctrajectory!(state,[correction_gate])
        # no_Z_error2 = compare_states(final_state, code_params.target_state, network_specs)
        # if !no_Z_error
        #     println("Z error registered")
        # end
        # @assert no_Z_error == no_Z_error2 "The analytic computation of logical errors did not yield the same result as the state tomography following the application of a correction gate."
        

        if maximum(logical_failures) >= 500# n_samples*noise.ps[floor(end:2)] 

            # For a given p, under linear scaling (~upper bound) we roughly expect n_samples*p errors. Thus, we take a mid-range p and limit the execution for all larger p, once this
            # number of logical failures has been recordedFor the upper half of physical noise values, we 
            # In the given noise regime, we expect logical error rates on the order of 1e-4 to 5e-3, which corresponds to 50 to 1000 logical errors for the noisiest Z logical observable over 5e5 executions.
            # Thus, collecting 200 errors (MQT QECC:500) at most should give us a good estimate of the true error rate.
            # This helps to decreae overall runtime, since it reduces the execution time particularly in the high-noise regime, in which we are less interested.
            @info "For physical error rate $(noise.p) and Bell state error rate $(noise.p_bell), we collected $(maximum(logical_failures)) logical failures after $sample iterations."
            n_samples = sample
            break
        end

    
        
        #print(correction_gate)
        #correction_gate = DecoderCorrectionGate(decoder, 1:7, 1:6)
        #apply!(circuit, correction_gate)
        #state_decoded, stats_decoded = mctrajectory!(Register(state.stab, network_specs.num_registers*(network_specs.num_registers-1)+length(syndrome_bits)+length(logical_Z_bits)), [correction_gate, syndrome_circ, logical_Z_circ ])

        # initial_state = Register(one(MixedDestabilizer,total_number_qubits),network_specs.num_registers*(network_specs.num_registers-1)+length(syndrome_bits)+length(logical_Z_bits))
        # st, stat = mctrajectory!(initial_state, copy(circuit))
        
        # syndrome = state_decoded.bits[syndrome_bits]
        # measured_logical_Z_bits = state_decoded.bits[logical_Z_bits] 

        # println(st.bits)
        #println("Syndrome: $syndrome")
        #println("Logical Zs: $measured_logical_Z_bits")
    end

    #println("DQC circuit length: $(length(filter!(g -> g isa NoiseOp, DQC_circuit)))")
    acceptance_ratio = 1-discarded_runs/n_samples
    logical_error_rate = mean(logical_failures)/(n_samples-discarded_runs)
    #println("Without decoding, the logical error rate is $(logical_failures_pre_decoding/noise.n_samples) ")
    #println("Over $(n_samples) runs, there were $discarded_runs discarded runs -> acceptance ratio: $acceptance_ratio, for the kept runs the the logical error rate (after decoding) is $logical_error_rate")
    return logical_error_rate, acceptance_ratio
end


# ----------- Helper functions for DQC logical evaluation --------------

function perfect_ancillary_paulimeasurement(p::PauliOperator, ancillary_index, bit_index, network_specs)
    circuit = AbstractOperation[]
    num_data_qubits = nqubits(p)
    @assert num_data_qubits == network_specs.num_data_qubits
    for qubit in 1:num_data_qubits
        # for the perfect ancillary measurement, qubits are back in their correct position already
        if p[qubit] == (1,0)
            push!(circuit, sXCX(qubit, ancillary_index)) # X-controlled X     
        elseif p[qubit] == (0,1)
            push!(circuit, sCNOT(qubit, ancillary_index)) # Z-controlled X
        elseif p[qubit] == (1,1)
            push!(circuit, sYCX(qubit, ancillary_index)) # Y-controlled X
        end
    end
    p.phase[] == 0 || push!(circuit, sX(ancillary_index))
    mz = sMRZ(ancillary_index, bit_index) #sMRZ may be used here for convenience
    push!(circuit, mz)

    return circuit
end

function syndrome_circuit(parity_check_tableau, ancillary_index, bit_index, network_specs)
    syndrome_circ = AbstractOperation[]
    ancillaries = 0
    bits = 0
    for check in parity_check_tableau
        append!(syndrome_circ, perfect_ancillary_paulimeasurement(check, ancillary_index+ancillaries, bit_index+bits, network_specs))
        ancillaries +=1
        bits +=1
    end

    #print("We consumed $bits bits for the nosiefre ancilla")

    return syndrome_circ, ancillaries, bit_index:bit_index+bits-1
end



# ---------------------------------------------------------------
# ------------ Construction of executable circuit ---------------
# ---------------------------------------------------------------

function construct_DQC_executable_circuit(data_circuit, verification_circuit, num_ancillas, ancilla_map, n::NetworkSpecifications, noise::NoiseSpecs, p_idle_telegate_layer::Float64)
    
    # Map, then add verification on mapped, then unmap

    circuit = Vector{QuantumClifford.AbstractOperation}()
   
    # in the permutation, [1,9,...] indicates that the 9th element gets permuted to second position, "9 is mapped to 2"
    # Apply the inverse permutation of the mapping by applying transpoistions of inverse perm in left action <-> transpoistions of perm in right action via reverse(mapping) [the mapping contains transposition derived from the permutation, implementing it in left action]
    for (i,j) in reverse(n.mapping_transpositions)
        push!(circuit, sSWAP(i,j))  # We could also use comm_perm_idx or comm_inv_perm_idx, since the relabeling based on the permutation conjugtes and thus fixes the permutation induces by the transposition SWAPS
    end
        
    #go through circuit and accumulate ASAP layers (we may assume that no further grouping for telegates can be done)
    layers_enc_circ = build_layers(data_circuit, n.num_data_qubits)

    # Add initialisatin noise to all qubits
    # Depolarising channel: https://github.com/QuantumSavory/QuantumClifford.jl/blob/74ee758e87f5d7b1255d6747b346cff15ee10cea/src/noise.jl#L63-73
    add_noise(circuit, [n.inv_map[data_q] for data_q in collect(1:n.num_data_qubits)], 3/2*noise.p) # initialisation noise p

    for layer in layers_enc_circ
    # for every layer, add gate noise for every normal gate and idling noise for any idle qubits (if there is telegates in the layer, we increase the probability of the noise channel)
    # also, insert the telegate gadgets in place (this includes the folloing noise: comm init noise, two qubit depolarising for each gate , measeuremtn noise, classical noise, one more gate noise for single quits )
    
        #find idling gates and apply idle_depolarising_noise
        affected_qubits = Set( Iterators.flatten( [affectedqubits(gate) for gate in layer] ) ) 
        idle_qubits = setdiff(1:n.num_data_qubits, affected_qubits)
        idle_qubits_DQC = [n.inv_map[idle_q] for idle_q in idle_qubits]
    
        telegates_layer = false
        tele_data_qubits = Vector{Int64}()

        for gate in layer
            T = typeof(gate)
            if T <: AbstractSingleQubitOperator
                qubit = affectedqubits(gate)[1]
                DQC_qubit = n.inv_map[qubit]
                push!(circuit, sHadamard(DQC_qubit))
                add_noise(circuit, [DQC_qubit], noise.p) # Single-qubit noise
                
            elseif T <: AbstractTwoQubitOperator
                control = affectedqubits(gate)[1]
                target = affectedqubits(gate)[2]
                DQC_control = n.inv_map[control]
                DQC_target = n.inv_map[target]
                control_register = n.register_lookup_array[DQC_control] 
                target_register = n.register_lookup_array[DQC_target] 
                

                if control_register == target_register 
                    push!(circuit, sCNOT(DQC_control, DQC_target))
                    add_noise(circuit, [DQC_control, DQC_target], noise.p; two_qubits = true) # two-qubit noise
                else
                    # Perform telegate between control and target qubit in different registers, 
                    circuit = add_telegate(circuit, DQC_control, DQC_target, control_register, target_register, n, noise, p_idle_telegate_layer)
                    telegates_layer = true
                    push!(tele_data_qubits, DQC_control)
                    push!(tele_data_qubits, DQC_target)
                end
            else
                throw("Circuit contains gates that have not been classified as Single- or Two-Qubit gate so far.")
            end
        end
        telegates_layer ? add_noise(circuit, setdiff(1:n.num_data_qubits, Set(tele_data_qubits)), p_idle_telegate_layer) : add_noise(circuit, idle_qubits_DQC, noise.p) # idling noise (if telegate: to all non-telegate qubits, no matter if idle or not (since telegate is much longer), otherwise: idling noise p on idle qubits)

    end

    
    # Add the verification circuit before reversing the virtual mapping
    # telegates between ancillas and data qubits can use the comm qubits of the register, likewise for telegates between ancillas and flags
    layers_ver_circ = build_layers(verification_circuit, n.num_data_and_comm_qubits + num_ancillas)

    #@info "Mapping: $ancilla_map"

    all_qubits = 1:(n.num_data_and_comm_qubits + num_ancillas)
    ancilla_qubits = setdiff(all_qubits, 1:n.num_data_and_comm_qubits)

    #add_noise(circuit, [data_q for data_q in collect(1:n.num_data_qubits)], noise.idle_depolarising_noise) # since noise is applied to all qubits, we don't need to worry about mapping
    add_noise(circuit, [data_q for data_q in ancilla_qubits], 3/2*noise.p) #  init noise on ancilla qubits

    for layer in layers_ver_circ
        # analogous to raw encoding circuit
        affected_qubits = Set( Iterators.flatten( [affectedqubits(gate) for gate in layer] ) ) 
        idle_qubits = setdiff(all_qubits, union(Set(n.comm_qubits), affected_qubits))
        #println("all qubits: $all_qubits, idle: $idle_qubits, data: $(n.data_qubits)")
        idle_qubits_DQC = Vector{Int}()
        for idle_q in idle_qubits
            if idle_q in n.data_qubits
                push!(idle_qubits_DQC, n.inv_map[idle_q] )
            elseif idle_q in ancilla_qubits
                push!(idle_qubits_DQC, idle_q)
            else 
                throw("Communication qubits are not idling since they will be re-initialised before usage")
            end
        end

        telegates_layer = false
        tele_qubits = Vector{Int64}()

        for gate in layer
            T = typeof(gate)
            if T <: AbstractSingleQubitOperator
                qubit = affectedqubits(gate)[1]
                # In the verifciation circuit, Hadamard gates are ONLY applied to ancillas, which sit at their regular index
                push!(circuit, sHadamard(qubit))
                add_noise(circuit, [qubit], noise.p) # single-qubit noise, ancilla qubits experience the same sort of noise, since they are of the same physical type
            elseif T <: AbstractTwoQubitOperator
                control = affectedqubits(gate)[1]
                target = affectedqubits(gate)[2]
                DQC_control = -1
                DQC_target = -1
                control_register = -1
                target_register = -1

                if control in ancilla_qubits 
                    DQC_control = control
                    control_register = ancilla_map[control-n.num_data_and_comm_qubits] 
                else
                    DQC_control = n.inv_map[control]
                    control_register = n.register_lookup_array[DQC_control] 
                end
                if target in ancilla_qubits
                    DQC_target = target
                    target_register = ancilla_map[target-n.num_data_and_comm_qubits] 
                else
                    DQC_target = n.inv_map[target]
                    target_register = n.register_lookup_array[DQC_target] 
                end
                
                @assert DQC_control > 0 
                @assert DQC_target > 0 
                @assert control_register > 0 
                @assert target_register > 0 

                if control_register == target_register 
                    push!(circuit, sCNOT(DQC_control, DQC_target))
                    add_noise(circuit, [DQC_control, DQC_target], noise.p; two_qubits = true) # two-qubit noise
                else
                    circuit = add_telegate(circuit, DQC_control, DQC_target, control_register, target_register, n, noise, p_idle_telegate_layer)
                    telegates_layer = true
                    push!(tele_qubits, DQC_control)
                    push!(tele_qubits, DQC_target)
                end
            elseif T <: sMZ
                add_noise(circuit, affectedqubits(gate), 3/2*noise.p) # measurement noise the ancilla is of the same type, thus we have the same measurement noise
                push!(circuit, gate) # only ancillas are ever measured, so we don't need a remapping
            else
                throw("Circuit contains gates that have not been classified as Single- or Two-Qubit gate so far.")
            end
        end
       # println("idle qubits: $idle_qubits_DQC")
        
        #telegates_layer ? add_noise(circuit, idle_qubits_DQC, noise.p_idle_telegate_layer) : add_noise(circuit, idle_qubits_DQC, noise.p) # idling noise
        telegates_layer ? add_noise(circuit, setdiff(vcat(1:n.num_data_qubits, n.num_data_and_comm_qubits+1:length(all_qubits)), Set(tele_qubits)), p_idle_telegate_layer) : add_noise(circuit, idle_qubits_DQC, noise.p) # idling noise (if telegate: to all non-telegate qubits (data + ancilla), no matter if idle or not (since telegate is much longer), otherwise: idling noise p on idle qubits)

    end

    # Revert swapping for measurement of target state
    for (i,j) in n.mapping_transpositions # reverse-reverse -> revert back to original permutation
        push!(circuit, sSWAP(i, j)) 
    end

    return circuit
end


# ------------ ASAP Layers for construction of executable circuit ---------------

function build_layers(gates,num_qubits)
    # adapted from MQT circuit_utils.py file
    gate_list = copy(gates)
    layers = Vector{Vector{AbstractOperation}}()
    while !isempty(gate_list)
        idx = 1
        layer = Vector{AbstractOperation}()
        qubit_used_in_layer = falses(num_qubits) 
        del_gates = Vector{Int}()
        while idx <= length(gate_list)
            gate = gate_list[idx]
            if !( any(qubit_used_in_layer[affectedqubits(gate)]) ) 
                push!(layer, gate)
                push!(del_gates, idx)
            end
            qubit_used_in_layer[affectedqubits(gate)] .= true    
            idx += 1
        end
        push!(layers, layer)
        deleteat!(gate_list, del_gates)
    end
    return layers
end


# ------------ Telegates for construction of executable circuit ---------------


function add_telegate(circuit, DQC_control, DQC_target, control_register, target_register, n, noise, p_idle_telegate_layer)
    if control_register == target_register
        throw("Ooops, this is not a telegate.")
    end
    
    # control_comm_idx is num_data_qubits + control_q_per_core*number of control_registers that came before the active register + target register offset
    control_comm_index = n.num_data_qubits+ (n.num_comm_qubits_per_register * (control_register-1)) + (control_register < target_register ? target_register-1 : target_register ) 
    target_comm_index = n.num_data_qubits+  (n.num_comm_qubits_per_register * (target_register-1)) + (target_register < control_register ? control_register-1 : control_register )
    
    # -------- EJPP Protocol --------
    
    # ---- I. Bell pair creation ----

    # This H-CNOT only mimics the way Bell state entanglement is created; in reality, this is achieved via beam splitters or such.
    # In order to account for the differing physical circumstance, we apply a two-qubit depolarising channel with a specific noise probability afterwards
    push!(circuit, sHadamard(control_comm_index))
    push!(circuit, sCNOT(control_comm_index, target_comm_index))
    add_noise(circuit, [control_comm_index, target_comm_index], noise.p_bell; two_qubits = true) # Bell state initialisation noise
    add_noise(circuit, [DQC_control], p_idle_telegate_layer) # idling noise on data qubit (was excluded in the telegate layer idling noise channel)
    add_noise(circuit, [DQC_target], p_idle_telegate_layer) # "

    # ---- II. CNOT(control, comm_c) + CNOT(comm_t, target) ----

    push!(circuit, sCNOT(DQC_control, control_comm_index))
    push!(circuit, sCNOT(target_comm_index, DQC_target))
    add_noise(circuit, [DQC_control, control_comm_index], noise.p; two_qubits = true) # mixed-species noise
    add_noise(circuit, [target_comm_index, DQC_target], noise.p; two_qubits = true) # "


    # ---- III. comm1 + comm2 Measurement ----
 
    #pauli_string_control = build_pauli_string_measurement(n.num_data_and_comm_qubits, [control_comm_index])
    classical_register_index_control = control_comm_index - n.num_data_qubits #sum(n.register_sizes[1:control_register])
    #meas_control = PauliMeasurement(pauli_string_control, classical_register_index_control)
    meas_control = sMRZ(control_comm_index, classical_register_index_control ) 
    # The restoration can probably be neglected since in reality, the Bell pair will be created anew via photonic beam splitters. 
    # For the sake of simulation however, we assume lossless restoration. This is handled by using sMRZ rather than sMZ.
    add_noise(circuit, [control_comm_index], 3/2*noise.p) # Comm qubit measurement noise
    push!(circuit, meas_control)

    classical_register_index_target = target_comm_index - n.num_data_qubits#sum(n.register_sizes[1:target_register])
    meas_target = sMRZ(target_comm_index, classical_register_index_target)#PauliMeasurement(pauli_string_target, classical_register_index_target)
    push!(circuit, sHadamard(target_comm_index)) # measuring in the X-basis requires us to Hadamard-transform and then measure in the Z basis
    add_noise(circuit, [target_comm_index], noise.p) # single-gate noise
    add_noise(circuit, [target_comm_index], 3/2*noise.p) # Comm qubit measurement noise
    push!(circuit, meas_target)

    add_noise(circuit, [DQC_control], noise.p) # usual idle depolarising_noise
    add_noise(circuit, [DQC_target], noise.p) # "

    # ---- IV. Conditional Operations on data qubits ----
    # perform conditional operations, conditioned on the measurement bit: if state.bits[op.controlbit] is true, the measurment yielded eigenvalue -1, if it is false, it yielded +1
    push!(circuit, Types.ConditionalGate(sX(DQC_target),sId1(DQC_target), meas_control.bit))  
    push!(circuit, Types.ConditionalGate(sZ(DQC_control),sId1(DQC_control), meas_target.bit))

    add_noise(circuit, [DQC_control], noise.p) # single-qubit gate noise
    add_noise(circuit, [DQC_target], noise.p)  # "

    # Ideally, we would add noise conditional on whether or not we apply a gate. However, it is acceptable to simply assume that we apply an identity gate
    #add_noise(circuit, [DQC_control], noise.classical_comm_noise) # no matter if X or I applied, we assume some classical communication noise
    #add_noise(circuit, [DQC_target], noise.classical_comm_noise) # 

    
    #push!(circuit, Types.ConditionalGate(sX(control_comm_index),sId1(control_comm_index), meas_control.bit)) # restore the |0> state in the control comm qubit
    #push!(circuit, Types.ConditionalGate(sX(target_comm_index),sId1(target_comm_index), meas_target.bit))  # restore the |0> state in the target comm qubit

    return circuit

end

# function build_pauli_string_measurement(num_qubits::Int, qubits::Vector{Int})
#     pauli = I # we can safely assume that the first qubit is a data qubit, since this is only false whenever there are zero qubits
#     @inbounds for i in 2:(num_qubits) # traverses all data and comm qubits 
#         pauli = (i in qubits) ? pauli⊗Z : pauli⊗I
#     end
#     return pauli
# end


# ------------ Noise functions for executable circuit ---------------

function add_noise(circuit, qubits::Vector{Int}, prob::Float64; two_qubits = false) 
    """Depolarising noise on a set of qubits"""
    if prob<0 || prob > 1
        throw("Please provide a valid noise probabilty in [0,1]")
    end
    if two_qubits
        @assert length(qubits) == 2 "Trying to apply a two-qubit channel to a system of size $(length(qubits))"
        noise = NoiseOp(TwoQubitDepolarisingNoise(prob),qubits);
    else
        noise = NoiseOp(UnbiasedUncorrelatedNoise(prob),qubits);
    end

    push!(circuit, noise)
    return circuit
end


struct TwoQubitDepolarisingNoise{T} <: AbstractNoise
    p::T
end
TwoQubitDepolarisingNoise(p::Integer) = TwoQubitDepolarisingNoise(float(p))

function applynoise!(s::AbstractStabilizer, noise::TwoQubitDepolarisingNoise, indices::Tuple{Int, Int})
    infid = noise.p/15
    i,j = indices[1], indices[2]
    r = rand()
    if r<infid
        apply_single_x!(s,i)
    elseif r<2infid
        apply_single_z!(s,i)
    elseif r<3infid
        apply_single_y!(s,i)
    elseif r<4infid
        apply_single_x!(s,j)
    elseif r<5infid
        apply_single_z!(s,j)
    elseif r<6infid
        apply_single_y!(s,j)
    elseif r<7infid
        apply_single_x!(s,i)
        apply_single_x!(s,j)
    elseif r<8infid
        apply_single_x!(s,i)
        apply_single_z!(s,j)
    elseif r<9infid
        apply_single_x!(s,i)
        apply_single_y!(s,j)
    elseif r<10infid
        apply_single_z!(s,i)
        apply_single_x!(s,j)
    elseif r<11infid
        apply_single_z!(s,i)
        apply_single_y!(s,j)
    elseif r<12infid
        apply_single_z!(s,i)
        apply_single_z!(s,j)
    elseif r<13infid
        apply_single_y!(s,i)
        apply_single_x!(s,j)
    elseif r<14infid
        apply_single_y!(s,i)
        apply_single_y!(s,j)
    elseif r<15infid
        apply_single_y!(s,i)
        apply_single_z!(s,j)
    end
    s
end


#We can also use NoisyGate: https://github.com/QuantumSavory/QuantumClifford.jl/blob/74ee758e87f5d7b1255d6747b346cff15ee10cea/docs/src/noisycircuits_ops.md

function apply!(state::Register, op::Types.ConditionalGate)
    if state.bits[op.controlbit]
        apply!(state, op.truegate)
    else
        apply!(state, op.falsegate)
    end
    return state
end

function affectedqubits(op::Types.ConditionalGate)
    qs = Int[]
    append!(qs, collect(affectedqubits(op.truegate)))
    if op.falsegate !== nothing
       append!(qs, collect(affectedqubits(op.falsegate)))
    end
    return unique(qs)
end

function affectedqubits(op::AbstractSingleQubitOperator)
    qs = Int[]
    append!(qs, op.q)
    return qs
end

function affectedqubits(op::sMZ)
    qs = Int[]
    append!(qs, op.qubit)
    return qs
end

function affectedqubits(op::AbstractTwoQubitOperator)
    qs = Int[]
    append!(qs, op.q1)
    append!(qs, op.q2)
    return qs
end



end

