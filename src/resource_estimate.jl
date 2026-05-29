module ResourceEstimation


using ..Types
using ..Helper
using ..DQCodeSimulator: build_layers

using QuantumClifford: Stabilizer, AbstractOperation, sCNOT, AbstractSingleQubitOperator, AbstractTwoQubitOperator, AbstractMeasurement
using Serialization
using CSV, DataFrames

export estimate_resources_encoding_circuit, estimate_resources_measurement_based_encoding

function estimate_resources_encoding_circuit(folder, qpu_sizes, ancilla_map)

    # --------- FT Encoding Circuit ---------------
    # Retrieve the data + verification circuit and analyse how many single-, two-qubit and telegates + measurements are needed 
    
    #SIMPLY RETRIEVE DATA FROM DATA file here

    #verification_circ = deserialize(joinpath(folder, "verification_circuit.jls"))
    #data_circ = deserialize(joinpath(folder, "data_circuit.jls") )
    df = CSV.read(joinpath(folder, "dqc_sim_data.csv"), DataFrame)
    depth_cx_layers = df.depth_cx_layers[1]
    depth_telegate_layers = df.depth_telegate_layers[1]
    total_single_qubit_count = df.total_single_qubit_count[1]
    total_two_qubit_count = df.total_two_qubit_count[1]
    total_telegate_count = df.total_telegate_count[1]
    #total_gate_counts = df.gate_counts[1]
    @info "total gate counts:$total_single_qubit_count,$total_two_qubit_count $total_telegate_count"
    num_meas = df.num_meas[1]
    #println(data_circ)
    #@info "Verification circuit: \n $verification_circ"
    #println(network_specs)
    #data_circ_gc = gate_counts(data_circ, network_specs)
    #@info "Gate counts data circuit: $data_circ_gc"
    #verif_circ_gc, num_meas = gate_counts(verification_circ, network_specs, ancilla_map)
    #@info "Gate counts verification circuit: $verif_circ_gc"
    #total_gate_counts = data_circ_gc + verif_circ_gc
    #@info "total gate counts: $total_gate_counts"
    #total_number_measurements = num_meas
    @info "The FT encoding circuit requires $(length(ancilla_map)) ancillas, $total_single_qubit_count single qubit gates, $total_two_qubit_count CX gates, $total_telegate_count telegates as well as $num_meas measurements."
    ancilla_distr = zeros(Int, length(qpu_sizes))
    for core in ancilla_map # iterate over all entries in the ancilla map list
        ancilla_distr[core] += 1 # increase the respective core in ancilla_dist count of the  count(core == idx, core in ancilla_map)
    end
    qpu_core_sizes = copy(qpu_sizes) .+ ancilla_distr 
    @info "The total qubit counts per core is $qpu_core_sizes "
    return length(ancilla_map), (depth_cx_layers, depth_telegate_layers), (total_single_qubit_count, total_two_qubit_count, total_telegate_count), num_meas, qpu_core_sizes
end
   
    
   # ---------- Measurement-based initialisation -------------

function estimate_resources_measurement_based_encoding(network_specs, code_params, qpu_sizes)

    # building the whole circ incl cores, then building layers asap, then executing them based on comm qubits 
    #TODO:only X checks: we map the qubits based on all stabilisers, but for fairness of comparison, we actually only need to measure the X checks  for the encoding part, since we may assume that we start in the all zero state (and logical Zs and Z checks stabilise that state already)
    # count depth and telegates based on same networking assumptions as encoding circuit (here, having two overlapping telegates is unlikely since this requires that two ancillas are mapped to the same core)
    gate_counts = [0,0,0]
    num_meas = 0
    qpu_core_sizes = copy(qpu_sizes)
    stabilisers = Stabilizer(code_params.qec_code)
    @info "$stabilisers"
    lZs = code_params.logical_Zs

    measurements = collect(stabilisers)
    append!(measurements, lZs)
    #vcat(stabilisers, lZs)
    @info "measurements: $measurements"

    # count gate count in this loop, then do a separate one for depth (with build layers)
    measurement_circ = Vector{AbstractOperation}()
    num_ancillas = 0
    ancilla_map = Vector{Int}()

    for (stab_meas_idx, check) in enumerate(measurements)
        msm_cores =  Dict{Int, Int}() # will track which core has the most interactions (can be done per check since each check consumes one ancilla)
        # we assume all positive phases       
        pauli = check
        pauli_support = Vector{Int}()
        pauli_type = "X"
        for idx in collect(1:length(pauli)) 
            qubit = pauli[idx]
            if qubit != (false, false)
                core = network_specs.register_lookup_array[network_specs.inv_map[idx]]
                msm_cores[core] = get!(msm_cores, core, 0) + 1
                push!(pauli_support, idx)
                if qubit == (false, true) # this can never be falsely overwritten, since in CSS codes the stabiliser checks (and (Z) logicals! ) don't mix (at least not in the canonical choice we use)
                    pauli_type = "Z"
                end
            end
        end
        @info "$check: $msm_cores"
        if pauli_type=="X"
            gate_counts[1] += 2 # we need to sandwich hadamards
            num_meas += 1 # and measure the ancilla
            num_ancillas += 1
            # cnots/telegates will be counted later
        else # we only apply the X type checks
            continue
        end
        mapped_core = first(keys(msm_cores))
        mapped_count = msm_cores[mapped_core] 
        for (c, cnt) in msm_cores
            if cnt > mapped_count # ties are broken in order of appearance in dict 
                mapped_core = c
                mapped_count = cnt
            end
        end
        qpu_core_sizes[mapped_core] +=1
        push!(ancilla_map, mapped_core)

        for idx in pauli_support
            mapped_q = network_specs.inv_map[idx]
            if network_specs.register_lookup_array[mapped_q] == mapped_core
                gate_counts[2] += 1
                @info "CX($idx, $(network_specs.num_data_qubits+stab_meas_idx)) is intra-core"
            else 
                gate_counts[3] += 1
                @info "CX($idx, $(network_specs.num_data_qubits+stab_meas_idx)) is telegate"
            end
            push!(measurement_circ, sCNOT(mapped_q, network_specs.num_data_qubits+stab_meas_idx))
        end
        
        @info "Current number of gate counts is $gate_counts"
        # now flag the qubits and target core ancilla qubit that have been used, if there is any overlap 
        
    end
    @info "Circuit is $measurement_circ"
    # MEASUREMTN CIRC CONTAINS THE MAPPED QUBITS ALREADY
    layers_measurement_initialisation = build_layers(measurement_circ, network_specs.num_data_qubits+num_ancillas)
    for layer in layers_measurement_initialisation
        @info "Layers are $layer"
    end
        #Next, determine layers: (depth_cx_layers, depth_telegate_layers)
    depth = [0,0]
    gate_counts_meas = [0,0,0]

    @info "ancilla_map: $ancilla_map"
    for layer in layers_measurement_initialisation
    
        num_telegates_layer = 0
        telegate_pairs = Set{Tuple{Int,Int}}()
        #telegate_qubits = Set{Int64}()

        for gate in layer
            T = typeof(gate)
            if T <: AbstractSingleQubitOperator
                #qubit = affectedqubits(gate)[1]
                #DQC_qubit = n.inv_map[qubit]
                #push!(circuit, sHadamard(DQC_qubit))
                gate_counts_meas[1] += 1
                #add_noise(circuit, [DQC_qubit], noise.p_single) # Single-qubit noise
                
            elseif T <: AbstractTwoQubitOperator
                control = gate.q1 #affectedqubits(gate)[1]
                target =  gate.q2#affectedqubits(gate)[2]

                # the control is always a data qubit, the target always an ancilla
                #DQC_control = network_specs.inv_map[control]
                control_register = network_specs.register_lookup_array[control] 
                
                
                target_register = ancilla_map[target-network_specs.num_data_qubits] 
                
                # if control > network_specs.num_data_qubits 
                #     DQC_control = control
                #     control_register = ancilla_map[control-network_specs.num_data_qubits] 
                # else
                #     DQC_control = network_specs.inv_map[control]
                #     control_register = network_specs.register_lookup_array[DQC_control] 
                # end
                # if target  > network_specs.num_data_qubits 
                #     DQC_target = target
                #     target_register = ancilla_map[target-network_specs.num_data_qubits] 
                # else
                #     DQC_target = network_specs.inv_map[target]
                #     target_register = network_specs.register_lookup_array[DQC_target] 
                # end
                
                # @assert DQC_control > 0 
                # @assert DQC_target > 0 
                # @assert control_register > 0 
                # @assert target_register > 0 



                if control_register == target_register 
                    gate_counts_meas[2] += 1
                    @info "CX($control, $target) is intra-core"
                else
                    if num_telegates_layer == 0
                        num_telegates_layer += 1 # if this is the first telegate in the layer, we need to increase the telegate count; for any of the following telegates, we increase the number of layers only when the comm pair has already been used (since then we must wait)
                    end
                    #telegates_layer = true
                    if (control_register, target_register) ∈ telegate_pairs || (target_register, control_register) ∈ telegate_pairs # if the communication qubits for this register pair have been used, we need to increase the number of telegates for this layer, and apply as much noise as there have been previous telegates in this layer
                        #for _ in 1:num_telegates_layer
                        #add_noise(circuit, [DQC_control, DQC_target], 1-(1-noise.p_idle_telegate_layer)^num_telegates_layer ) # current telegate qubit obtain noise from previous telegate_layers (where they were not yet part of the telegate qubits by construction)
                        #add_noise(circuit, collect(telegate_qubits), noise.p_idle_telegate_layer) # other telegate qubits (used so far) incur noise
                        #end
                        num_telegates_layer +=1
                    end 
                    #circuit = add_telegate(circuit, DQC_control, DQC_target, control_register, target_register, n, noise)
                    #telegates_layer = true
                    #num_telegates_layer +=1
                    gate_counts_meas[3] += 1
                    @info "CX($control, $target) is telegate"
                    push!(telegate_pairs, (control_register, target_register))
                    #push!(telegate_qubits, DQC_control)
                    #push!(telegate_qubits, DQC_target)
                end
            elseif T <: AbstractMeasurement
                #add_noise(circuit, affectedqubits(gate), noise.p) # measurement noise the ancilla is of the same type, thus we have the same measurement noise
                #push!(circuit, gate) # only ancillas are ever measured, so we don't need a remapping
                num_meas += 1
            else
                throw("Circuit contains gates that have not been classified as Single- or Two-Qubit gate so far.")
            end
        end

        @info "After this layer, the gate count is $gate_counts_meas"
        if num_telegates_layer == 0
            #add_noise(circuit, idle_qubits_DQC, noise.p_idle) # idling noise (if telegate: to all non-telegate qubits, no matter if idle or not (since telegate is much longer), otherwise: idling noise p on idle qubits)
            depth[1] +=1
        else
            # the telegate qubits are excluded from this noise application on the genuine idle qubits (here, noise is also applied on the ancillas)
            #add_noise(circuit, setdiff(vcat(1:n.num_data_qubits, n.num_data_and_comm_qubits+1:length(all_qubits)), Set(telegate_qubits)), 1-(1-noise.p_idle_telegate_layer)^num_telegates_layer) 
            # the set diff targets all non-telegate qubits (since noise is uncorrelated here, DQC vs. usually indices does not matter)
            depth[2] += num_telegates_layer 
            #@info "Number of telegates in this layer: $layer is $num_telegates_layer"
        end
        #telegates_layer ? add_noise(circuit, idle_qubits_DQC, noise.p_idle_telegate_layer) : add_noise(circuit, idle_qubits_DQC, noise.p) # idling noise
#        telegates_layer ? add_noise(circuit, setdiff(vcat(1:n.num_data_qubits, n.num_data_and_comm_qubits+1:length(all_qubits)), Set(tele_qubits)), noise.p_idle_telegate_layer) : add_noise(circuit, idle_qubits_DQC, noise.p_idle) # idling noise (if telegate: to all non-telegate qubits (data + ancilla), no matter if idle or not (since telegate is much longer), otherwise: idling noise p on idle qubits)

    end

    @info gate_counts
    @info gate_counts_meas

    @assert gate_counts_meas[2] == gate_counts[2] && gate_counts_meas[3] == gate_counts[3]


    return num_ancillas, depth, gate_counts, num_meas, qpu_core_sizes

end



   

end

