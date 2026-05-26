module ResourceEstimation


using ..Types
using ..Helper

using QuantumClifford: Stabilizer
using Serialization

export estimate_resources_encoding_circuit, estimate_resources_measurement_based_encoding

function estimate_resources_encoding_circuit(folder, network_specs, qpu_sizes, ancilla_map)

    # --------- FT Encoding Circuit ---------------
    # Retrieve the data + verification circuit and analyse how many single-, two-qubit and telegates + measurements are needed 
    
    verification_circ = deserialize(joinpath(folder, "verification_circuit.jls"))
    data_circ = deserialize(joinpath(folder, "data_circuit.jls") )

    println(data_circ)
    @info "Verification circuit: \n $verification_circ"
    println(network_specs)
    data_circ_gc = gate_counts(data_circ, network_specs)
    @info "Gate counts data circuit: $data_circ_gc"
    verif_circ_gc, num_meas = gate_counts(verification_circ, network_specs, ancilla_map)
    @info "Gate counts verification circuit: $verif_circ_gc"
    total_gate_counts = data_circ_gc + verif_circ_gc
    total_number_measurements = num_meas
    @info "The FT encoding circuit requires $(length(ancilla_map)) ancillas, $(total_gate_counts[1]) single qubit gates, $(total_gate_counts[2]) CX gates, $(total_gate_counts[3]) telegates as well as $total_number_measurements measurements."
    ancilla_distr = zeros(Int, length(qpu_sizes))
    for core in ancilla_map # iterate over all entries in the ancilla map list
        ancilla_distr[core] += 1 # increase the respective core in ancilla_dist count of the  count(core == idx, core in ancilla_map)
    end
    qpu_core_sizes = qpu_sizes .+ ancilla_distr 
    @info "The total qubit counts per core is $qpu_core_sizes "
    return length(ancilla_map), total_gate_counts, total_number_measurements, qpu_core_sizes
end
   
    
   # ---------- Measurement-based initialisation -------------

function estimate_resources_measurement_based_encoding(network_specs, code_params, qpu_sizes)

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
    num_ancillas = length(measurements)
    for check in measurements
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
            gate_counts[1] += 2
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
        for idx in pauli_support
            if network_specs.register_lookup_array[network_specs.inv_map[idx]] == mapped_core
                gate_counts[2] += 1
            else 
                gate_counts[3] += 1
            end
        end
        num_meas += 1
        @info "Current number of gate counts is $gate_counts"
    end
        
    
    return num_ancillas, gate_counts, num_meas, qpu_core_sizes
end
   

end

