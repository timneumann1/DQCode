module Parameters

using ..Types
using ..Simulation
using ..Helper

using GLMakie
GLMakie.activate!()

using BenchmarkTools

export run_parameter_sweep


function retrieve_parameters(depolarising_noise_time)
    params = SimulationParameters(
        [5,4], #register sizes
        #12000.0,#T1
        depolarising_noise_time, #4200 T2
        20e-6, # Execution Time of a single qubit gate   #20e^-6
        200e-6,   # Two-qubit gates 200e^-6
        1e-5,  #1e^-5  projective measurement time
        1e-2,  # 1e^-2 classical comm time
        1,#0.9689, # Bell state fidelity
        1.41e-4, # Bell state generation,from [Main, 2025]    success probability 
        1.168e-9 #1.168e-9  # attempt time

    #TODO: Define units and insert realistic values
    #TODO: Add state preparation fidelity and single-shot readout of 99.93% [Harty], single-qubit gate fidelity of 99.99916%, two-qubit fidelity of 99.97% [Löschnauer]
    # characteristic_time = 1000
    # p = 1-exp(-1/characteristic_time) # define probability for Pauli Noise application (Poisson point process)
    #TODO: Define single qubit error rate
    )
    return params
end

#TODO: sweep over parameters and get fidelities

function steane_encoding_circuit(params)

    circuit = Circuit(sum(params.register_sizes), 9)   # params.register_sizes rows (qubits) and 8 columns (time steps)
    
    circuit.gates[2,1] = HadamardGate()
    circuit.gates[3,1] = HadamardGate()
    circuit.gates[5,1] = HadamardGate()

    circuit.gates[2,2] = circuit.gates[4,2] = CNOT_Gate(2,4)
    circuit.gates[5,2] = circuit.gates[8,2] = CNOT_Gate(5,8)

    circuit.gates[3,3] = circuit.gates[9,3] = CNOT_Gate(3,9)

    circuit.gates[2,4] = circuit.gates[7,4] = CNOT_Gate(2,7)

    circuit.gates[5,5] = circuit.gates[9,5] = CNOT_Gate(5,9)

    circuit.gates[3,6] = circuit.gates[8,6] = CNOT_Gate(3,8)

    circuit.gates[2,7] = circuit.gates[9,7] = CNOT_Gate(2,9)

    circuit.gates[3,8] = circuit.gates[4,8] = CNOT_Gate(3,4)
    circuit.gates[5,8] = circuit.gates[7,8] = CNOT_Gate(5,7)

    circuit.gates[2,9] = HadamardGate()
    circuit.gates[3,9] = HadamardGate()
    circuit.gates[4,9] = HadamardGate()
    circuit.gates[5,9] = HadamardGate()
    circuit.gates[7,9] = HadamardGate()
    circuit.gates[8,9] = HadamardGate()
    circuit.gates[9,9] = HadamardGate()

    return circuit
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
    register_lookup_array = create_lookup_array(params)      # create lookup array
    circuit = steane_encoding_circuit(params)                # build initial circuit
    #TODO: block all communication qubit layers! Can be done via row check != comm_qubits,
    #TODO: Include check for no overlaps within one layer
    
    depolarising_times = collect(1e13:1:1e13+1)
    state_fidelities = Float64[]

    #TODO: Add standard deviation 
    #TODO: sweep over other parameters as well
    for depolarising_time in depolarising_times
        params = retrieve_parameters(depolarising_time)  # retrieve parameters, changing bell state fidelity in every sweep
        fid_depol = 0
        num_runs = 1
        for _ in 1:num_runs
            sim_fid = run_simulation(params, circuit, register_lookup_array)
            #@btime run_simulation(params, circuit, register_lookup_array) #sim_fid = run_simulation(params, circuit, register_lookup_array)
            fid_depol += sim_fid.fidelity
        end
        fid_depol = fid_depol/num_runs
        push!(state_fidelities, fid_depol)
        @info "For depolarising tau $depolarising_time, we obtain final state fidelity $(fid_depol)"
    end

    plot_sweep(depolarising_times, state_fidelities, params) # !! takes the last params iteration
end


end
