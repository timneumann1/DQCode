include("helper_functions.jl") 
using GLMakie
GLMakie.activate!()

# System Parameters
sizes = [6,3]        # Number of qubits in each register
T1 = 1000.00         # T1 relaxation time of all qubits
T2 = 100000.0         # T2 dephasing time of all qubits
F = 1              # Fidelity of the raw Bell pairs

# TODO: Define units and insert realistic values
entangler_wait_time = 0.1  # How long to wait if all qubits are busy before retrying entangling
entangler_busy_time = 1.0  # How long it takes to establish a newly entangled pair
swapper_wait_time = 0.1    # How long to wait if all qubits are unavailable for swapping
swapper_busy_time = 0.15   # How long it takes to swap two qubits
purifier_wait_time = 0.15  # How long to wait if there are no pairs to be purified
purifier_busy_time = 0.2   # How long the purification circuit takes to execute
single_qubit_gate_execution_time = 1  # Execution Time of a single qubit gate
two_qubit_gate_execution_time = 3     # Execution Time of a two qubit gate
projective_measurement_time = 0.5     # Time to peform a measurement
noisy_pair = noisy_pair_func(F)       # Retrieving a (noisy) Bell pair with fidelity F
init_time = 1                         # Time to initialise the system (e.g., in the all-zero state)

@resumable function apply_single_qubit_ops(sim, network, register_array, qubit_array, ops_array)
    # Accepts an array of single-qubit operations and executes them in parallel (during on time slice)
    num_ops = length(register_array)
    @assert num_ops === length(qubit_array) === length(ops_array)
    for i in 1:num_ops
        r_i = register_array[i]
        q_i = qubit_array[i]
        op_i = ops_array[i]
        @yield request(network[r_i,q_i]) 
        apply!(network[r_i,q_i], op_i)
        unlock(network[r_i, q_i])
    end
    @yield timeout(sim, single_qubit_gate_execution_time)
end

@resumable function apply_CNOTs(sim, network, register_array, control_qubits, target_qubits)
    # Accepts an array of control and target qubits for CNOT applications
    # that can be executed in parallel (during on time slice) within one core
    num_ops = length(register_array)
    @assert num_ops === length(control_qubits) === length(target_qubits)
    for i in 1:num_ops
        r_i = register_array[i]
        control_i = control_qubits[i]
        target_i = target_qubits[i]
        @yield request(network[r_i,control_i]) 
        @yield request(network[r_i,target_i]) 
        apply!( (network[r_i,control_i], network[r_i,target_i]), CNOT)
        unlock(network[r_i, control_i])
        unlock(network[r_i, target_i])
    end
    @yield timeout(sim, two_qubit_gate_execution_time)
end

@resumable function apply_remote_CNOTs(sim, network, ri, rj, q1, q2)
    # Accepts control and target qubit for remote CNOT application across cores
    @assert network[ri, 1] !== nothing "Missing qubit ri,1"
    @assert network[rj, 1] !== nothing "Missing qubit rj,1"
    @assert network[ri, q1] !== nothing "Missing control qubit"
    @assert network[rj, q2] !== nothing "Missing target qubit"
    @assert q1 !== 1 "First qubit is not a memory qubit"
    @assert q2 !== 1 "First qubit is not a memory qubit"
        
    # Create Bell Pair between first qubits of registers i and j
    eprot = EntanglerProt(sim, network, ri, rj; pairstate=noisy_pair, chooseslotA=1, chooseslotB=1, rounds=5, success_prob=1.)
    @process eprot()
    @yield timeout(sim, entangler_busy_time)
    
    # Perform CNOT telegate between qubits q1 and q2 (TODO: extend to application of more gates, cf. EJPP protocol)
    @yield request(network[ri,q1]) 
    @yield request(network[rj,q2])
    @yield request(network[ri,1]) 
    @yield request(network[rj,1])
    @simlog sim "Executing Telegate-CNOT between ($ri, $q1) and ($rj, $q2)"
    apply!((network[ri,q1], network[ri,1]), CNOT)  
    @yield timeout(sim, two_qubit_gate_execution_time)
    m1 = project_traceout!(network[ri,1], Z)
    @yield timeout(sim, projective_measurement_time)
    # If the result is '1' (m1 == 2), flip qubit 1 in rj
    # TODO: Add communication channel
    if m1 == 2
        apply!(network[rj,1], X)
        @yield timeout(sim, single_qubit_gate_execution_time)
    end
    apply!((network[rj,1], network[rj,q2]), CNOT)
    @yield timeout(sim, two_qubit_gate_execution_time)
    apply!(network[rj,1], H)
    @yield timeout(sim, single_qubit_gate_execution_time)
    m2 = project_traceout!(network[rj,1], Z)
    @yield timeout(sim, projective_measurement_time)
    # If the result is '1' (m1 == 2), apply Z to q1 in ri
    if m2 == 2
        apply!(network[ri,q1], Z)
        @yield timeout(sim, single_qubit_gate_execution_time)
    end
    unlock(network[ri,q1])
    unlock(network[rj,q2])
    unlock(network[ri,1])
    unlock(network[rj,1])
end

@resumable function circuit(sim, network)
    # Initialize memory qubits in zero state
    initialize!([network[i,j] for i in 1:length(sizes) for j in 2:sizes[i]], SProjector(Z1⊗Z1⊗Z1⊗Z1⊗Z1⊗Z1⊗Z1)  ) 

    # TODO: For each layer, identify the longest execution time among parallel operations
    @yield @process apply_single_qubit_ops(sim, network, [1,1,1],[2,3,5],[H,H,H])
    @yield @process apply_CNOTs(sim, network, 1, 2, 4)
    @yield @process apply_remote_CNOTs(sim, network, 1, 2, 5, 2)
    @yield @process apply_remote_CNOTs(sim, network, 1, 2, 3, 3)
    @yield @process apply_CNOTs(sim, network, 1, 2, 6)
    @yield @process apply_remote_CNOTs(sim, network, 1, 2, 5, 3)
    @yield @process apply_remote_CNOTs(sim, network, 1, 2, 3, 2)
    @yield @process apply_remote_CNOTs(sim, network, 1, 2, 2, 3)
    @yield @process apply_CNOTs(sim, network, 1, 3, 4)
    @yield @process apply_CNOTs(sim, network, 1, 5, 6)
end

# Initialise the network
sim, network = simulation_setup(sizes, T2; representation = QuantumOpticsRepr)  
@process circuit(sim, network)

step_ts = range(0, 75, step=0.2)

# Plotting capabilities
fig = Figure(size=(800,400))
_,ax,_,obs = registernetplot_axis(fig[1,1],network)
ts = Observable(Float64[0])
fidelities = Observable(Float64[0])
ax_fid = Axis(fig[1,2][1,1], xlabel="time", ylabel="Steane-7 Stabilizer State Fidelity")
lZ1 = stairs!(ax_fid,ts,fidelities,label="X")
xlims!(0, nothing)
ylims!(-1.05, 1.05)
Legend(fig[1,2][2,1],[lZ1],["Fidelity"],
            orientation = :horizontal, tellwidth = false, tellheight = true)

display(fig)

steane_7_state = StabilizerState("ZIZIZIZ XIXIXIX IZZIIZZ IXXIIXX IIIZZZZ IIIXXXX ZZZZZZZ") 

record(fig, "steane-7-fidelity.mp4", step_ts, framerate=1, visible=true) do t
    run(sim, t)
    
    fidelity = real(observable(vcat( [network[1,i] for i in 2:6],[network[2,j] for j in 2:3]), SProjector(steane_7_state)))
    print("Fidelity: $fidelity at time $t \n")
    

    push!(fidelities[],fidelity)
    push!(ts[],t)
    ax.title = "t=$(t)"
    notify(obs)
    notify(ts)
    xlims!(ax_fid, 0, t+0.5)
end
