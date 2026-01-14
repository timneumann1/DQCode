# For convenient graph data structures
using Graphs

# For discrete event simulation
using ResumableFunctions
using ConcurrentSim

# Useful for interactive work
# Enables automatic re-compilation of modified codes
using Revise

# The workhorse for the simulation
using QuantumSavory

using QuantumSavory.ProtocolZoo: EntanglerProt

bell = StabilizerState("XX ZZ")
noisy_pair_func(F) = (1-F)*MixedState(bell) + F*projector(bell) # TODO: Is this general enough?

"""Creates the datastructures representing the simulated network"""
function simulation_setup(
    sizes, # Array giving the number of qubits in each node
    T2 # T2 dephasing times for the qubits
    ;
    representation = QuantumOpticsRepr # Representation to use for the qubits
    )
    R = length(sizes) # Number of registers

    # All of the quantum register we will be simulating
    registers = Register[]
    for s in sizes
        traits = [Qubit() for _ in 1:s]
        repr = [representation() for _ in 1:s]
        bg = [T2Dephasing(T2) for _ in 1:s]
        push!(registers, Register(traits,repr, bg)) 
    end

    # A graph structure defining the connectivity among registers
    # It is not necessary to use such a structure, however, it is a convenient way to
    # store data about the simulation (and we have created helper plotting functions
    # expecting such a structure).
    graph = grid([R])
    network = RegisterNet(graph, registers) # A graphs with extra "meta data"

    # The scheduler datastructure for the discrete event simulation
    sim = get_time_tracker(network)

    # Add a register datastructures and event locks to each node.
    for v in vertices(network)
        # Create an array specifying whether a qubit is entangled with another qubit
        network[v,:enttrackers] = Any[nothing for i in 1:sizes[v]]
    end

    sim, network
end


@resumable function entangler(
    sim::Environment,   # The scheduler for all simulation events
    network,            # The graph of quantum nodes
    nodea, nodeb,       # The two nodes which we will be entangling
    noisy_pair,         # A raw entangled pair
    entangler_wait_time,# The wait time in case all qubits are "busy"
    entangler_busy_time # How long it takes to establish entanglement
    )
    while true
        ia = findfreequbit(network, nodea)
        ib = findfreequbit(network, nodeb)
        if isnothing(ia) || isnothing(ib)
            @yield timeout(sim, entangler_wait_time)
            continue
        end
        slota = network[nodea,ia]
        slotb = network[nodeb,ib]
        @yield request(slota) & request(slotb)
        registera = network[nodea]
        registerb = network[nodeb]
        @yield timeout(sim, entangler_busy_time)
        initialize!((registera[ia],registerb[ib]),noisy_pair; time=now(sim))
        network[nodea,:enttrackers][ia] = (node=nodeb,slot=ib)
        network[nodeb,:enttrackers][ib] = (node=nodea,slot=ia)
        @simlog sim "entangled node $(nodea):$(ia) and node $(nodeb):$(ib)"
        unlock(slota)
        unlock(slotb)
    end
end

"""Find an uninitialized unlocked qubit on a given node"""
function findfreequbit(network, node)
    register = network[node]
    regsize = nsubsystems(register)
    findfirst(i->!isassigned(register,i) && !islocked(register[i]), 1:regsize)
end