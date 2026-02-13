module Helper

export create_lookup_array, comm_qubits_array

# lookup arrays needs to be created only once before executing the genetic search
#TODO: replace with cleaner version
function create_lookup_array(register_sizes)
    register_lookup_array = Vector{Int}(undef, sum(register_sizes))
    register_start_indices = Vector{Int}(undef, length(register_sizes))
    register = 1
    register_start_index = 1
    
    for i in eachindex(register_sizes)
        j = register_sizes[i]
        register_lookup_array[register_start_index:register_start_index+j-1] .= register
        register_start_indices[register] = register_start_index
        register_start_index+=j
        register +=1
    end
    #print(register_lookup_array, register_start_indices)
    return register_lookup_array, register_start_indices
end

#function create_lookup_array(params)
# register_lookup_array = Int[]
# for (register, size) in enumerate(params.register_sizes)
#     append!(register_lookup_array, fill(register, size))
# end
#end

function comm_qubits_array(params)
    comm_qubits_array = Vector{Int}(undef,length(params.register_sizes))
    index = 1
    for i in eachindex(params.register_sizes)
        j = params.register_sizes[i]
        comm_qubits_array[i] = index
        index += j
    end
    return comm_qubits_array
end


end