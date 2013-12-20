export ## Types
       ## Methods
       unitary_propagator

CompositeQSystem() = CompositeQSystem(QSystem[], Interaction[], ParametricInteraction[], Vector{Vector{Int}}[], Vector{Vector{Int}}[], Dissipation[])

function getindex(c::CompositeQSystem, key::String)
    for s in c.subSystems
        if label(s) == key
            return s
        end
    end
    throw(KeyError(key))
end

+(s1::QSystem, s2::QSystem) = (c = CompositeQSystem(); c += s1; c += s2; c)

function +(c::CompositeQSystem, q::QSystem)
    append!(c.subSystems, [q])
    update_expansion_indices!(c)
    return c
end

function +(c::CompositeQSystem, i::Interaction)
    append!(c.interactions, [i])
    update_expansion_indices!(c)
    return c
end

function +(c::CompositeQSystem, pi::ParametricInteraction)
    append!(c.parametericInteractions, [pi])
    return c
end

function +(c::CompositeQSystem, d::Dissipation)
    append!(c.dissipators, [d])
end

function update_expansion_indices!(c::CompositeQSystem)
    c.subSystemExpansions = [Vector{Int}[] for _ = 1:length(c.subSystems)]
    for (ct, sys) in enumerate(c.subSystems)
        c.subSystemExpansions[ct] = expand_indices([ct], dims(c))
    end

    c.interactionExpansions = [Vector{Int}[] for _ = 1:length(c.interactions)]
    for (ct, i) in enumerate(c.interactions)
        c.interactionExpansions[ct] = expand_indices(find_subsystem_pos(c, i), dims(c))
    end
end

labels(c::CompositeQSystem) = [label(s) for s in c.subSystems]
dims(c::CompositeQSystem) = [dim(s) for s in c.subSystems]
dim(c::CompositeQSystem) = prod([dim(s) for s in c.subSystems])

function find_subsystem_pos(c::CompositeQSystem, s::QSystem)
    @assert s in c.subSystems "Oops! Subsystem not found."
    findin(c.subSystems, [s])
end

function find_subsystem_pos(c::CompositeQSystem, i::Interaction)
    #Field-system interactions are one-body terms
    if isa(i, SemiClassicalDipole)
        return find_subsystem_pos(c, i.system2)
    elseif isa(i, FluxTransmon)
        return find_subsystem_pos(c, i.transmon)
    else
        return [find_subsystem_pos(c, i.system1), find_subsystem_pos(c, i.system2)]
    end
end

function hamiltonian(c::CompositeQSystem, t::Float64)

    #Initialize Hamiltonian
    Htot = zeros(Complex128, dim(c), dim(c))

    #Add in all the terms    
    hamiltonian_add!(Htot, c, t)

    return Htot
end
hamiltonian(c::CompositeQSystem) = hamiltonian(c, 0.0)

function hamiltonian_add!(Ham::Matrix{Complex128}, c::CompositeQSystem, t::Float64)
    #Fast system hamiltonian calculator with total Hamiltonian preallocated
    
    #Zero the Hamiltonian memory
    Ham[:] = zero(Complex128)

    #Update the subsystems with the parameteric interactions
    for pi in c.parametericInteractions 
        update_params(c, pi, t)
    end

    #Add together subsystem Hamiltonians
    for (subsys, expander) in zip(c.subSystems, c.subSystemExpansions)
        expand_add!(Ham, hamiltonian(subsys, t), expander)
    end

    #Add interactions
    for (i, expander) in zip(c.interactions, c.interactionExpansions)
        expand_add!(Ham, hamiltonian(i,t), expander)
    end
end

function expand(m::Matrix, actingOn::Vector, dims::Vector)
    #Expand an operator onto a larger Hilbert space
    # m: matrix form of  operator
    # actingOn: array of which subsystem index the operator should be acting on 
    # dims: array of dimensions of all the subsystems

    @assert size(m, 1) == prod(dims[actingOn]) "Oops! Dimensions of matrix do not match dims argument."

    #Create the large matrix by tensoring on identity
    l = length(dims)
    eyeIndices = filter(x->!(x in actingOn), 1:l)
    M = kron(m, eye(eltype(m), prod(dims[eyeIndices])))  

    #Reshape into multi-dimensional array given by subsystem dimensions
    #Since we have a matrix we repeat for rows then columns
    M = reshape(M, tuple([dims, dims]...))

    #Permute magic 
    forwardPerm = [actingOn, eyeIndices]
    reversePerm = invperm(forwardPerm)
    #Handle the way tensor product indices work (last subsystem is fastest)
    reversePerm = reverse(l+1-reversePerm)
    M = permutedims(M, [reversePerm, reversePerm+l])

    #Reshape back
    return reshape(M, prod(dims), prod(dims))
end

function expand(m::Matrix, indices::Vector, sizeM::Int )
    M = zeros(eltype(m), (sizeM, sizeM))
    for (ct, inds) in enumerate(indices)
        M[inds] = m[ct]
    end
    return M
end

function expand_add!(M::Matrix, m::Matrix, indices::Vector )
    #Add to certain indices of M with terms from m according to expansion indices.
    for ct=1:length(indices)
        # M[indices[ct]] += m[ct]
        for idx = indices[ct]
            M[idx] += m[ct]
        end
    end
end

function expand_indices(actingOn::Vector, dims::Vector)
    #Calculate the indices for expansion
    actingOnDim = prod(dims[actingOn])
    sm = (actingOnDim, actingOnDim)
    lenm = actingOnDim^2;
    M = expand(reshape(1:lenm, sm), actingOn, dims)
    return [find(M .== x) for x in 1:lenm]
end
