export ## Types
       ## Methods
       unitary_propagator

CompositeQSystem() = CompositeQSystem(QSystem[],
                                      Interaction[],
                                      ParametricInteraction[],
                                      Tuple{Vector{IndexSet},Vector{IndexSet},Vector{IndexSet}}[],
                                      Tuple{Vector{IndexSet},Vector{IndexSet},Vector{IndexSet}}[],
                                      Tuple{Vector{IndexSet},Vector{IndexSet},Vector{IndexSet}}[],
                                      Dissipation[])

function getindex(c::CompositeQSystem, key::String)
    for s in c.subSystems
        if label(s) == key
            return s
        end
    end
    throw(KeyError(key))
end

+(s1::QSystem, s2::QSystem) = (c = CompositeQSystem(); c += s1; c += s2; c)
+(s::QSystem, i::Interaction) = (c = CompositeQSystem(); c += s; c += i; c)

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
    update_expansion_indices!(c)
    return c
end

function update_expansion_indices!(c::CompositeQSystem)
    ddims = [dims(c); dims(c)]
    subsystems = length(c.subSystems)

    c.subSystemExpansions = [(IndexSet[],IndexSet[],IndexSet[]) for _ = 1:length(c.subSystems)]
    for (ct, sys) in enumerate(c.subSystems)
        c.subSystemExpansions[ct] = (expand_indices([ct], dims(c)),            # operator
                                     expand_indices([ct], ddims),              # superoperator right
                                     expand_indices([ct.+subsystems], ddims))   # superoperator left
    end

    c.interactionExpansions = [(IndexSet[],IndexSet[],IndexSet[]) for _ = 1:length(c.interactions)]
    for (ct, i) in enumerate(c.interactions)
        pos_list = find_subsystem_pos(c,i)
        c.interactionExpansions[ct] = (expand_indices(pos_list, dims(c)), # operator
                                       expand_indices(pos_list, ddims),   # superoperator right
                                       expand_indices(map(x->x.+subsystems,pos_list), ddims)) #superoperator left
    end

    c.dissipatorExpansions = [(IndexSet[], IndexSet[], IndexSet[]) for _ = 1:length(c.dissipators)]
    for (ct, d) in enumerate(c.dissipators)
        subsys = find_subsystem_pos(c, d)
        # for efficiency, we need expansion for an operator acting on the left,
        # another for one acting on the right, and one for operators acting on both sides
        c.dissipatorExpansions[ct]= (expand_indices( [subsys.+subsystems;], ddims),  # left
                                     expand_indices( [subsys;], ddims),             # right
                                     expand_indices( [subsys; subsys.+subsystems], ddims)) # bilateral
    end
end

labels(c::CompositeQSystem) = [label(s) for s in c.subSystems]
dims(c::CompositeQSystem) = [dim(s) for s in c.subSystems]
dim(c::CompositeQSystem) = prod([dim(s) for s in c.subSystems])

function find_subsystem_pos(c::CompositeQSystem, s::QSystem)
    @assert s in c.subSystems "Oops! Subsystem not found."
    findin(c.subSystems, [s])
end

find_subsystem_pos(c::CompositeQSystem, i::Interaction) = [find_subsystem_pos(c, i.system1); find_subsystem_pos(c, i.system2)]
find_subsystem_pos(c::CompositeQSystem, i::FluxTransmon) = find_subsystem_pos(c, i.transmon)
find_subsystem_pos(c::CompositeQSystem, i::SemiClassicalDipole) = find_subsystem_pos(c, i.system2)
find_subsystem_pos(c::CompositeQSystem, i::RotatingSemiClassicalDipole) = find_subsystem_pos(c, i.system2)

function find_subsystem_pos(c::CompositeQSystem, d::Dissipation)
    @assert d in c.dissipators "Oops! Dissipator not found in composite system."
    findin(c.subSystems, [d.system])
end

function hamiltonian(c::CompositeQSystem, t::Float64=0.0)

    if length(c.subSystems) != 0
        #Initialize Hamiltonian
        Htot = zeros(Complex128, dim(c), dim(c))

        #Add in all the terms
        hamiltonian_add!(Htot, c, t)

        return Htot
    else
        error("No systems added to composite system yet.")
    end
end

function hamiltonian_add!{T<:Number}(Ham::AbstractMatrix{T}, c::CompositeQSystem, t::Float64)
    #Fast system hamiltonian calculator with total Hamiltonian preallocated

    #Zero the Hamiltonian memory
    if issparse(Ham)
        # TODO: check: is it safe to assume that we will not make liouv denser and denser?
        rows,cols,_ = findnz(Ham)
        for i = 1:length(rows)
            Ham[rows[i],cols[i]] = 0.0
        end
    else
        Ham[:] = 0.0
    end

    #Update the subsystems with the parameteric interactions
    for pi in c.parametericInteractions
        update_params(c, pi, t)
    end

    #Add together subsystem Hamiltonians
    for (subsys, expander) in zip(c.subSystems, c.subSystemExpansions)
        expand_add!(Ham, hamiltonian(subsys, t), expander[1])
    end

    #Add interactions
    for (subsys, expander) in zip(c.interactions, c.interactionExpansions)
        expand_add!(Ham, hamiltonian(subsys,t), expander[1])
    end
end

function liouvillian_dual_add!(liouv::Matrix{Complex128}, c::CompositeQSystem, t::Float64 )
    #Fast system superoperator calculator with liouvillian preallocated

    #Zero the preallocated operators
    liouv[:] = 0.0

    #Update the subsystems with the parameteric interactions
    for pi in c.parametericInteractions
        update_params(c, pi, t)
    end

    #Add together subsystem Hamiltonians
    for (subsys, expander) in zip(c.subSystems, c.subSystemExpansions)
        expand_add!(liouv, transpose(hamiltonian(subsys, t)), expander[2], mult=-1im) # superoperator right
        expand_add!(liouv,           hamiltonian(subsys, t),  expander[3], mult= 1im) # superoperator left
    end

    #Add interactions
    for (subsys, expander) in zip(c.interactions, c.interactionExpansions)
        expand_add!(liouv, transpose(hamiltonian(subsys, t)), expander[2], mult=-1im) # superoperator right
        expand_add!(liouv,           hamiltonian(subsys, t),  expander[3], mult= 1im) # superoperator left
    end

    # Add the Liouvillian for the dissipators
    for (subsys, expander) in zip(c.dissipators, c.dissipatorExpansions)
        expand_add!(liouv, liouvillian_left(subsys,t),  expander[1]) # left
        expand_add!(liouv, liouvillian_right(subsys,t), expander[2]) # right
        expand_add!(liouv, liouvillian_bilat(subsys,t), expander[3]) # bilateral
    end
end

function liouvillian_add!{T<:Number}(liouv::AbstractMatrix{T}, c::CompositeQSystem, t::Float64 )
    #Fast system superoperator calculator with liouvillian preallocated

    #Zero the preallocated operators.
    if issparse(liouv)
        # TODO: check: is it safe to assume that we will not make liouv denser and denser?
        rows,cols,_ = findnz(liouv)
        for i = 1:length(rows)
            liouv[rows[i],cols[i]] = 0.0
        end
    else
        liouv[:] = 0.0
    end

    #Update the subsystems with the parameteric interactions
    for pi in c.parametericInteractions
        update_params(c, pi, t)
    end

    #Add together subsystem Hamiltonians
    for (subsys, expander) in zip(c.subSystems, c.subSystemExpansions)
        expand_add!(liouv, transpose(hamiltonian(subsys, t)), expander[2], mult= 1im) # superoperator right
        expand_add!(liouv,           hamiltonian(subsys, t),  expander[3], mult=-1im) # superoperator left
    end

    #Add interactions
    for (subsys, expander) in zip(c.interactions, c.interactionExpansions)
        expand_add!(liouv, transpose(hamiltonian(subsys, t)), expander[2], mult= 1im) # superoperator right
        expand_add!(liouv,           hamiltonian(subsys, t),  expander[3], mult=-1im) # superoperator left
    end

    # Add the Liouvillian for the dissipators
    for (subsys, expander) in zip(c.dissipators, c.dissipatorExpansions)
        expand_add!(liouv, liouvillian_left(subsys,t),  expander[1]) # left
        expand_add!(liouv, liouvillian_right(subsys,t), expander[2]) # right
        expand_add!(liouv, liouvillian_bilat(subsys,t)', expander[3]) # bilateral
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
    M = isempty(eyeIndices) ? m : kron(m, eye(eltype(m), prod(dims[eyeIndices])))

    #Reshape into multi-dimensional array given by subsystem dimensions
    #Since we have a matrix we repeat for rows then columns
    M = reshape(M, tuple([dims; dims]...))

    #Permute magic
    forwardPerm = [actingOn; eyeIndices]
    reversePerm = invperm(forwardPerm)
    #Handle the way tensor product indices work (last subsystem is fastest)
    reversePerm = reverse((l+1) .- reversePerm)
    M = permutedims(M, [reversePerm; reversePerm .+ l])

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

function expand_add!{T<:Number,U<:Number}(M::AbstractMatrix{T}, m::AbstractMatrix{U}, indices::Vector; mult=1.0 )
    #Add to certain indices of M with terms from m according to expansion indices.
    for ct=1:length(indices)
        # M[indices[ct]] += m[ct]
        for idx = indices[ct]
            M[idx] += mult*m[ct]
        end
    end
end

function expand_indices(actingOn::Vector, dims::Vector)
    #Calculate the indices for expansion
    actingOnDim = prod(dims[actingOn])
    sm = (actingOnDim, actingOnDim)
    lenm = actingOnDim^2;
    M = expand(reshape(1:lenm, sm), actingOn, dims)
    return IndexSet[find(M .== x) for x in 1:lenm]
end

function expm_eigen(A::Matrix, t)
    #Calculates expm(t*A) via eigenvalue decomposition and assuming Hermitian matrix
    F = eigfact(Hermitian(A))

    # V * diagm(exp(t*D)) * V'
    return scale(F[:vectors], exp(t*F[:values])) * F[:vectors]'
end

function unitary_propagator(sys::CompositeQSystem, timeStep::Float64, startTime::Float64, endTime::Float64)

    #Preallocate Hamiltonian memory
    Ham = zeros(Complex128, (dim(sys), dim(sys)))
    Uprop = @parallel (*) for time = startTime:timeStep:endTime
        #a *= b expands to a = a*b
        hamiltonian_add!(Ham, sys, time)
        expm_eigen(Ham, 1im*timeStep)
    end
    return Uprop'
end
