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
    return c
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

function hamiltonian(c::CompositeQSystem, t::Float64=0.0)

    #Initialize Hamiltonian
    Htot = zeros(Complex128, dim(c), dim(c))

    #Add in all the terms    
    hamiltonian_add!(Htot, c, t)

    return Htot
end

function hamiltonian_add!(Ham::Matrix{Complex128}, c::CompositeQSystem, t::Float64)
    #Fast system hamiltonian calculator with total Hamiltonian preallocated
    
    #Zero the Hamiltonian memory
    Ham[:] = 0.0

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

# [todo] - need a liouvillian counterpart to hamiltonian_add!
## due to how parallel multiplies matrices, we need the dual to the liouvillian generators
# function dual_liouvillian_add!(Liou::Matrix{Complex128}, c::CompositeQSystem, t::Float64)
#     #Fast toal system dual Liouvillian calculator with Liouvillian preallocated

#     #Zero the Liouvillian memory
#     Liou[:] = zero(Complex128)

#     #Update the subsystems with the parameteric interactions
#     for pi in c.parametericInteractions 
#         update_params(pi, t)
#     end

#     #Add together subsystem Hamiltonian liouvillians
#     for (ct, s) in enumerate(c.subSystems)
#         expand_add!(Liou, (hamiltonian(s, t),), c.subSystemExpansions[ct])
#     end

#     #Add interactions
#     for (ct, i) in enumerate(c.interactions)
#         expand_add!(Liou, (hamiltonian(i,t),), c.interactionExpansions[ct])
#     end

#     #Add dissipators
#     for (ct, d) in enumerate(c.dissipators)
#         for g in generator(d,t)
#             expand_add!(Liou, generator(d, t), c.dissipatorExpansions[ct])
#         end
#     end
# end

function lindbladian(c::CompositeQSystem, t::Real=0.)
  local H, Gsop, id, d
  H = QSimulator.hamiltonian(c, t)
  d = size(H,1)
  id = eye(d)

  Gsop = QIP.hamiltonian(H)

  for (ct, diss) in enumerate(c.dissipators)
      for sop in generator(diss,t)
          add!(Gsop, 
               QIP.liou(expand(sop[1],c.dissipatorExpansions[ct],d),
                        expand(sop[2],c.dissipatorExpansions[ct],d)))
      end
  end
  return Gsop
end

# function liouville_propagator(sys::CompositeQSystem, timeStep::Float64, startTime::Float64, endTime::Float64)
#     local lind, liouprop
#     #Preallocate Liouvillian memory
#     lind = zeros(Complex128, (dim(sys)^2, dim(sys)^2))
#     liouprop = @parallel (*) for time = startTime:timeStep:endTime
#         dual_liouvillian_add!(lind, sys, time)
#         expm(lind, timeStep)
#     end
#     return lind'
# end

