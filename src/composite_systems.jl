export CompositeQSystem,
       hamiltonian

# tensor products of quantum systems
mutable struct CompositeQSystem
    subsystems::Vector
    fixed_Hs::Vector # tuple of Matrix and exansion indices
    parametric_Hs::Vector # tuple of Functions and expansion indices
    dim::Int
end

CompositeQSystem(qs) = CompositeQSystem(qs, [], [], prod(dim(q) for q in qs))

""" Calculate the drift or natural Hamiltonian of a CompositeQSystem """
function hamiltonian(cqs::CompositeQSystem)
    ham = zeros(Complex128, (dim(cqs), dim(cqs)))
    for h = cqs.fixed_Hs
        expand_add!(ham, new_ham, expand_idxs)
    end
end


""" In place addition of an operator embedded into a larger Hilbert space given a set of expansion indices"""
function expand_add!(op, new_op, expand_idxs)
    for ct = 1:length(expand_idxs)
        for idx = expand_idxs[ct]
            op[idx] += new_op[ct]
        end
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
    M = expand(reshape([1:lenm;], sm), actingOn, dims)
    return IndexSet[find(M .== x) for x in 1:lenm]
end
