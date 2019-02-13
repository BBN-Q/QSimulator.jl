using LinearAlgebra: I

export CompositeQSystem, add_hamiltonian!, add_lindblad!

# tensor products of quantum systems
mutable struct CompositeQSystem
    subsystems::Vector
    fixed_Hs::Vector{Tuple} # tuple of Matrix and expansion indices
    parametric_Hs::Vector{Tuple} # tuple of Functions and expansion indices
    fixed_Ls::Vector{Tuple} # tuple of Matrix and expansion indices for collapse operators
    parametric_Ls::Vector{Tuple} # tuple of Matrix, expansion indices, and time-dependent function for collapse operators
    dimension::Int
end

CompositeQSystem(qs) = CompositeQSystem(qs, [], [], [], [], prod(dimension(q) for q in qs))

# helper functions for CompositeQSystems
dimension(cqs::CompositeQSystem) = cqs.dimension
dim(cqs::CompositeQSystem) = dimension(cqs)

find_indices(cqs::CompositeQSystem, s::QSystem) = findall([sub == s for sub in cqs.subsystems])
find_indices(cqs::CompositeQSystem, s_label::AbstractString) = findall([label(sub) == s_label for sub in cqs.subsystems])
function find_indices(cqs::CompositeQSystem, s::Union{Vector{<:QSystem}, Vector{<:AbstractString}})
    return vcat([find_indices(cqs, _s) for _s in s]...)
end

""" Add a subsystem static Hamiltonian matrix to a CompositeQSystem """
function add_hamiltonian!(cqs::CompositeQSystem, ham::AbstractMatrix{T}, acting_on::Union{Q, Array{Q}}) where {T<:Number, Q<:QSystem}
    idxs = embed_indices(cqs, acting_on)
    push!(cqs.fixed_Hs, (ham, idxs))
end

""" Add a subystem Hamiltonian to a CompositeQSystem """
add_hamiltonian!(cqs::CompositeQSystem, qs::QSystem) = add_hamiltonian!(cqs, hamiltonian(qs), qs)
add_hamiltonian!(cqs::CompositeQSystem, qs::QSystem, ϕ::Real) = add_hamiltonian!(cqs, hamiltonian(qs, ϕ), qs)

# TODO how to do this dispatch vs adding a fixed Hamiltonian - Jameson says not to do this https://discourse.julialang.org/t/functions-and-callable-methods/2983/3
""" Add a time parameterized subsystem Hamiltonian to a CompositeQSystem """
function add_hamiltonian!(cqs::CompositeQSystem, ham::Function, acting_on::Union{Q, Array{Q}}) where Q<:QSystem
    idxs = embed_indices(find_indices(cqs, acting_on), [dimension(s) for s in cqs.subsystems])
    push!(cqs.parametric_Hs, (ham, idxs))
end

""" Add the parametric hamiltonians from a CompositeQSystem to a Hamiltonian matrix, in place. """
function add_parametric_hamiltonians!(ham::Matrix{<:Number}, cqs::CompositeQSystem, t::Real)
    for (ham_func, idxs) = cqs.parametric_Hs
        embed_add!(ham, ham_func(t), idxs)
    end
end

""" Add a subystem collapse static operator matrix to a CompositeQSystem """
function add_lindblad!(cqs::CompositeQSystem, lind_op::AbstractMatrix{T}, acting_on::Union{Q, Array{Q}}) where {T<:Number, Q<:QSystem}
    idxs = embed_indices(cqs, acting_on)
    push!(cqs.fixed_Ls, (lind_op, idxs))
end

""" Add a subystem time dependent collapse operator matrix to a CompositeQSystem """
function add_lindblad!(cqs::CompositeQSystem, lind_op::Function, acting_on::Union{Q, Array{Q}}) where Q<:QSystem
    idxs = embed_indices(cqs, acting_on)
    push!(cqs.parametric_Ls, (lind_op, idxs))
end

""" In place addition of an operator embedded into a larger Hilbert space given a set of expansion indices"""
function embed_add!(op, added_op, expand_idxs)
    # TODO explore performance of @inbounds
    for ct = 1:length(expand_idxs)
        for idx = expand_idxs[ct]
            op[idx] += added_op[ct]
        end
    end
end

"""
    embed(m::Matrix, acting_on::Vector, dims::Vector)

Embed a subsystem operator `m` acting on subsystem indices `acting_on` into a larger tensor product
space with subsystem dimensions `dims`. This could be implemented by starting with a tensor product
ordering of `m` ⊗ I and then swapping the subsystems of `m` into the desired order. E.g. in a three
qubit system CNOT13 can be created with Id⊗SWAP * CNOT⊗Id * Id⊗SWAP.

In practice we use the fact that a tensor product structure is a statement about the strides bewteen
various subsystem elements and using reshaping the array to match the subsystem dimensions to expose
those strides and make the swaps a `permutedims`.
"""
function embed(m::Matrix, acting_on::Vector, dims::Vector)

    @assert size(m, 1) == prod(dims[acting_on]) "Oops! Dimensions of operator do not match dims argument."

    # create the large matrix by tensoring on an identity term to the operator`m`
    d_rest = prod(dims) ÷ size(m,1)
    M = kron(m, Matrix{eltype(m)}(I, d_rest, d_rest))

    # Reshape into a multi-dimensional array with the dimensions given by swapped (i.e. `m` operator
    # first) subsystem dimensions. Because the of the way tensor product indices work (last
    # subsystem is fastest) we need to reverse the order of the dimensions.
    # e.g
    # julia> reshape(["a"; "b"] ⊗ ["a"; "b"; "c"] ⊗ ["a"; "b"; "c"; "d"], (4,3,2))
    # 4×3×2 Array{String,3}:
    # [:, :, 1] =
    #  "aaa"  "aba"  "aca"
    #  "aab"  "abb"  "acb"
    #  "aac"  "abc"  "acc"
    #  "aad"  "abd"  "acd"
    #
    # [:, :, 2] =
    #  "baa"  "bba"  "bca"
    #  "bab"  "bbb"  "bcb"
    #  "bac"  "bbc"  "bcc"
    #  "bad"  "bbd"  "bcd"

    # For a matrix we repeat for rows then columns
    l = length(dims)
    identity_idxs = filter(x -> !(x in acting_on), 1:l)
    swapped_dims = reverse([dims[acting_on]; dims[identity_idxs]])
    M = reshape(M, tuple([swapped_dims; swapped_dims]...))

    # now we want to permute the dimensions to restore the correct subsystem ordering start from the
    # permutation to get the current ordering with the `m` operator first; then invert and again
    # handle the way tensor product indices work e.g. if we wanted the permutation [2 1 3] (i.e.
    # swap 1st and 2nd subsystem) then we actually need to swap the 2nd and 3rd dimensions
    reverse_perm = (l+1) .- reverse(invperm([acting_on; identity_idxs]))
    M = permutedims(M, [reverse_perm; reverse_perm .+ l])

    # reshape back
    return reshape(M, prod(dims), prod(dims))
end

"""
    embed_indices(acting_on::Vector, dims::Vector)

Return the linear indices that map a subystem operation acting on subystem indices `acting_on` into
a larger tensor product space with subsystem dimensions `dims`.
"""
function embed_indices(acting_on::Vector, dims::Vector)
    # Strategy:
    # 1. create a subystem matrix with elements equal to the linear index e.g. [1 2; 3 4]
    dim_acting_on = prod(dims[acting_on])
    m = reshape(collect(1:dim_acting_on^2), (dim_acting_on,dim_acting_on))
    # 2. embed the matrix
    M = embed(m, acting_on, dims)
    # 3. find where the linear indices got embedded
    f(A) = (LinearIndices(A))[findall(A)]
    return [f(M .== x) for x in 1:dim_acting_on^2]
end

embed_indices(cqs::CompositeQSystem, acting_on) = embed_indices(find_indices(cqs, acting_on), [dimension(s) for s in cqs.subsystems])

""" Calculate the static part of the Hamiltonian of a CompositeQSystem """
function hamiltonian(cqs::CompositeQSystem)
    ham = zeros(ComplexF64, (dimension(cqs), dimension(cqs)))
    for (new_ham, idxs) in cqs.fixed_Hs
        embed_add!(ham, new_ham, idxs)
    end
    return ham
end
