using Iterators: product

export photons_to_index, index_to_photons, photons_to_state, index_to_state, tensor_product_states
export unique_tol

"""
    tensor_product_states(dims::Vector{Int})

Enumerate all basis states in a tensor product space in the canonical order
used by the Kronecker product.

## args
* `dims`: an array with the dimension of each mode.

## returns
An array of arrays where each array indicates the number of photons
in each mode.

## example
`tensor_product_states([2,2]) == [[0,0], [0,1], [1,0], [1,1]]`.
"""
function tensor_product_states(dims::Vector{Int})
   return [reverse(collect(prod_inds).-1) for prod_inds in product([1:dim for dim in reverse(dims)]...)]
end

"""
    photons_to_index(photons::Vector{Int}, dims::Vector{Int})

Convert number of photons in each subsystem of a tensor product space to an index in the full space.

## args
* `photons`: a vector with the number of photons in each subsystem.
* `dims`: a vector of dimensions for each subsystem.

## returns
The index of the state according to the indexing rules of the kronecker product,
starting at 1.

## example
Consider a tensor product of two qubits so that `dims = [2,2]`. Then states
`[0,0], [0,1], [1,0], [1,1]` are numbered `1, 2, 3, 4` so that
`photons_to_index([1,0], [2,2]) == 3`.
"""
function photons_to_index(photons::Vector{Int}, dims::Vector{Int})
    @assert length(photons) == length(dims)
    for (photon_num, dim) in zip(photons, dims)
        if !(0 <= photon_num < dim)
            return nothing
        end
    end
    return LinearIndices(tuple(reverse(dims)...))[reverse(photons.+1)...]
end

"""
    index_to_photons(index::Int, dims::Vector{Int})

Convert an index in a tensor product space to the number of photons in each subsystem.

## args
* `index`: the index of the state in the tensor product space.
* `dims`: a vector of dimensions for each subsystem.

## returns
The an array with the number of photons in each subsystem.

## example
Consider a tensor product of two qubits so that `dims = [2,2]`. Then states
`[0,0], [0,1], [1,0], [1,1]` are numbered `1, 2, 3, 4` so that
`index_to_photons(3, [2,2]) == [1,0]`.
"""
function index_to_photons(index::Int, dims::Vector{Int})
    if !(1 <= index <= prod(dims))
        return nothing
    end
    return collect(reverse(Tuple(CartesianIndices(tuple(reverse(dims)...))[index])).-1)
end

"""
    photons_to_state(photons::Vector{Int}, dims::Vector{Int})

Create the state corresponding to a given number of photons in each mode.

## args
* `photons`: an array with the number of photons in each mode.
* `dims`: a array with the dimension of each mode.

## returns
The standard basis state corresponding to the given number of photons
in each mode.
"""
function photons_to_state(photons::Vector{Int}, dims::Vector{Int})
    index = photons_to_index(photons, dims)
    return index_to_state(index, prod(dims))
end

"""
    index_to_state(index::Int, dim::Int)

Create a standard basis vector.

## args
* `index`: the index of the standard basis vector.
* `dim`: the dimension of the vector space.

## returns
The standard basis vector.
"""
function index_to_state(index::Int, dim::Int)
    ans = zeros(ComplexF64, dim)
    ans[index] = 1.0
    return ans
end

"""
    unique_tol(ts::Vector{<:Real}, dt::Real)

Compute an array `a` of values where each is within `dt/2` of some value in `ts`
and no two differ by less than `dt`. Further return a vector of indices `inds`
such that `a[inds]` is an approximation of `ts` up to `dt`.

## args
* `ts`: an array of numbers.
* `dt`: a small number within which errors are unimportant.

## returns
An array of the unique values and an array of indices.
"""
function unique_tol(ts::Vector{<:Real}, dt::Real)
    vals = round.(Int, ts ./ dt) .* dt
    unique_vals = unique(vals)
    unique_inds = [findfirst(unique_vals .== v) for v in vals]
    return unique_vals, unique_inds
end
