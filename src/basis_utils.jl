using Base.Iterators: product
import Base: vec, getindex

export TensorProductBasis, TensorProductBasisState, basis_states

"""
A simple covenience wrapper for a tuple of subsystem dimensions for a tensor product space.
"""
struct TensorProductBasis
    dims::Tuple{Vararg{Int}}
end

# helper contructor to convert from dimensions given in AbstractVector form
TensorProductBasis(dims::AbstractVector{Int}) = TensorProductBasis(Tuple(dims))

abstract type BasisState end

"""
A basis element of a tensor product space.
"""
struct TensorProductBasisState <: BasisState
    basis::TensorProductBasis
    states::Tuple{Vararg{Int}}

    function TensorProductBasisState(b::TensorProductBasis, states::Tuple{Vararg{Int}})
        # check that given subsystem basis states are compatible with TensorProductBasis dimensions
        @assert length(states) == length(b.dims)
        @assert all(states .>= 0) && all(states .< b.dims)
        return new(b, states)
    end
end

# helper contructor to convert from dimensions given in AbstractVector form
TensorProductBasisState(b::TensorProductBasis, states::AbstractVector{Int}) = TensorProductBasisState(b, Tuple(states))

"""
    vec(state::TensorProductBasisState)

Create the complex basis vector correspoding to a given TensorProductBasisState.

## args
* `state`: a TensorProductBasisState to convert to a complex state vector.

## returns
The complex basis state vector.
"""
function vec(state::TensorProductBasisState)
    state_vector = zeros(ComplexF64, prod(state.basis.dims))
    state_vector[index(state)] = 1.0
    return state_vector
end


"""
    basis_states(b::TensorProductBasis)

Enumerate all basis states in a tensor product space in the canonical order used by the Kronecker
product.

## args
* `b`: a TensorProductBasis.

## returns
A vector of TensorProductBasisState.

## example
`[bs.states for bs in basis_states(TensorProductBasis((2,2)))] == [(0,0), (0,1), (1,0), (1,1)]]`.
"""
function basis_states(b::TensorProductBasis)
    return [TensorProductBasisState(b, reverse(x .- 1))
            for x in vec(collect(product([1:dim for dim in reverse(b.dims)]...)))]
end


"""
    index(bs::TensorProductBasisState)

Convert basis state of a tensor product space to an index into an enumerate of all basis states from
the canonical order used by the Kronecker product.

## args
* `bs`: a TensorProductBasisState.

## returns
The linear index of the state into the canonical tensor product basis states.

## example
Consider a tensor product of two qubits `b = TensorProductBasis((2,2))`. Then states `(0,0), (0,1),
(1,0), (1,1)` are numbered `1, 2, 3, 4` so that `index(TensorProductBasisState(b, (1,0)) == 3`.
"""
function index(bs::TensorProductBasisState)
    return LinearIndices(tuple(reverse(bs.basis.dims)...))[reverse(bs.states.+1)...]
end


"""
    getindex(b::TensorProductBasis, i::Int)

Convert an index into a tensor product space to the basis state. Note that
this function can be called with the notation `b[i]`.

## args
* `b`: the TensorProductBasis into which to index.
* `i`: the index of the basis state in the tensor product space.

## returns
The indexed TensorProductBasisState.

## example
Consider a tensor product of two qubits so that `b = TensorProductBasis((2,2))`. Then the basis
states `(0,0), (0,1), (1,0), (1,1)` are numbered `1, 2, 3, 4` so that
`b[3] == TensorProductBasisState(b, (1,0))`.
"""
function getindex(b::TensorProductBasis, i::Int)
    states = reverse(Tuple(CartesianIndices(tuple(reverse(b.dims)...))[i])).-1
    return TensorProductBasisState(b, states)
end
