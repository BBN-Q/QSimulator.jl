using Base.Iterators: product

import Base: vec, getindex


"""
A simple covenience wrapper for a tuple of subsystem dimensions for a tensor product space.
"""
struct TensorProductBasis
    dims::Tuple{Vararg{Int}}
end

# helper contructor to convert from dimensions given in Vector form
TensorProductBasis(dims::Vector{Int}) = TensorProductBasis(Tuple(dims))

abstract type BasisState end;

"""
A basis element of a tensor product space
"""
struct TensorProductBasisState <: BasisState
    basis::TensorProductBasis
    states::Tuple{Vararg{Int}}
end

# helper contructor to convert from dimensions given in Vector form
TensorProductBasis(b::TensorProductBasis, states::Vector{Int}) = TensorProductBasis(b, Tuple(states))

"""
    vec(bs::TensorProductBasisState)

Create the complex basis vector correspoding to a given TensorProductBasisState basis state

## args
* `bs`: a TensorProductBasisState to conver to a complex state vector

## returns
The complex basis state vector
"""
function vec(state::TensorProductBasisState)
    state_vector = zeros(ComplexF64, prod(state.basis.dims))
    state_vector[index(state)] = 1.0
    state_vector
end


"""
    basis_states(b::TensorProductBasis)

Enumerate all basis states in a tensor product space in the canonical order used by the Kronecker
product.

## args
* `b`: a TensorProductBasis

## returns
A vector of TensorProductBasisState

## example
`[bs.states for bs in basis_states(TensorProductBasis((2,2))] == [(0,0), (0,1), (1,0), (1,1)]]`.
"""
function basis_states(b::TensorProductBasis)
    return [
        TensorProductBasisState(b, reverse(x .- 1))
            for x in vec(collect(product([1:dim for dim in reverse(b.dims)]...)))
            ]
end


"""
    index(bs::TensorProductBasisState)

Convert basis state of a tensor product space to an index into an enumerate of all basis states from
the canonical order used by the Kronecker product.

## args
* `bs`: a TensorProductBasisState

## returns
The 1-index of the state into the cannonical tensor product basis states.

## example
Consider a tensor product of two qubits `b = TensorProductBasis((2,2))`. Then states `(0,0), (0,1),
(1,0), (1,1)` are numbered `1, 2, 3, 4` so that `index(TensorProductBasisState(b, (1,0)) == 3`.
"""
function index(bs::TensorProductBasisState)
    return LinearIndices(tuple(reverse(bs.basis.dims)...))[reverse(bs.states.+1)...]
end


"""
    getindex(b::TensorProductBasis, i)

Convert an index into a tensor product space to the basis state

## args
* `b` : the TensorProductBasis to index into
* `i`: the index of the basis state in the tensor product space.

## returns
The indexed TensorProductBasisState.

## example
Consider a tensor product of two qubits so that ``b = TensorProductBasis((2,2))``. Then the basis
states `(0,0), (0,1), (1,0), (1,1)` are numbered `1, 2, 3, 4` so that
`b[3] == TensorProductBasisState(b, (1,0))`.
"""
function getindex(b::TensorProductBasis, i)
    states = reverse(Tuple(CartesianIndices(tuple(reverse(b.dims)...))[i])).-1
    TensorProductBasisState(b, states)
end
