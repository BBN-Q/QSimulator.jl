using Test, QSimulator

@testset "basis state utilities" begin
    # check that constructors from AbstractVector work
    basis = TensorProductBasis((2,2))
    @test basis == TensorProductBasis([2,2])
    @test TensorProductBasisState(basis, (0,0)) == TensorProductBasisState(basis, [0,0])

    # simple tensor product expansion
    basis = TensorProductBasis((2,2))
    @test [bs.states for bs in basis_states(basis)] == [(0,0), (0,1), (1,0), (1,1)]
    basis = TensorProductBasis((2,3))
    @test [bs.states for bs in basis_states(basis)] == [(0,0), (0,1), (0,2), (1,0), (1,1), (1,2)]

    # go from state -> index
    basis = TensorProductBasis((2,3,4,5))
    all_states = basis_states(basis)
    @test all(QSimulator.index(state) == i for (i, state) in enumerate(all_states))

    # and from index to state
    @test all(basis[i] == state for (i,state) in enumerate(all_states))

    # check that TensorProductBasisState construction catches errors
    basis = TensorProductBasis((2,2))
    # wrong number of subsystems
    @test_throws AssertionError TensorProductBasisState(basis, (0,0,0))
    # state index exceeds dimension
    @test_throws AssertionError TensorProductBasisState(basis, (2,0))
end
