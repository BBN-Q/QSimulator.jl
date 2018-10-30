using Test, QSimulator

@testset "basis state utilities" begin
    # simple tensor product expansion
    basis = QSimulator.TensorProductBasis((2,2))
    @test [bs.states for bs in QSimulator.basis_states(basis)] == [(0,0), (0,1), (1,0), (1,1)]
    basis = QSimulator.TensorProductBasis((2,3))
    @test [bs.states for bs in QSimulator.basis_states(basis)] == [(0,0), (0,1), (0,2), (1,0), (1,1), (1,2)]

    # go from state -> index
    basis = QSimulator.TensorProductBasis((2,3,4,5))
    all_states = QSimulator.basis_states(basis)
    @test all(QSimulator.index(state) == i for (i,state) in enumerate(all_states))

    # and from index to state
    @test all(basis[i] == state for (i,state) in enumerate(all_states))

    # check that TensorProductBasisState construction catches errors
    basis = QSimulator.TensorProductBasis((2,2))
    # wrong number of subsystems
    @test_throws AssertionError QSimulator.TensorProductBasisState(basis, (0,0,0))
    # state index exceeds dimension
    @test_throws AssertionError QSimulator.TensorProductBasisState(basis, (2,0))
end
