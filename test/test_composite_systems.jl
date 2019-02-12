using Test, QSimulator
using LinearAlgebra: I

@testset "CompositeQSystem hamiltonian" begin
    # check that we correctly embed the subystem Hamiltonians
    r1 = Resonator("r1", 3, ResonatorSpec(4.0))
    r2 = Resonator("r2", 2, ResonatorSpec(5.0))
    r3 = Resonator("r3", 2, ResonatorSpec(6.0))
    qs = [r1, r2, r3]
    cqs = CompositeQSystem(qs)
    for q in qs; add_hamiltonian!(cqs, q); end
    m1 = hamiltonian(r1) ⊗ Matrix{Int}(I, 2, 2) ⊗ Matrix{Int}(I, 2, 2)
    m2 = Matrix{Int}(I, 3, 3) ⊗ hamiltonian(r2) ⊗ Matrix{Int}(I, 2, 2)
    m3 = Matrix{Int}(I, 3, 3) ⊗ Matrix{Int}(I, 2, 2) ⊗ hamiltonian(r3)
    @test hamiltonian(cqs) == m1 + m2 + m3
end

@testset "embed" begin
    CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]
    SWAP = [1 0 0 0; 0 0 1 0;0 1 0 0; 0 0 0 1]

    # trivial embedding should return itself
    @test QSimulator.embed(CNOT, [1,2], [2,2]) == CNOT

    # we can use embed to swap control target
    @test QSimulator.embed(CNOT, [2,1], [2,2]) == SWAP*CNOT*SWAP

    # use embed to create a CNOT13 in a three qubit system
    Id = Matrix{Int}(I, 2, 2)
    @test QSimulator.embed(CNOT, [1,3], [2,2,2]) == (Id⊗SWAP) * (CNOT⊗Id) * (Id⊗SWAP)

    # check an arbitrary matrix with unequal subsystem dimensions
    @test QSimulator.embed([1 2;3 4], [2], [3, 2, 4]) == Matrix{Int}(I, 3, 3) ⊗ [1 2; 3 4] ⊗ Matrix{Int}(I, 4, 4)
end
