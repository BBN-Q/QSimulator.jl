using Test, QSimulator
using LinearAlgebra: I

@testset "failing test" begin
    r1 = Resonator("r1", 3, ResonatorSpec(4.0))
    r2 = Resonator("r2", 2, ResonatorSpec(5.0))
    r3 = Resonator("r3", 2, ResonatorSpec(6.0))
    qs = [r1, r2, r3]
    cqs = CompositeQSystem(qs)
    for q in qs; add_hamiltonian!(cqs, q); end
    m1 = hamiltonian(r1) ⊗ Matrix(I, 2, 2) ⊗ Matrix(I, 2, 2)
    m2 = Matrix(I, 3, 3) ⊗ hamiltonian(r2) ⊗ Matrix(I, 2, 2)
    m3 = Matrix(I, 3, 3) ⊗ Matrix(I, 2, 2) ⊗ hamiltonian(r3)
    @test hamiltonian(cqs) == m1 + m2 + m3
end

@testset "passes if all are dimension 3" begin
    r1 = Resonator("r1", 3, ResonatorSpec(4.0))
    r2 = Resonator("r2", 3, ResonatorSpec(5.0))
    r3 = Resonator("r3", 3, ResonatorSpec(6.0))
    qs = [r1, r2, r3]
    cqs = CompositeQSystem(qs)
    for q in qs; add_hamiltonian!(cqs, q); end
    m1 = hamiltonian(r1) ⊗ Matrix(I, 3, 3) ⊗ Matrix(I, 3, 3)
    m2 = Matrix(I, 3, 3) ⊗ hamiltonian(r2) ⊗ Matrix(I, 3, 3)
    m3 = Matrix(I, 3, 3) ⊗ Matrix(I, 3, 3) ⊗ hamiltonian(r3)
    @test hamiltonian(cqs) == m1 + m2 + m3
end

@testset "or if there is only two resonators" begin
    r1 = Resonator("r1", 3, ResonatorSpec(4.0))
    r2 = Resonator("r2", 2, ResonatorSpec(5.0))
    qs = [r1, r2]
    cqs = CompositeQSystem(qs)
    for q in qs; add_hamiltonian!(cqs, q); end
    m1 = hamiltonian(r1) ⊗ Matrix(I, 2, 2)
    m2 = Matrix(I, 3, 3) ⊗ hamiltonian(r2)
    @test hamiltonian(cqs) == m1 + m2
end
