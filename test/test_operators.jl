using Test, QSimulator

using LinearAlgebra: I
@testset "Primitives" begin
    r = Resonator("test", 3, QSimulator.ResonatorSpec(5.5))

    # test some trivial cases of the ladder operators
    @test raising(r) == [0 0 0;1 0 0;0 √2 0]
    @test lowering(r) == [0 1 0;0 0 √2;0 0 0]

    r = Resonator("test", 7, QSimulator.ResonatorSpec(5.5))

    # test harmonic oscillator operator identity a†a = N
    @test raising(r) * lowering(r) ≈ number(r)

    # test X = a + a†
    @test X(r) ≈ raising(r) + lowering(r)

    # test phase rotation of 0.25 on X → Y
    @test X(r, 0.25) ≈ Y(r)

    # test commutator relation [a, a†] = 1; because of truncation it doesn't hold for last row/columns
    commutator = lowering(r)*raising(r)- raising(r)*lowering(r)
    @test commutator[1:end-1,1:end-1] ≈ Matrix{eltype(commutator)}(I, dimension(r)-1, dimension(r)-1)

end
