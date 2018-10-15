using QSimulator
using LsqFit

using Test

# test harmonic oscillator operator identity a†a = N
r = Resonator("test", 3, QSimulator.ResonatorSpec(5.5))
@test isapprox(raising(r) * lowering(r), number(r))

include("test_evolution.jl")
include("test_fit.jl")
include("test_fourier.jl")
