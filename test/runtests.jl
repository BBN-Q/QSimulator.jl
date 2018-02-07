using QSimulator

using Base.Test

# test harmonic oscillator operator identity a†a = N
r = Resonator("test", 5.5, 3)
@test raising(r) * lowering(r) == number(r)
