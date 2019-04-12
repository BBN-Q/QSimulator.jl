using QSimulator
using LinearAlgebra

include("unitary-tests.jl")
include("dissipation-tests.jl")

include("speed_test.jl")
speed_test(5)
sim_setup(3, 10)
