module QSimulator

export ⊗

const ⊗ = kron

include("systems.jl")
include("composite_systems.jl")
include("nicometrics.jl")
using .Nicometrics

include("hamiltonians.jl")
include("parametric_hamiltonians.jl")
include("time_evolution.jl")
end
