module QSimulator

export ⊗

const ⊗ = kron

include("systems.jl")
include("composite_systems.jl")
include("parametric_hamiltonians.jl")
include("time_evolution.jl")
include("hamiltonians.jl")
end
