module QSimulator

export ⊗

const ⊗ = kron

include("systems.jl")
include("composite_systems.jl")
include("time_evolution.jl")
end
