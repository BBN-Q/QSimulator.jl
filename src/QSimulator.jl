module QSimulator

export ⊗

const ⊗ = kron

include("systems.jl")
include("composite_systems.jl")
include("hamiltonians.jl")
include("parametric_hamiltonians.jl")
include("time_evolution.jl")
<<<<<<< HEAD
include("hamiltonians.jl")
=======
>>>>>>> e26ba3fc74b52e8863e9c9c562f84a8378e695cd
end
