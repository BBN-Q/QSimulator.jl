module QSimulator

using Grid,
      NumericExtensions,
      QIP

import QIP.hamiltonian

include("base-types.jl")

include("systems.jl")

include("controls.jl")

include("interactions.jl")

include("composite-systems.jl")

include("dissipation.jl")

# [todo] - Do away with the early evolution routines
include("evolution.jl")

end