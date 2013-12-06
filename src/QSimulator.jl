module QSimulator

using Grid,
      NumericExtensions

include("base-types.jl")

include("systems.jl")

include("controls.jl")

include("interactions.jl")

include("composite-systems.jl")

#include("dissipation.jl")

include("evolution.jl")

include("spectra.jl")

end