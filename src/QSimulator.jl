module QSimulator

using Grid,
      NumericExtensions,
      Iterators

include("base-types.jl")

include("systems.jl")

include("controls.jl")

include("interactions.jl")

include("composite-systems.jl")

#include("dissipation.jl")

include("evolution.jl")

include("read_APS_file.jl")

include("paulis.jl")

end