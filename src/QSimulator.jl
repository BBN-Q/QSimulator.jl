VERSION >= v"0.4.0-dev+6521" && __precompile__()

module QSimulator

using Interpolations,
      Iterators

import Base: getindex, +

include("base-types.jl")
include("systems.jl")
include("controls.jl")
include("interactions.jl")
include("composite-systems.jl")
include("dissipation.jl")
include("evolution.jl")
include("read_APS_file.jl")

end
