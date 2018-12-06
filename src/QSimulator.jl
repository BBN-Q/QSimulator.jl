module QSimulator

using Interpolations

import Base: getindex
import Base.+
using LinearAlgebra
using SparseArrays


eye(dim::Integer) = Matrix{Float64}(LinearAlgebra.I, dim, dim)
eye(::Type{T}, dim::Integer) where T = Matrix{T}(LinearAlgebra.I, dim, dim)
speye(dim::Integer) = sparse(one(Float64) * LinearAlgebra.I, dim, dim)
speye(::Type{T}, dim::Integer) where T = sparse(one(T) * LinearAlgebra.I, dim, dim)

# `findin(a, b)` is deprecated, use `findall((in)(b), a)`
findin(a, b) = findall((in)(b), a) # :)

include("base-types.jl")
include("systems.jl")
include("controls.jl")
include("interactions.jl")
include("composite-systems.jl")
include("dissipation.jl")
include("evolution.jl")
include("read_APS_file.jl")

end
