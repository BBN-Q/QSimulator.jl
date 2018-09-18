export raising, lowering, number, X, Y, X_Y, rotating_operator,
       decay, dephasing

######################################################
# Primitives
######################################################

raising(q::QSystem, ϕ::Real=0.0) = diagm(sqrt.(1:(dim(q)-1)), -1) * exp(1im*2π*ϕ)
lowering(q::QSystem, ϕ::Real=0.0) = diagm(sqrt.(1:(dim(q)-1)), 1) * exp(-1im*2π*ϕ)
number(q::QSystem) = diagm(collect(Complex128, 0:dim(q)-1))

X(q::QSystem, ϕ::Real=0.0) = raising(q, ϕ) + lowering(q, ϕ)
X(qs::Vector{<:QSystem}, ϕs::Vector{<:Real}) = reduce(⊗, [X(q, ϕ) for (q, ϕ) in zip(qs, ϕs)])
X(qs::Vector{<:QSystem}) = reduce(⊗, [X(q) for q in qs])

Y(q::QSystem, ϕ::Real=0.0) = 1im*(raising(q, ϕ) - lowering(q, ϕ))
Y(qs::Vector{<:QSystem}, ϕs::Vector{<:Real}) = reduce(⊗, [Y(q, ϕ) for (q, ϕ) in zip(qs, ϕs)])
Y(qs::Vector{<:QSystem}) = reduce(⊗, [Y(q) for q in qs])

X_Y(qs::Vector{<:QSystem}, ϕs::Vector{<:Real}) = X(qs, ϕs) + Y(qs, ϕs)
X_Y(qs::Vector{<:QSystem}) = X(qs) + Y(qs)


"""
    decay(qs::QSystem, γ:Real)

T1 decay for a QSystem.

## args
* `qs`: a QSystem.
* `γ`: a decay rate in frequency units. Note T1 = 1/(2πγ).

## returns
The lindblad operator for decay.
"""
decay(qs::QSystem, γ::Real) = sqrt(γ) * lowering(qs)

"""
    dephasing(qs::QSystem, γ::Real)

Dephasing for a QSystem.

## args
* `qs`: a QSystem.
* `γ`: a decay rate in frequency units. Note Tϕ = 1/(2πγ).

## returns
The lindblad operator for decay.
"""
dephasing(qs::QSystem, γ::Real) = sqrt(2γ) * number(qs)


######################################################
# Backwards compatibility
######################################################

export dipole, flip_flop, XY, rotating_flip_flop

function dipole(a::QSystem, b::QSystem)
    warn("Deprecation warning: dipole.")
    return X([a, b])
end

function XY(a::QSystem, b::QSystem; ϕ::Real=0.0)
    warn("Deprecation warning: XY.")
    return .5 * X_Y([a, b], [ϕ, 0.0])
end

function flip_flop(a::QSystem, b::QSystem; ϕ::Real=0.0)
    warn("Deprecation warning: flip_flop.")
    return .5 * X_Y([a, b], [ϕ, 0.0])
end

function create(q::QSystem)
    warn("Deprecation warning: create.")
    return raising(q)
end

function destroy(q::QSystem)
    warn("Deprecation warning: destroy.")
    return lowering(q)
end
