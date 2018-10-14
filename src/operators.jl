using LinearAlgebra: diagm

export raising, lowering, number, X, Y, X_Y, rotating_operator,
       decay, dephasing, dipole_drive, flux_drive

######################################################
# Primitives
######################################################

raising(q::QSystem, ϕ::Real=0.0) = diagm(-1 => sqrt.(1:(dim(q)-1))) * exp(1im*2π*ϕ)
lowering(q::QSystem, ϕ::Real=0.0) = diagm(1 => sqrt.(1:(dim(q)-1))) * exp(-1im*2π*ϕ)
number(q::QSystem) = diagm(0 => collect(ComplexF64, 0:dim(q)-1))

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
function decay(qs::QSystem, γ::Real)
    return sqrt(γ) * lowering(qs)
end

"""
    dephasing(qs::QSystem, γ::Real)

Dephasing for a QSystem.

## args
* `qs`: a QSystem.
* `γ`: a decay rate in frequency units. Note Tϕ = 1/(2πγ).

## returns
The lindblad operator for decay.
"""
function dephasing(qs::QSystem, γ::Real)
    return sqrt(2γ) * number(qs)
end

"""
    dipole_drive(qs::QSystem, drive::Function)

Given some function of time, return a function applying a time dependent
dipole Hamiltonian.

## args
* `qs`: a QSystem.
* `drive`: a function of time returning a real or complex value. The real
    part couples to X and the imaginary part couples to Y.

## returns
A function of time.
"""
function dipole_drive(qs::QSystem, drive::Function)
    x_ham = X(qs)
    y_ham = Y(qs)
    function ham(t)
        pulse = drive(t)
        return real(pulse) * x_ham + imag(pulse) * y_ham
    end
    return ham
end


"""
    flux_drive(qs::QSystem, drive::Function)

Given some function of time, return a function applying a
time dependent Hamiltonian.

## args
* `qs`: a QSystem with a method of `hamiltonian` accepting a function of time.
* `drive`: a function of time returning a real value.

## returns
A function of time.
"""
function flux_drive(qs::QSystem, drive::Function)
    ham(t) = hamiltonian(qs, drive(t))
    return ham
end

######################################################
# Backwards compatibility
######################################################

export microwave_drive, dipole, flip_flop, XY, rotating_flip_flop

function microwave_drive(q::QSystem, drive::Function)
    @warn "Deprecation warning: microwave_drive."
    return dipole_drive(q, drive)
end

function dipole(a::QSystem, b::QSystem)
    @warn "Deprecation warning: dipole."
    return X([a, b])
end

function XY(a::QSystem, b::QSystem; ϕ::Real=0.0)
    @warn "Deprecation warning: XY."
    return .5 * X_Y([a, b], [ϕ, 0.0])
end

function flip_flop(a::QSystem, b::QSystem; ϕ::Real=0.0)
    @warn "Deprecation warning: flip_flop."
    return .5 * X_Y([a, b], [ϕ, 0.0])
end

function rotating_flip_flop(a::QSystem, b::QSystem, strength::Real, freq::Real)
    @warn "Deprecation warning: rotating_flip_flop."
    op(t) = .5 * strength * X_Y([a, b], [freq, 0.0] * t)
    return op
end
