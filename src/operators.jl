using LinearAlgebra: diagm

export raising, lowering, number, X, Y, X_Y,
       decay, dephasing, dipole_drive, parametric_drive

######################################################
# Primitives
######################################################

raising(q::QSystem, ϕ::Real=0.0) = diagm(-1 => sqrt.(1:(dimension(q)-1))) * exp(1im*2π*ϕ)
lowering(q::QSystem, ϕ::Real=0.0) = diagm(1 => sqrt.(1:(dimension(q)-1))) * exp(-1im*2π*ϕ)
number(q::QSystem) = diagm(0 => collect(ComplexF64, 0:dimension(q)-1))

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
    dipole_drive(qs::QSystem, drive::Function, rotation_rate::Real=0.0)

Given some function of time, return a function applying a time dependent
dipole Hamiltonian.

## args
* `qs`: a QSystem.
* `drive`: a function of time returning a real or complex value. The real
    part couples to X and the imaginary part couples to Y.
* `rotation_rate`: the rotation rate of a rotating frame.

## returns
A function of time.
"""
function dipole_drive(qs::QSystem, drive::Function, rotation_rate::Real=0.0)
    ham(t) = drive(t) * X(qs, rotation_rate * t)
    return ham
end

"""
    parametric_drive(qs::QSystem, drive::Function)

Given some function of time, return a function applying a
time dependent Hamiltonian.

## args
* `qs`: a QSystem with a method of `hamiltonian` accepting a function of time.
* `drive`: a function of time returning a real value.

## returns
A function of time.
"""
function parametric_drive(qs::QSystem, drive::Function)
    ham(t) = hamiltonian(qs, drive(t))
    return ham
end
