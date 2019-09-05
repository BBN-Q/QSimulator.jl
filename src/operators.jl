using LinearAlgebra: diagm, lmul!

export raising, lowering, number, X, Y, X_Y, XY, flip_flop,
       decay, dephasing, dipole_drive, parametric_drive, rwa_dipole

######################################################
# Primitives
######################################################

# ladder operators
"""
    raising(q::QSystem, ϕ::Real=0.0)

raising/creation/a ladder operator. Given an eignenstate of the number operator `|n⟩` on a
`QSystem q`, `raising(q)|n⟩ = √(n+1)|n+1⟩`. Optionally an additional phase exp(2πiϕ) is applied.

# Examples
```jldoctest
julia> q = DuffingTransmon("q", 3, DuffingSpec(5, -0.2))
DuffingTransmon("q", 3, DuffingSpec(5, -0.2))

julia> raising(q)
3×3 Array{Float64,2}:
 0.0  0.0      0.0
 1.0  0.0      0.0
 0.0  1.41421  0.0
```
"""
raising(q::QSystem, ϕ::Real=0.0) = raising(dimension(q), ϕ)

# custom diagm because it's twice as fast for small matrices
function raising(dim::Integer,  ϕ::Real=0.0)
    m = zeros(typeof(complex(ϕ)), dim, dim)
    fac = exp(1im * 2π * ϕ)
    for i in 1:(dim -1 )
        m[i + 1, i] = sqrt(i) * fac
    end
    return m
end

"""
    lowering(q::QSystem, ϕ::Real=0.0, factor::Real=1))

lowering/destruction/a† ladder operator. Given an eignenstate of the number operator `|n⟩` on a
`QSystem q`, `lowering(q)|n⟩ = √(n+1)|n+1⟩`. Optionally apply an additional phase ϕ and scaling `factor`.

# Examples
```jldoctest
julia> q = DuffingTransmon("q", 3, DuffingSpec(5, -0.2))
DuffingTransmon("q", 3, DuffingSpec(5, -0.2))

julia> lowering(q)
3×3 Array{Float64,2}:
 0.0  1.0  0.0
 0.0  0.0  1.41421
 0.0  0.0  0.0
```
"""
lowering(q::QSystem, ϕ::Real=0, factor::Real=1) = lowering(dimension(q), ϕ, factor)

# custom diagm because it's twice as fast for small matrices
function lowering(dim::Integer,  ϕ::Real=0, factor::Real=1)
    m = zeros(typeof(complex(float(ϕ))), dim, dim)
    efac = exp(-1im * 2π * ϕ) * factor
    for i in 1:(dim - 1)
        m[i, i + 1] = sqrt(i) * efac
    end
    return m
end

"""
    number(q::QSystem, factor::Real=1)

Number operator (a†a) on a QSystem `q` with an optional scaling factor which gives a dephasing
matrix if `factor != 1`.
"""
number(q::QSystem, factor::Real=1) = number(dimension(q), factor)

# custom diagm because it's twice as fast for small matrices
function number(dim::Integer, factor::Real=1)
    m = zeros(factor |> float |> complex |> typeof, dim, dim)
    for i in 2:dim
        m[i, i] = (i - 1) * factor
    end
    return m
end


"""
    X(q::QSystem)

X (a + a†) operator on a QSystem `q`.
"""
function X(q::QSystem)
    # Rather than just returning raising(q) + lowering(q) we specialize this implementation because
    # it more than twice as fast with less than half the allocations. See a4add5ca0 for details.
    diag_elements = sqrt.(1:(dimension(q)-1))
    diagm(-1 => diag_elements, +1 => diag_elements)
end

"""
Apply a phase ϕ (units of τ) to an X operator

# Examples
```jldoctest
julia> q = DuffingTransmon("q", 2, DuffingSpec(5, -0.2))
DuffingTransmon("q", 2, DuffingSpec(5, -0.2))

julia> X(q, 0.25)
2×2 Array{Complex{Float64},2}:
         0.0+0.0im  6.12323e-17-1.0im
 6.12323e-17+1.0im          0.0+0.0im

julia> X(q, 0.25) ≈ Y(q)
true
```
"""
X(q::QSystem, ϕ::Real) = raising(q, ϕ) + lowering(q, ϕ)

X(qs::Vector{<:QSystem}, ϕs::Vector{<:Real}) = reduce(⊗, [X(q, ϕ) for (q, ϕ) in zip(qs, ϕs)])
X(qs::Vector{<:QSystem}) = reduce(⊗, [X(q) for q in qs])

"""
    Y(q::QSystem)

Y (-ia + ia†) operator on a QSystem `q`.
"""
function Y(q::QSystem)
    diag_elements = sqrt.(1:(dimension(q)-1))
    diagm(-1 => 1im*diag_elements, +1 => -1im*diag_elements)
end

Y(q::QSystem, ϕ::Real) = 1im*(raising(q, ϕ) - lowering(q, ϕ))

Y(qs::Vector{<:QSystem}, ϕs::Vector{<:Real}) = reduce(⊗, [Y(q, ϕ) for (q, ϕ) in zip(qs, ϕs)])
Y(qs::Vector{<:QSystem}) = reduce(⊗, [Y(q) for q in qs])

X_Y(qs::Vector{<:QSystem}, ϕs::Vector{<:Real}) = X(qs, ϕs) + Y(qs, ϕs)
X_Y(qs::Vector{<:QSystem}) = X(qs) + Y(qs)


"""
    flip_flop(a::QSystem, b::QSystem)

Bilinear photon exchange Hamiltonian on `a` and `b`: `ab† + a†b`. For qubits corresponds to Pauli operator
`σ⁺σ⁻ + σ⁻σ⁺`
"""
flip_flop(a::QSystem, b::QSystem) = raising(a)⊗lowering(b) + lowering(a)⊗raising(b)

"""
    flip_flop(a::QSystem, b::QSystem, ϕ::Real)

Apply an additional relative phase ϕ: `exp(2πiϕ)ab† + exp(-2πiϕ)a†b`
"""
flip_flop(a::QSystem, b::QSystem, ϕ::Real) = exp(-1im*2π*ϕ)*raising(a)⊗lowering(b) + exp(1im*2π*ϕ)*lowering(b)⊗raising(b)

"""
    XY(a::QSystem, b::QSystem)

Bilinear XY Hamiltonian on `a` and `b` which is (up to a scale) equivalent to a "flip-flop"
Hamiltonian. XY(a,b) = XᵃXᵇ + YᵃYᵇ = 2*(ab† + a†b)
"""
XY(a::QSystem, b::QSystem) = lmul!(2.0, flip_flop(a,b))

"""
    XY(a::QSystem, b::QSystem, ϕ::Real)

Apply an additional phase rotation to the XY Hamiltonian 2*(exp(2πiϕ)ab† + exp(-2πiϕ)a†b)
"""
XY(a::QSystem, b::QSystem, ϕ::Real) =  lmul!(2.0, flip_flop(a,b,ϕ))


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
    return lowering(qs, 0, sqrt(γ))
end

"""
    dephasing(qs::QSystem, γ::Real)

Dephasing for a QSystem.

## args
* `qs`: a QSystem.
* `γ`: a decay rate in frequency units. Note Tϕ = 1/(2πγ).

## returns
The lindblad operator for dephasing.
"""
function dephasing(qs::QSystem, γ::Real)
    return number(qs, sqrt(2γ))
end

"""
    dipole_drive(qs::QSystem, drive::Function, rotation_rate::Real=0.0)

Given some function of time, return a function applying a time dependent
dipole Hamiltonian. Note that this does not use the rotating wave approximation
and therefore requires a real valued drive. See also `rwa_dipole`

## args
* `qs`: a QSystem.
* `drive`: a function of time returning a real value.
* `rotation_rate`: the rotation rate of a rotating frame.

## returns
A function of time.
"""
function dipole_drive(qs::QSystem, drive::Function, rotation_rate::Real=0.0)
    function ham(t)
        pulse::Real = drive(t)
        return pulse * X(qs, rotation_rate * t)
    end
    return ham
end

"""
    rwa_dipole(qs::QSystem, drive::Function)

Given some function of time, return a function applying a time dependent
dipole Hamiltonian under the rotating wave approximation (RWA).

## args
* `qs`: a QSystem.
* `drive`: a function of time returning a real or complex value. The real
    part couples to X and the imaginary part couples to Y.

## returns
A function of time.
"""
function rwa_dipole(qs::QSystem, drive::Function)
    x_ham = X(qs)
    y_ham = Y(qs)
    function ham(t)
        pulse = drive(t)
        return real(pulse) * x_ham + imag(pulse) * y_ham
    end
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
