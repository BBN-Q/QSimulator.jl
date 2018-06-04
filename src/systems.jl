using Optim

export QSystem,
       Resonator,
       TunableTransmon,
       FixedTransmon,
       FixedDuffingTransmon,
       TunableDuffingTransmon,
       MathieuTransmon,
       label,
       dim

export raising,
       create,
       lowering,
       destroy,
       number,
       X,
       Y,
       flip_flop,
       XY,
       dipole

""" atomic  quantum systems such as resonators or transmons """
abstract type QSystem end

label(q::QSystem) = q.label
dim(q::QSystem) = q.dim

# some basic operators
raising(q::QSystem) = diagm(sqrt.(1:(dim(q)-1)), -1)
const create = raising
lowering(q::QSystem) = diagm(sqrt.(1:(dim(q)-1)), 1)
const destroy = lowering
number(q::QSystem) = raising(q) * lowering(q)
X(q::QSystem) = raising(q) + lowering(q)
Y(q::QSystem) = 1im*(raising(q) - lowering(q))

""" linear resonator """
mutable struct Resonator <: QSystem
    label::AbstractString
    frequency::Float64
    dim::Int
end

""" full Transmon """
mutable struct TunableTransmon <: QSystem
    label::AbstractString
    E_C::Float64
    E_J::Float64 # sum of junction E_J's
    d::Float64 # juntion asymmetry parameter
    dim::Int
end

""" Scale the effective E_J given an asymmetry paramter and flux threading the loop  in Φ_0 """
function scale_EJ(E_J, flux, d)
    flux_rad = π*flux
    E_J * sqrt(cos(flux_rad)^2 + d^2*(sin(flux_rad)^2))
end

mutable struct FixedTransmon <: QSystem
    label::AbstractString
    E_C::Float64
    E_J::Float64 # sum of junction E_J's
    dim::Int
end

mutable struct FixedDuffingTransmon <: QSystem
    label::AbstractString
    frequency::Float64
    anharmonicity::Float64
    dim::Int
end

mutable struct TunableDuffingTransmon <: QSystem
    label::AbstractString
    E_C::Float64
    E_J::Float64 #sum of junction E_J's
    d::Float64 #asymmetry parameter
    dim::Int
end

mutable struct MathieuTransmon <: QSystem
    label::AbstractString
    E_C::Float64
    E_J::Float64 #sum of junction E_J's
    d::Float64 #asymmetry parameter
    dim::Int
end

"""
Fit  E_C and E_J for a fixed frequency transmon given f_01 and anharmonicity
"""
function fit_fixed_transmon(f_01, α)
    # helper function for least squares difference between measured and calculated
    function f(params)
      E_C = params[1]
      E_J  = params[2]
      t = FixedTransmon("dummy", E_C, E_J, 3)
      levels = sort(real(eigvals(hamiltonian(t))))
      test_f_01 = levels[2] - levels[1]
      test_f_12 = levels[3] - levels[2]
      return ( test_f_01 - f_01)^2 + ( test_f_12 -  test_f_01 - α)^2
    end
    res = optimize(f, [abs(α), (f_01+abs(α))^2 / (8 * abs(α))] )
    return res.minimizer[1], res.minimizer[2]
end

"""
Fit  E_C, E_J, and d for either the TunableTransmon or TunableDuffingTransmon
model given the qubit frequencies at 0 and 1/2 flux, and the anharmonicity at 0 flux.
"""
function fit_tunable_transmon{T<:QSystem}(f_01_max, f_01_min, α_max, model::Type{T})
    # helper function for least squares difference between measured and calculated
    function f(params)
      E_C = params[1]
      E_J  = params[2]
      d = params[3]
      t = model("dummy", E_C, E_J, d, 3)
      levels_fmax = sort(real(eigvals(hamiltonian(t, 0))))
      levels_fmin = sort(real(eigvals(hamiltonian(t, 0.5))))
      test_f_01_max = levels_fmax[2] - levels_fmax[1]
      test_f_12_max = levels_fmax[3] - levels_fmax[2]
      test_f_01_min = levels_fmin[2] - levels_fmin[1]
      (test_f_01_max - f_01_max)^2 + (test_f_12_max -  test_f_01_max - α_max)^2 +
      (test_f_01_min - f_01_min)^2
    end

    EC_guess = abs(α_max)
    EJ_guess = (f_01_max+abs(α_max))^2 / (8 * abs(α_max))
    d_guess = (f_01_min+abs(α_max))^2 / (8 * EJ_guess * abs(α_max))

    res = optimize(f, [EC_guess, EJ_guess, d_guess])
    return res.minimizer[1], res.minimizer[2], res.minimizer[3]
end

# two body Hamiltonian terms

"""
    flip_flop(a::QSystem, b::QSystem; ϕ=0.0)

FlipFlop or XY interaction between two QSystems. Optional phase ϕ given in units of τ.
"""
function flip_flop(a::QSystem, b::QSystem; ϕ=0.0)
    # TODO: confirm sign of phase term
    phase = exp(1im*2π*ϕ)
    return phase * raising(a) ⊗ lowering(b) + conj(phase) * lowering(a) ⊗ raising(b)
end

const XY = flip_flop

function dipole(a::QSystem, b::QSystem)
    return X(a) ⊗ X(b)
end
