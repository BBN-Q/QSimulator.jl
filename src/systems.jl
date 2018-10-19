using Optim: optimize

export QSpec, TransmonSpec, DuffingSpec, ResonatorSpec, HermitianSpec
export QSystem, label, dim, spec
export LiteralHermitian, Resonator, DuffingTransmon, PerturbativeTransmon
export ChargeBasisTransmon, CHARGE_NUM_TERMS
export hamiltonian, duffing_from_transmon, fit_transmon
export scale_EJ

######################################################
# QSpec - the physical parameters defining a quantum system
######################################################
abstract type QSpec end

struct TransmonSpec <: QSpec
   EC::Real
   EJ1::Real
   EJ2::Real
end

struct DuffingSpec <: QSpec
   frequency::Real
   anharmonicity::Real
end

struct ResonatorSpec <: QSpec
   frequency::Real
end

struct HermitianSpec <: QSpec
   matrix::Matrix{<:Number}
   function HermitianSpec(matrix::Matrix{<:Number})
       @assert ishermitian(matrix)
       return new(matrix)
   end
end

######################################################
# QSystem - a quantum system
######################################################
abstract type QSystem end

# required functions
label(q::QSystem) = q.label
dim(q::QSystem) = q.dim
spec(q::QSystem) = q.spec

######################################################
# LiteralHermitian, an arbitrary Hamiltonian
######################################################

struct LiteralHermitian <: QSystem
    label::AbstractString
    dim::Int
    spec::HermitianSpec

    LiteralHermitian(label::AbstractString, hermitian_spec::HermitianSpec) = new(label, size(hermitian_spec.matrix, 1), hermitian_spec)
end

hamiltonian(q::LiteralHermitian) = spec(q).matrix

######################################################
# Resonator, a linear resonator
######################################################

struct Resonator <: QSystem
    label::AbstractString
    dim::Int
    spec::ResonatorSpec
end

hamiltonian(r::Resonator) = diagm([spec(r).frequency * n for n in 0:dim(r)-1])

######################################################
# DuffingTransmon
######################################################

function duffing_hamiltonian(ω::Real, η::Real, dimension::Int)
    n = collect(0:dimension-1)
    return diagm((ω - 0.5 * η) * n + 0.5 * η * n.^2)
end

struct DuffingTransmon <: QSystem
    label::AbstractString
    dim::Int
    spec::DuffingSpec
end

hamiltonian(t::DuffingTransmon) = duffing_hamiltonian(spec(t).frequency, spec(t).anharmonicity, dim(t))

######################################################
# PerturbativeTransmon
######################################################

struct PerturbativeTransmon <: QSystem
    label::AbstractString
    dim::Int
    spec::TransmonSpec
end

function hamiltonian(t::PerturbativeTransmon, ϕ::Real=0.0)
    s = spec(t)
    ω = perturbative_transmon_freq(s.EC, s.EJ1, s.EJ2, ϕ, num_terms=PERTURBATIVE_NUM_TERMS)
    η = perturbative_transmon_anharm(s.EC, s.EJ1, s.EJ2, ϕ, num_terms=PERTURBATIVE_NUM_TERMS)
    return duffing_hamiltonian(ω, η, dim(t))
end

function duffing_from_transmon(s::TransmonSpec, ϕ::Real=0.0)
    ω = perturbative_transmon_freq(s.EC, s.EJ1, s.EJ2, ϕ, num_terms=PERTURBATIVE_NUM_TERMS)
    η = perturbative_transmon_anharm(s.EC, s.EJ1, s.EJ2, ϕ, num_terms=PERTURBATIVE_NUM_TERMS)
    return DuffingSpec(ω, η)
end

######################################################
# Charge basis transmon
######################################################

struct ChargeBasisTransmon <: QSystem
    label::AbstractString
    dim::Int
    spec::TransmonSpec
    function ChargeBasisTransmon(label::AbstractString, dim::Int, spec::TransmonSpec)
        @assert mod(dim, 2) == 1
        return new(label, dim, spec)
    end
end

function hamiltonian(t::ChargeBasisTransmon, ϕ::Real=0.0)
    N = floor(Int, dim(t)/2)
    s = spec(t)
    EJ = sqrt(s.EJ1^2 + s.EJ2^2 + 2 * s.EJ1 * s.EJ2 * cos(2π * ϕ))
    EC = s.EC
    charging_term = 4 * EC * diagm((-N:N).^2)
    tunneling_term = -0.5 * EJ * (diagm(ones(dim(t)-1),-1) + diagm(ones(dim(t)-1),1))
    return charging_term + tunneling_term
end

# TODO: update raising and lowering to be correct for ChargeBasisTransmon


######################################################
# Fit transmons
######################################################
"""
    fit_transmon(f_max::Real, f_min::Real, η_max::Real, model::Type{T}, num_terms::Int) where {T<:QSystem}

Fit a transmon to given frequency and anharmonicity at fmax and the frequency
at fmin. The model can be either PerturbativeTransmon or ChargeBasisTransmon.
If `f_max == ω_min` then `EJ2 == 0` is enforced.

## args
* `f_max`: the qubit maximum frequency.
* `f_min`: the qubit minimum frequency.
* `η_max`: the qubit maximum anharmonicity.
* `model`: either PerturbativeTransmon or ChargeBasisTransmon.
*   `dim`: the dimension of the model.

## returns
A TransmonSpec.
"""
function fit_transmon(f_max::Real, f_min::Real, η_max::Real, model::Type{T}, dim::Int) where {T<:QSystem}
    function f_fixed(params)
      EC, EJ = params
      t = model("", dim, TransmonSpec(EC, EJ, 0.0))
      levels = real(sort(eigvals(hamiltonian(t))))
      test_f_01 = levels[2] - levels[1]
      test_f_12 = levels[3] - levels[2]
      return abs(test_f_01 - f_max) + abs(test_f_12 - test_f_01 - η_max)
    end

    function f_tunable(params)
        t = model("", dim, TransmonSpec(params...))
        levels_fmax = real(sort(eigvals(hamiltonian(t, 0.0))))
        levels_fmin = real(sort(eigvals(hamiltonian(t, 0.5))))
        test_f_01_max = levels_fmax[2] - levels_fmax[1]
        test_f_12_max = levels_fmax[3] - levels_fmax[2]
        test_f_01_min = levels_fmin[2] - levels_fmin[1]
        return abs(test_f_01_max - f_max) + abs(test_f_12_max - test_f_01_max - η_max) + abs(test_f_01_min - f_min)
    end

    EC_guess = abs(η_max)
    EJ_max_guess = (f_max + EC_guess)^2 / (8 * EC_guess)
    EJ_min_guess = (f_min + EC_guess)^2 / (8 * EC_guess)
    EJ1_guess = .5 * (EJ_max_guess + EJ_min_guess)
    EJ2_guess = .5 * abs(EJ_max_guess - EJ_min_guess)
    if f_max == f_min
        res = optimize(f_fixed, [EC_guess, EJ1_guess])
        EC, EJ = res.minimizer
        return TransmonSpec(EC, EJ, 0.0)
    else
        res = optimize(f_tunable, [EC_guess, EJ1_guess, EJ2_guess])
        return TransmonSpec(res.minimizer...)
    end
end

######################################################
# Backwards compatibility
######################################################

export TunableTransmon, FixedTransmon, FixedDuffingTransmon, TunableDuffingTransmon
export fit_fixed_transmon, fit_tunable_transmon
export MathieuTransmon, create, destroy

function Resonator(label::AbstractString, frequency::Real, dim::Int)
    warn("Deprecation warning: Resonator.")
    return Resonator(label, dim, ResonatorSpec(frequency))
end

function asymmetry_to_EJs(EJ::Real, d::Real)
    EJ1 = .5 * (d + 1) * EJ
    EJ2 = EJ - EJ1
    return EJ1, EJ2
end

function EJs_to_asymmetry(EJ1::Real, EJ2::Real)
    EJΣ = EJ1 + EJ2
    EJΔ = EJ1 - EJ2
    return EJΣ, abs(EJΔ/EJΣ)
end

function TunableTransmon(label::AbstractString,
    E_C::Real,
    E_J::Real, # sum of junction E_J's
    d::Real, # juntion asymmetry parameter
    dim::Int)
    warn("Deprecation warning: TunableTransmon.")
    EJ1, EJ2 = asymmetry_to_EJs(E_J, d)
    return ChargeBasisTransmon(label, dim, TransmonSpec(E_C, EJ1, EJ2))
end

function FixedTransmon(label::AbstractString,
    E_C::Real,
    E_J::Real, # sum of junction E_J's
    dim::Int)
    warn("Deprecation warning: FixedTransmon.")
    EJ1, EJ2 = E_J, 0.0
    return ChargeBasisTransmon(label, dim, TransmonSpec(E_C, EJ1, EJ2))
end

function FixedDuffingTransmon(label::AbstractString,
    frequency::Real,
    anharmonicity::Real,
    dim::Int)
    warn("Deprecation warning: FixedDuffingTransmon.")
    return DuffingTransmon(label, dim, DuffingSpec(frequency, anharmonicity))
end

function TunableDuffingTransmon(label::AbstractString,
    E_C::Real,
    E_J::Real, #sum of junction E_J's
    d::Real, #asymmetry parameter
    dim::Int)
    warn("Deprecation warning: TunableDuffingTransmon.")
    EJ1, EJ2 = asymmetry_to_EJs(E_J, d)
    return PerturbativeTransmon(label, dim, TransmonSpec(E_C, EJ1, EJ2))
end

function MathieuTransmon(label::AbstractString,
    E_C::Real,
    E_J::Real, #sum of junction E_J's
    d::Real, #asymmetry parameter
    dim::Int)
    warn("Deprecation warning: MathieuTransmon.")
    EJ1, EJ2 = asymmetry_to_EJs(E_J, d)
    return PerturbativeTransmon(label, dim, TransmonSpec(E_C, EJ1, EJ2))
end

function fit_fixed_transmon(ω, η, dim)
    warn("Deprecation waring: fit_fixed_transmon.")
    t = fit_transmon(ω, ω, η, ChargeBasisTransmon, dim)
    return t.EC, t.EJ1 + t.EJ2
end

function fit_tunable_transmon(f_max, f_min, η_max, dim, model::Type{T}) where {T<:QSystem}
    warn("Deprecation waring: fit_tunable_transmon.")
    t = fit_transmon(f_max, f_min, η_max, model, dim)
    EJ, d = EJs_to_asymmetry(t.EJ1, t.EJ2)
    return t.EC, EJ, d
end

function scale_EJ(E_J, flux, d)
    flux_rad = π*flux
    E_J * sqrt(cos(flux_rad)^2 + d^2*(sin(flux_rad)^2))
end
