using Optim: optimize
using LinearAlgebra: diag, diagm, eigvals

export QSpec, TransmonSpec, DuffingSpec, ResonatorSpec, HermitianSpec
export QSystem, label, dim, spec
export LiteralHermitian, Resonator, DuffingTransmon, PerturbativeTransmon, ChargeBasisTransmon
export hamiltonian, fit_transmon

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

# required functions with defaults
label(q::QSystem) = q.label
dim(q::QSystem) = q.dim
spec(q::QSystem) = q.spec

######################################################
# LiteralHermitian, an arbitrary Hamiltonian
######################################################

struct LiteralHermitian <: QSystem
    label::AbstractString
    spec::HermitianSpec
end

dim(h::LiteralHermitian) = size(h.spec.matrix, 1)

hamiltonian(q::LiteralHermitian) = spec(q).matrix

######################################################
# Resonator, a linear resonator
######################################################

struct Resonator <: QSystem
    label::AbstractString
    dim::Int
    spec::ResonatorSpec
end

hamiltonian(r::Resonator) = diagm(0 => [spec(r).frequency * n for n in 0:dim(r)-1])

######################################################
# DuffingTransmon
######################################################

function duffing_hamiltonian(freq::Real, anharm::Real, dimension::Int)
    n = collect(0:dimension-1)
    return diagm(0 => (freq - 0.5 * anharm) * n + 0.5 * anharm * n.^2)
end

struct DuffingTransmon <: QSystem
    label::AbstractString
    dim::Int
    spec::DuffingSpec
end

hamiltonian(t::DuffingTransmon) = duffing_hamiltonian(spec(t).frequency, spec(t).anharmonicity, dim(t))

hamiltonian(t::DuffingTransmon, frequency::Real) = duffing_hamiltonian(frequency, spec(t).anharmonicity, dim(t))

######################################################
# PerturbativeTransmon
######################################################

struct PerturbativeTransmon <: QSystem
    label::AbstractString
    dim::Int
    spec::TransmonSpec
    num_terms::Int

    PerturbativeTransmon(label::AbstractString, dim::Int, spec::TransmonSpec; num_terms::Int=PERTURBATIVE_NUM_TERMS) =
                new(label::AbstractString, dim::Int, spec::TransmonSpec, PERTURBATIVE_NUM_TERMS)
end

function DuffingSpec(t::TransmonSpec, ϕ::Real=0.0, num_terms::Int=PERTURBATIVE_NUM_TERMS)
    freq = perturbative_transmon_freq(t.EC, t.EJ1, t.EJ2, ϕ, num_terms=num_terms)
    anharm = perturbative_transmon_anharm(t.EC, t.EJ1, t.EJ2, ϕ, num_terms=num_terms)
    return DuffingSpec(freq, anharm)
end

function hamiltonian(t::PerturbativeTransmon, ϕ::Real=0.0)
    s = DuffingSpec(spec(t), ϕ, t.num_terms)
    return duffing_hamiltonian(s.frequency, s.anharmonicity, dim(t))
end

######################################################
# Charge basis transmon
######################################################

const CHARGE_NUM_TERMS = 101

struct ChargeBasisTransmon <: QSystem
    label::AbstractString
    dim::Int
    spec::TransmonSpec
    num_terms::Int
    function ChargeBasisTransmon(label::AbstractString, dim::Int, spec::TransmonSpec; num_terms::Int=CHARGE_NUM_TERMS)
        @assert mod(num_terms, 2) == 1
        @assert dim <= num_terms
        return new(label, dim, spec, num_terms)
    end
end


function hamiltonian(t::ChargeBasisTransmon, ϕ::Real=0.0)
    d = t.num_terms
    N = floor(Int, d/2)
    s = spec(t)
    EJ = sqrt(s.EJ1^2 + s.EJ2^2 + 2 * s.EJ1 * s.EJ2 * cos(2π * ϕ))
    EC = s.EC
    charging_term = 4 * EC * diagm(0 => (-N:N).^2)
    tunneling_term = -0.5 * EJ * (diagm(-1 => ones(d-1), 1 => ones(d-1)))
    # since the Hamiltonian is Hermitian the eigenvalues should already be sorted by the LAPACK
    # solver; however,  since that implementation is not guaranteed by Julia,  belts and suspenders
    # style we sort again
    return diagm(0 => sort(real(eigvals(charging_term + tunneling_term)))[1:dim(t)])
end
# TODO: update raising and lowering to be correct for ChargeBasisTransmon


######################################################
# Fit transmons
######################################################
"""
    fit_transmon(freq_max::Real, freq_min::Real, anharm_max::Real, model::Type{T}, num_terms::Int) where {T<:QSystem}

Fit a transmon to given frequency and anharmonicity at fmax and the frequency
at fmin. The model can be either PerturbativeTransmon or ChargeBasisTransmon.
If `freq_max == freq_min` then `EJ2 == 0` is enforced.

## args
* `freq_max`: the qubit maximum frequency.
* `freq_min`: the qubit minimum frequency.
* `anharm_max`: the qubit maximum anharmonicity.
* `model`: either PerturbativeTransmon or ChargeBasisTransmon.
* `num_terms`: the number of terms used in the model.

## returns
A TransmonSpec.
"""
function fit_transmon(freq_max::Real, freq_min::Real, anharm_max::Real, model::Type{T}, num_terms::Int) where {T<:QSystem}
    function f_fixed(params)
      EC, EJ = params
      t = model("", 3, TransmonSpec(EC, EJ, 0.0), num_terms=num_terms)
      levels = diag(hamiltonian(t))
      test_f_01 = levels[2] - levels[1]
      test_f_12 = levels[3] - levels[2]
      return abs(test_f_01 - freq_max) + abs(test_f_12 - test_f_01 - anharm_max)
    end

    function f_tunable(params)
        t = model("", 3, TransmonSpec(params...), num_terms=num_terms)
        levels_fmax = diag(hamiltonian(t, 0.0))
        levels_fmin = diag(hamiltonian(t, 0.5))
        test_f_01_max = levels_fmax[2] - levels_fmax[1]
        test_f_12_max = levels_fmax[3] - levels_fmax[2]
        test_f_01_min = levels_fmin[2] - levels_fmin[1]
        return abs(test_f_01_max - freq_max) + abs(test_f_12_max - test_f_01_max - anharm_max) + abs(test_f_01_min - freq_min)
    end

    EC_guess = -anharm_max
    EJ_max_guess = -(freq_max - anharm_max)^2 / (8 * anharm_max)
    EJ_min_guess = -(freq_min - anharm_max)^2 / (8 * anharm_max)
    EJ1_guess = .5 * (EJ_max_guess + EJ_min_guess)
    EJ2_guess = .5 * abs(EJ_max_guess - EJ_min_guess)
    if freq_max == freq_min
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
    @warn "Deprecation warning: Resonator."
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
    @warn "Deprecation warning: TunableTransmon."
    EJ1, EJ2 = asymmetry_to_EJs(E_J, d)
    return ChargeBasisTransmon(label, dim, TransmonSpec(E_C, EJ1, EJ2), num_terms=dim)
end

function FixedTransmon(label::AbstractString,
    E_C::Real,
    E_J::Real, # sum of junction E_J's
    dim::Int)
    @warn "Deprecation warning: FixedTransmon."
    EJ1, EJ2 = E_J, 0.0
    return ChargeBasisTransmon(label, dim, TransmonSpec(E_C, EJ1, EJ2), num_terms=dim)
end

function FixedDuffingTransmon(label::AbstractString,
    frequency::Real,
    anharmonicity::Real,
    dim::Int)
    @warn "Deprecation warning: FixedDuffingTransmon."
    return DuffingTransmon(label, dim, DuffingSpec(frequency, anharmonicity))
end

function TunableDuffingTransmon(label::AbstractString,
    E_C::Real,
    E_J::Real, #sum of junction E_J's
    d::Real, #asymmetry parameter
    dim::Int)
    @warn "Deprecation warning: TunableDuffingTransmon."
    EJ1, EJ2 = asymmetry_to_EJs(E_J, d)
    return PerturbativeTransmon(label, dim, TransmonSpec(E_C, EJ1, EJ2), num_terms=2)
end

function MathieuTransmon(label::AbstractString,
    E_C::Real,
    E_J::Real, #sum of junction E_J's
    d::Real, #asymmetry parameter
    dim::Int)
    @warn "Deprecation warning: MathieuTransmon."
    EJ1, EJ2 = asymmetry_to_EJs(E_J, d)
    return PerturbativeTransmon(label, dim, TransmonSpec(E_C, EJ1, EJ2))
end

function create(q::QSystem)
    @warn "Deprecation warning: create."
    return raising(q)
end

function destroy(q::QSystem)
    @warn "Deprecation warning: destroy."
    return lowering(q)
end

function fit_fixed_transmon(freq, anharm, num_terms)
    @warn "Deprecation waring: fit_fixed_transmon."
    t = fit_transmon(freq, freq, anharm, ChargeBasisTransmon, num_terms)
    return t.EC, t.EJ1 + t.EJ2
end

function fit_tunable_transmon(freq_max, freq_min, anharm_max, num_terms, model::Type{T}) where {T<:QSystem}
    @warn "Deprecation waring: fit_tunable_transmon."
    t = fit_transmon(freq_max, freq_min, anharm_max, model, num_terms)
    EJ, d = EJs_to_asymmetry(t.EJ1, t.EJ2)
    return t.EC, EJ, d
end
