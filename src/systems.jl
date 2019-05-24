using Optim: optimize
using LinearAlgebra: diag, diagm, eigvals, ishermitian

export QSpec, TransmonSpec, DuffingSpec, ResonatorSpec, HermitianSpec
export QSystem, label, dimension, spec
export LiteralHermitian, Resonator, DuffingTransmon
export PerturbativeTransmon, ChargeBasisTransmon, DiagonalChargeBasisTransmon, RotatingFrameSystem
export hamiltonian, fit_transmon

######################################################
# Types
######################################################

HermitianMatrix = Matrix{<:Number}
abstract type QSpec end
abstract type QSystem end

# Required QSystem functions with defaults
label(q::QSystem)::AbstractString = q.label
dimension(q::QSystem)::Int64 = q.dimension
spec(q::QSystem)::QSpec = q.spec

################################################################################
# QSpecs
# QSpecs are the physical parameters that define a quantum system
################################################################################

# Structs

struct HermitianSpec <: QSpec
   matrix::HermitianMatrix
   function HermitianSpec(matrix::Matrix{<:Number})
       @assert ishermitian(matrix)
       return new(matrix)
   end
end

"""
Resonator, specified by its transition frequency.
"""
struct ResonatorSpec <: QSpec
   frequency::Real
end

"""
Transmon, specified by its charge and junctions energies, EC and EJ1/EJ2,
respectively.
"""
struct TransmonSpec <: QSpec
   EC::Real
   EJ1::Real
   EJ2::Real
end

"""
Alternative Transmon specifiication by its (base linear) transition frequency
and anharmonicity.
"""
struct DuffingSpec <: QSpec
   frequency::Real
   anharmonicity::Real
end

# Functions

"""
Construct a TransmonSpec (EC, EJ1, EJ2) from a DuffingSpec (frequency,
anharmonicity), and a model type to perform the conversion.

An optimization routine, fit_transmon, must be run to convert
frequency/anharmonicity to EC/EJ.
"""
function TransmonSpec(t::DuffingSpec, model::Type{T}=PerturbativeTransmon,
    num_terms::Int=PERTURBATIVE_NUM_TERMS) where {T<:QSystem}
    return fit_transmon(t.frequency, t.frequency, t.anharmonicity, model, num_terms)
end

"""
Construct a DuffingSpec (freqency, anharmonicity) from a TransmonSpec (EC, EJ1,
EJ2), at a given flux quantum ϕ. May optionally supply a customer number of
perturbative terms.
"""
function DuffingSpec(t::TransmonSpec, ϕ::Real=0.0, num_terms::Int=PERTURBATIVE_NUM_TERMS)
    freq = perturbative_transmon_freq(t.EC, t.EJ1, t.EJ2, ϕ, num_terms=num_terms)
    anharm = perturbative_transmon_anharm(t.EC, t.EJ1, t.EJ2, ϕ, num_terms=num_terms)
    return DuffingSpec(freq, anharm)
end

################################################################################
# QSystems
# QSystems are Structs comprised of a QSpec, a label and a dimension.
################################################################################

# Structs

struct LiteralHermitian <: QSystem
    label::AbstractString
    spec::HermitianSpec
end

dimension(h::LiteralHermitian)::Int64 = size(h.spec.matrix, 1)
hamiltonian(q::LiteralHermitian)::HermitianMatrix = spec(q).matrix

struct Resonator <: QSystem
    label::AbstractString
    dimension::Int
    spec::ResonatorSpec
end

struct DuffingTransmon <: QSystem
    label::AbstractString
    dimension::Int
    spec::DuffingSpec
end

"""
The PerturbativeTransmon type is a DuffingTransmon with the added `num_terms`
attribute, specifying to how many terms the transmon should be calculated.
"""
struct PerturbativeTransmon <: QSystem
    label::AbstractString
    dimension::Int
    spec::TransmonSpec
    num_terms::Int

    PerturbativeTransmon(
        label::AbstractString, dimension::Int, spec::TransmonSpec;
        num_terms::Int=PERTURBATIVE_NUM_TERMS
    ) =
        new(label::AbstractString, dimension::Int, spec::TransmonSpec, num_terms)
end

struct ChargeBasisTransmon <: QSystem
    label::AbstractString
    dimension::Int
    spec::TransmonSpec
    function ChargeBasisTransmon(label::AbstractString, dimension::Int, spec::TransmonSpec)
        @assert mod(dimension, 2) == 1
        return new(label, dimension, spec)
    end
end

const CHARGE_NUM_TERMS = 101

struct DiagonalChargeBasisTransmon <: QSystem
    label::AbstractString
    dimension::Int
    spec::TransmonSpec
    num_terms::Int
    function DiagonalChargeBasisTransmon(label::AbstractString, dimension::Int, spec::TransmonSpec; num_terms::Int=CHARGE_NUM_TERMS)
        @assert mod(num_terms, 2) == 1
        @assert dimension <= num_terms
        return new(label, dimension, spec, num_terms)
    end
end

struct RotatingFrameSystem <: QSystem
    label::AbstractString
    system::QSystem
    rotation_rate::Real
end

# contruct from an existing QSystem without a different label
RotatingFrameSystem(q::QSystem, rotation_rate::Real) = RotatingFrameSystem(label(q), q, rotation_rate)

# defer `dimension` and `spec` to underlying QSystem
dimension(r::RotatingFrameSystem) = dimension(r.system)
spec(r::RotatingFrameSystem) = spec(r.system)

################################################################################
# Hamiltonians
################################################################################

# Resonator

"""
Construct a resonator hamiltonian from energy levels equally spaced
by the resonator frequency, along the diagonal of a square matrix.
"""
hamiltonian(r::Resonator)::HermitianMatrix =
    diagm(0 =>[spec(r).frequency * n for n in 0:dimension(r)-1])

# Duffing Transmon

"""
Construct the hamiltonian of a duffing transmon from energy levels spaced by the
transmon's frequency minus its anharmonicity:

```
f_{n+1} - f_n = frequency - anharmonicity.
```
"""
function duffing_hamiltonian(
    freq::Real, anharm::Real, dimension::Int
)::HermitianMatrix
    diagm(0 => [(freq - 0.5 * anharm) * n + 0.5 * anharm * n^2
                for n in 0:dimension-1])
end

"""
Construct a hamiltonian for a Duffing transmon.
"""
hamiltonian(t::DuffingTransmon)::HermitianMatrix =
    duffing_hamiltonian(spec(t).frequency, spec(t).anharmonicity, dimension(t))

"""
Construct a hamiltonian for a Duffing transmon using a custom frequency.
"""
hamiltonian(t::DuffingTransmon, frequency::Real)::HermitianMatrix =
    duffing_hamiltonian(frequency, spec(t).anharmonicity, dimension(t))

# Perturbative Transmon

"""
Construct a hamiltonian from a PerturbativeTransmon at a given flux quantum
ϕ.
"""
function hamiltonian(t::PerturbativeTransmon, ϕ::Real=0.0)::HermitianMatrix
    s = DuffingSpec(spec(t), ϕ, t.num_terms)
    return duffing_hamiltonian(s.frequency, s.anharmonicity, dimension(t))
end

# Charge Basis Transmon

function hamiltonian(t::ChargeBasisTransmon, ϕ::Real=0.0)
    d = dimension(t)
    N = floor(Int, d/2)
    s = spec(t)
    EJ = sqrt(s.EJ1^2 + s.EJ2^2 + 2 * s.EJ1 * s.EJ2 * cos(2π * ϕ))
    EC = s.EC
    charging_term = 4 * EC * diagm(0 => (-N:N).^2)
    tunneling_term = -0.5 * EJ * (diagm(-1 => ones(d-1), 1 => ones(d-1)))
    return charging_term + tunneling_term
end

function hamiltonian(t::DiagonalChargeBasisTransmon, ϕ::Real=0.0)
    h = hamiltonian(ChargeBasisTransmon(label(t), t.num_terms, spec(t)), ϕ)
    # since the Hamiltonian is Hermitian the eigenvalues should already be sorted by the LAPACK
    # solver; however,  since that implementation is not guaranteed by Julia, belts and suspenders
    # style we sort again
    return diagm(0 => sort(real(eigvals(h)))[1:dimension(t)])
end
# TODO: update raising and lowering to be correct for DiagonalChargeBasisTransmon

# Rotating Frame

function hamiltonian(r::RotatingFrameSystem)
    h = hamiltonian(r.system)
    subtract_number!(h, r.rotation_rate)
    return h
end

function hamiltonian(r::RotatingFrameSystem, t::Real)
    h = hamiltonian(r.system, t)
    subtract_number!(h, r.rotation_rate)
    return h
end

# Util

"""
Helper function that subtracts a real number off the diagonal elements of a Real
matrix.
"""
function subtract_number!(h::Matrix{<:Real}, r::Real)
    if r != 0
        for i in 1:size(h, 1)
            h[i,i] -= (i-1) * r
        end
    end
end

################################################################################
# Conversion fitting (EC/EJ <> Frequency/Anharmonicity)
################################################################################

"""
    fit_transmon(freq_max::Real, freq_min::Real, anharm_max::Real, model::Type{T}, num_terms::Int) where {T<:QSystem}

Fit a transmon to given frequency and anharmonicity at fmax and the frequency
at fmin. The model can be either PerturbativeTransmon or DiagonalChargeBasisTransmon.
If `freq_max == freq_min` then `EJ2 == 0` is enforced.

## args
* `freq_max`: the qubit maximum frequency.
* `freq_min`: the qubit minimum frequency.
* `anharm_max`: the qubit maximum anharmonicity.
* `model`: either PerturbativeTransmon or DiagonalChargeBasisTransmon.
* `num_terms`: the number of terms used in the model.

## returns
A TransmonSpec.
"""
function fit_transmon(
    freq_max::Real, freq_min::Real, anharm_max::Real, model::Type{T}, num_terms::Int
)::TransmonSpec where {T<:Union{PerturbativeTransmon,DiagonalChargeBasisTransmon}}
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
