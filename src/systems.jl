export QSystem,
       Resonator,
       TunableTransmon,
       FixedTransmon,
       FixedDuffingTransmon,
       TunableDuffingTransmon,
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

export hamiltonian

# atomic  quantum systems such as resonators or transmons
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

# linear resonator
mutable struct Resonator <: QSystem
    label::AbstractString
    frequency::Float64
    dim::Int
end
hamiltonian(r::Resonator) = r.frequency * number(r)

""" Calculate the drift or natural Hamiltonian of a CompositeQSystem """
function hamiltonian(cqs::CompositeQSystem)
    ham = zeros(Complex128, (dim(cqs), dim(cqs)))
    for (new_ham, idxs) = cqs.fixed_Hs
        embed_add!(ham, new_ham, idxs)
    end
    return ham
end

# full Transmon
mutable struct TunableTransmon <: QSystem
    label::AbstractString
    E_C::Float64
    E_J::Float64 # sum of junction E_J's
    d::Float64 # juntion asymmetry parameter
    dim::Int
end

""" Scale the effective E_J given an asymmetry  paramter and flux threading the loop  in Φ_0 """
function scale_EJ(E_J, flux, d)
    flux_rad = π*flux
    E_J * sqrt(cos(flux_rad)^2 + d^2*(sin(flux_rad)^2))
end

""" Transmon Hamiltonian in the charge basis """
function hamiltonian(t::TunableTransmon, flux)
  N = floor(Int, dim(t)/2)
  scaled_EJ = scale_EJ(t.E_J, flux, t.d)
  4 * t.E_C * diagm((-N:N).^2) - scaled_EJ  * 0.5 * (diagm(ones(dim(t)-1),-1) + diagm(ones(dim(t)-1),1))
end

mutable struct FixedTransmon <: QSystem
    label::AbstractString
    E_C::Float64
    E_J::Float64 # sum of junction E_J's
    dim::Int
end

""" Transmon Hamiltonian in the charge basis """
function hamiltonian(t::FixedTransmon)
  N = floor(Int, dim(t)/2)
  4 * t.E_C * diagm((-N:N).^2) - 0.5*t.E_J*(diagm(ones(dim(t)-1),-1) + diagm(ones(dim(t)-1),1))
end

mutable struct FixedDuffingTransmon <: QSystem
    label::AbstractString
    frequency::Float64
    anharmonicity::Float64
    dim::Int
end
hamiltonian(t::FixedDuffingTransmon) = (t.frequency - 0.5*t.anharmonicity)*number(t) + 0.5*t.anharmonicity * number(t)^2

mutable struct TunableDuffingTransmon <: QSystem
    label::AbstractString
    E_C::Float64
    E_J::Float64 #sum of junction E_J's
    d::Float64 #asymmetry parameter
    dim::Int
end

function hamiltonian(t::TunableDuffingTransmon, flux)
    scaled_EJ = scale_EJ(t.E_J, flux, t.d)
    ωₚ = sqrt(8*t.E_C*scaled_EJ)
    return diagm([(ωₚ - t.E_C / 2) * ct - t.E_C / 2 * ct^2 for ct in 0:(t.dim-1)])
end

"""
Fit  E_C and E_J for a fixed frequency transmon given f_01 and anharmonicity
"""
function fit_fixed_transmon(f_01, α, dim)
    # helper function for least squares difference between measured and calculated
    function f(params)
      E_C = params[1]
      E_J  = params[2]
      t = FixedTransmon("dummy", E_C, E_J, dim)
      levels = eigvals(hamiltonian(t))
      test_f_01 = levels[2] - levels[1]
      test_f_12 = levels[3] - levels[2]
      ( test_f_01 - f_01)^2 + ( test_f_12 -  test_f_01 - α)^2
    end
    res = optimize(f, [abs(α), (f_01+abs(α))^2 / (8 * abs(α))] )
    return res.minimizer[1], res.minimizer[2]
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
