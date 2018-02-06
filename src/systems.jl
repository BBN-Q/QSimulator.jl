export Resonator,
       TunableTransmon,
       FixedTransmon,
       FixedDuffingTransmon,
       TunableDuffingTransmon,
       hamiltonian,
       CompositeQSystem

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

# full Transmon
mutable struct TunableTransmon <: QSystem
    label::AbstractString
    E_C::Float64
    E_J::Float64 # sum of junction E_J's
    d::Float64 # juntion asymmetry parameter
    dim::Int
end

""" Scale the effective E_J given an asymmetry  paramter and flux threading the loop  in Φ_0 """
scale_EJ(flux, d) = sqrt(cos(π*flux)^2 + d^2*(sin(π*flux)^2))


""" Transmon Hamiltonian in the charge basis """
function hamiltonian(t::TunableTransmon, flux::Float64)
  N = floor(Int, dim(t)/2)
  scaled_EJ = t.E_J * scale_EJ(flux, t.d)
  4 * t.E_C * diagm((-N:N).^2) - scaled_EJ  * 0.5 * (diagm(ones(dim(t)-1),-1) + diagm(ones(dim(t)-1),1))
end


mutable struct FixedTransmon <: QSystem
    label::AbstractString
    E_C::Float64
    E_J::Float64 # sum of junction E_J's
    dim::Int
end

""" Transmon Hamiltonian in the charge basis """
function hamiltonian(t::FixedTransmon, flux::Float64)
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

function hamiltonian(t::TunableDuffingTransmon, flux::Float64)
    scaled_EJ = t.E_J * scale_EJ(flux, t.d)
    ωₚ = sqrt(8*t.E_C*scaled_EJ)
    omega = [(omega_p-tt.E_C/2)*ct - tt.E_C/2*ct^2 for ct in 0:(t.dim-1)]
    return diagm(omega)
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


# tensor products of quantum systems
mutable struct CompositeQSystem
    subsystems::Vector{Any}
    fixed_Hs::Vector{Any} # tuple of Matrix and exansion indices
    parametric_Hs::Vector{Any} # tuple of Functions and expansion indices
    dim::Int
end

CompositeQSystem(qs) = CompositeQSystem(qs, [], [], prod(dim(q) for q in qs))

""" Calculate the drift or natural Hamiltonian of a CompositeQSystem """
function hamiltonian(cqs::CompositeQSystem)
    ham = zeros(Complex128, (dim(cqs), dim(cqs)))
    for h = cqs.fixed_Hs
        expand_add!(ham, new_ham, expand_idxs)
    end
end


""" In place addition of an operator embedded into a larger Hilbert space given a set of expansion indices"""
function expand_add!(op, new_op, expand_idxs)
    for ct = 1:length(expand_idxs)
        for idx = expand_idxs[ct]
            op[idx] += new_op[ct]
        end
    end
end
