import Optim

export ## Types
       Duffing,
       Resonator,
       TunableDuffingTransmon,
       TunableFullTransmon,
       TunableTransmon,
       Qubit,
       ## Methods
       hamiltonian

#Resonator
mutable struct Resonator <: QSystem
    label::AbstractString
    freq::Float64
    dim::Int
end
hamiltonian(r::Resonator) = r.freq*number(r)

#Fixed frequency transmon
mutable struct Transmon <: QSystem
    label::AbstractString
    E_C::Float64
    E_J::Float64
    dim::Int
end

"""
Transmon(label::AbstractString, ν::Float64, α::Float64; dim=3)

#Arguments
* ν - 0 ↔ 1 frequency
* α - anharmonicity
"""
# function Transmon(label::AbstractString, ν::Float64, α::Float64; dim=3)
#     E_C = abs(α)
#     E_J = (ν + E_C)^2 / (8 * E_C)
#     return Transmon(label, E_C, E_J, dim)
# end

# hamiltonian(t::Transmon) = sqrt(8*t.E_J*t.E_C)*number(t) - t.E_C/12*(X(t)^4)
function hamiltonian(t::Transmon)
  N = floor(Int, dim(t)/2)
  4 * t.E_C * diagm(0 => (-N:N).^2) - t.E_J  * 0.5 * (diagm(-1 => ones(dim(t)-1)) + diagm(1 => ones(dim(t)-1)))
end

"""
Fit  E_C and E_J for a fixed frequency transmon given f_01 and anharmonicity
"""
function fit_fixed_transmon(f_01, α, dim)
    # helper function for least squares difference between measured and calculated
    function f(params)
      E_C = params[1]
      E_J  = params[2]
      t = QSimulator.Transmon("dummy", E_C, E_J, dim)
      levels = eigvals(hamiltonian(t))
      test_f_01 = levels[2] - levels[1]
      test_f_12 = levels[3] - levels[2]
      ( test_f_01 - f_01)^2 + ( test_f_12 -  test_f_01 - α)^2
    end
    res = Optim.optimize(f, [abs(α), (f_01+abs(α))^2 / (8 * abs(α))] )
    return res.minimizer[1], res.minimizer[2]
end

#Tunable transmon
abstract type TunableTransmon <: QSystem end

mutable struct TunableFullTransmon <: TunableTransmon
    label::AbstractString
    E_C::Float64
    E_J::Float64 #sum of junction E_J's
    d::Float64 #asymmetry parameter
    dim::Int
    flux::Float64 #total flux in units of Phi_0
end

#Helper function to calculate effective EJ for a transmon
# scale_EJ(flux::Float64, d::Float64) = cos(pi*flux)*sqrt(1 + d^2*(tan(pi*flux)^2))

scale_EJ(flux::Float64, d::Float64) = sqrt(cos(π*flux)^2 + d^2*(sin(π*flux)^2))


function hamiltonian(tt::TunableFullTransmon, t::Float64=0.0)
    myE_J = tt.E_J*scale_EJ(tt.flux, tt.d)
    return (sqrt(8*tt.E_C*myE_J)*number(tt) - (1.0/12)*tt.E_C*(X(tt)^4))
end

#Tunable Duffing approx. transmon
mutable struct TunableDuffingTransmon <: TunableTransmon
    label::AbstractString
    E_C::Float64
    E_J::Float64 #sum of junction E_J's
    d::Float64 #asymmetry parameter
    dim::Int
    flux::Float64 #total flux in units of Phi_0
end

function hamiltonian(tt::TunableDuffingTransmon, t::Float64=0.0)
    myE_J = tt.E_J*scale_EJ(tt.flux, tt.d)
    omega_p = sqrt(8*tt.E_C*myE_J)
    omega = [(omega_p-tt.E_C/2)*i-tt.E_C/2*i^2 for i in 0:(tt.dim-1)]
    return diagm(0 => omega)
end

#Basic two-level qubit
mutable struct Qubit <: QSystem
    label::AbstractString
    freq::Float64
end
dim(q::Qubit) = 2
hamiltonian(q::Qubit) = q.freq*number(q)

#Duffing oscillator
mutable struct Duffing <: QSystem
    label::AbstractString
    freq::Float64
    alpha::Float64
    dim::Int
end
hamiltonian(s::Duffing) = (s.freq - 0.5*s.alpha)*number(s) + 0.5*s.alpha * number(s)^2

"""
Fit E_C, E_J and d for a tunable transmon given f_max, f_min and anharmonicty
"""
function fit_tunable_transmon()
  0
end
