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
type Resonator <: QSystem
    label::String
    freq::Float64
    dim::Int
end 
hamiltonian(r::Resonator) = r.freq*number(r)

#Fixed frequency transmon 
type Transmon <: QSystem
    label::String
    E_C::Float64
    E_J::Float64
    dim::Int
end

function Transmon(label::String, nu::Float64, alpha::Float64; dim=3)
    E_C = abs(alpha)
    E_J = (nu + E_C)^2 / (8 * E_C)
    return Transmon(label, E_C, E_J, dim)
end

hamiltonian(t::Transmon) = sqrt(8*t.E_J*t.E_C)*number(t) - t.E_C/12*(X(t)^4)

#Tunable transmon 

abstract TunableTransmon <: QSystem

type TunableFullTransmon <: TunableTransmon
    label::String
    E_C::Float64
    E_J::Float64 #sum of junction E_J's
    d::Float64 #asymmetry parameter
    dim::Int
    fluxBias::Float64 # flux bias in units of Phi_0
    flux::Float64 #total flux in units of Phi_0
end
TunableFullTransmon(label::String, E_C::Float64, E_J::Float64, d::Float64, dim::Int, fluxBias::Float64) = 
  TunableFullTransmon(label, E_C, E_J, d, dim, fluxBias, fluxBias)

#Helper function to calculate effective EJ for a transmon
scale_EJ(flux::Float64, d::Float64) = cos(pi*flux)*sqrt(1 + d^2*(tan(pi*flux)^2))

function hamiltonian(tt::TunableFullTransmon, t::Float64=0.0)
    myE_J = tt.E_J*scale_EJ(tt.flux, tt.d)
    return (sqrt(8*tt.E_C*myE_J)*number(tt) - (1.0/12)*tt.E_C*(X(tt)^4))
end

#Tunable Duffing approx. transmon 
type TunableDuffingTransmon <: TunableTransmon
    label::String
    E_C::Float64
    E_J::Float64 #sum of junction E_J's
    d::Float64 #asymmetry parameter
    dim::Int
    fluxBias::Float64 # flux bias in units of Phi_0
    flux::Float64 #total flux in units of Phi_0
end
TunableDuffingTransmon(label::String, E_C::Float64, E_J::Float64, d::Float64, dim::Int, fluxBias::Float64) = TunableDuffingTransmon(label, E_C, E_J, d, dim, fluxBias, fluxBias)

#Helper function to calculate effective EJ for a transmon
scale_EJ(flux::Float64, d::Float64) = cos(pi*flux)*sqrt(1 + d^2*(tan(pi*flux)^2))

function hamiltonian(tt::TunableDuffingTransmon, t::Float64=0.0)
    myE_J = tt.E_J*scale_EJ(tt.flux, tt.d)
    omega_p = sqrt(8*tt.E_C*myE_J)
    omega = [(omega_p-tt.E_C/2)*i-tt.E_C/2*i^2 for i in 0:(tt.dim-1)]
    return diagm(omega)
end

#Basic two-level qubit
type Qubit <: QSystem
    label::String
    freq::Float64
end
dim(q::Qubit) = 2
hamiltonian(q::Qubit) = q.freq*number(q)

#Duffing oscillator
type Duffing <: QSystem
    label::String
    freq::Float64
    alpha::Float64
    dim::Int
end
hamiltonian(s::Duffing) = (s.freq - 0.5*s.alpha)*number(s) + 0.5*s.alpha * number(s)^2
