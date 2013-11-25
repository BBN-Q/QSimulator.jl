export FFTransmon,
       Resonator,
       TunableTransmon,
       Qubit

#Resonator 
type Resonator <: QSystem
    label::String
    freq::Float64
    dim::Int
end 
hamiltonian(r::Resonator) = 2*pi*r.freq*number(r)
hamiltonian(r::Resonator, t::Float64) = hamiltonian(r)

#Duffing approximation to transmons

#Fixed frequency transmon 
type FFTransmon <: QSystem
    label::String
    E_C::Float64
    E_J::Float64
    dim::Int
end 

#Tunable transmon 
type TunableTransmon <: QSystem
    label::String
    E_C::Float64
    E_J::Float64 #sum of junction E_J's
    d::Float64 #asymmetry parameter
    dim::Int
    fluxBias::Float64 # flux bias in units of Phi_0
    flux::Float64 #total flux in units of Phi_0
end

#Helper function to calculate effective EJ for a transmon
scale_EJ(flux::Float64, d::Float64) = cos(pi*flux)*sqrt(1 + d^2*(tan(pi*flux)^2))

function hamiltonian(tt::TunableTransmon, t::Float64=0.0)
    myE_J = tt.E_J*scale_EJ(tt.flux, tt.d)
    return 2*pi*(sqrt(8*tt.E_C*myE_J)*number(tt) - (1.0/12)*tt.E_C*(X(tt)^4))
end

#Basic two-level qubit
type Qubit <: QSystem
    label::String
    freq::Float64
end
dim(q::Qubit) = 2
hamiltonian(q::Qubit) = 2*pi*q.freq*number(q)
hamiltonian(q::Qubit, t::Float64) = hamiltonian(q)
