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
hamiltonian(r::Resonator) = r.freq*number(r)

#Duffing spectrum
function duffing_spectrum(dim::Int, E_J::Float64, E_C::Float64)
    plasma = sqrt(8*E_J*E_C)
    Float64[ (plasma-E_C/2)*j - E_C/2*j^2 for j in 0:(dim-1) ]
end

#Fixed frequency transmon 
type FFTransmon <: QSystem
    label::String
    E_C::Float64
    E_J::Float64
    dim::Int
    hamiltonian::Matrix{Float64}
    #Inner constructor needed to make sure the hamiltonian is properly innitalized
    function FFTransmon(label::String, E_C::Float64, E_J::Float64, dim::Int, spectrum::Function=duffing_spectrum)
        new(label, E_C, E_J, dim, diagm(spectrum(dim,E_J,E_C)))
    end
end 

function hamiltonian(tt::FFTransmon, t::Float64=0.0)
    return tt.hamiltonian
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
    spectrum::Function # function used to approximate the transmon spectrum
end

TunableTransmon(label::String, 
                E_C::Float64, 
                E_J::Float64, 
                d::Float64, 
                dim::Int, 
                fluxBias::Float64,
                spectrum::Function=duffing_spectrum) = TunableTransmon(label, E_C, E_J, d, dim, fluxBias, fluxBias, spectrum)

#Helper function to calculate effective EJ for a transmon
scale_EJ(flux::Float64, d::Float64) = cos(pi*flux)*sqrt(1 + d^2*(tan(pi*flux)^2))

function hamiltonian(tt::TunableTransmon, t::Float64=0.0)
    myE_J = tt.E_J*scale_EJ(tt.flux, tt.d)
    return diagm(tt.spectrum(tt.dim,myE_J,tt.E_C))
end

#Basic two-level qubit
type Qubit <: QSystem
    label::String
    freq::Float64
end
dim(q::Qubit) = 2
hamiltonian(q::Qubit) = q.freq*number(q)
