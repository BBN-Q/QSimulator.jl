export ## Types,
       FlipFlop,
       FluxTransmon,
       SemiClassicalDipole,
       ## Methods
       update_params

#Iteractions linearly add new Hamiltonians to the system
type FlipFlop <: Interaction
    system1::QSystem
    system2::QSystem
    strength::Float64
end
function hamiltonian(f::FlipFlop, t::Float64=0.0)
    return 2pi*f.strength*(kron(raising(f.system1), lowering(f.system2)) + kron(lowering(f.system1), raising(f.system2)))
end

type SemiClassicalDipole <: Interaction
    system1::Field
    system2::QSystem
    strength::Float64
end
function hamiltonian(scd::SemiClassicalDipole, t::Float64)
    return scd.strength*amplitude(scd.system1, t)*X(scd.system2)
end

# ParametricInteractions are interactions which modify 
# parameters of systems rather than linearly add new Hamiltonians

#Flux control 
type FluxTransmon <: ParametricInteraction
    flux::Field
    transmon::TunableTransmon
    strength::Float64
end
#Set the flux in the associated transmon
function update_params(ft::FluxTransmon, t::Float64)
    ft.transmon.flux = ft.transmon.fluxBias + ft.strength*amplitude(ft.flux,t)
end

