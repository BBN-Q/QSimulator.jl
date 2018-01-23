export ## Types,
       Dipole,
       FlipFlop,
       FluxTransmon,
       SemiClassicalDipole,
       RotatingSemiClassicalDipole,
       ## Methods
       update_params

#Iteractions linearly add new Hamiltonians to the system

type Dipole <: Interaction
    system1::QSystem
    system2::QSystem
    strength::Float64
end

function hamiltonian(f::Dipole)
    return f.strength*kron(X(f.system1), X(f.system2))
end

type FlipFlop <: Interaction
    system1::QSystem
    system2::QSystem
    strength::Float64
end
function hamiltonian(f::FlipFlop)
    return f.strength*(kron(raising(f.system1), lowering(f.system2)) + kron(lowering(f.system1), raising(f.system2)))
end

type SemiClassicalDipole <: Interaction
    system1::Field
    system2::QSystem
    strength::Float64
end
function hamiltonian(scd::SemiClassicalDipole, t::Float64)
    return scd.strength*amplitude(scd.system1, t)*X(scd.system2)
end

type RotatingSemiClassicalDipole <: Interaction
    field::Field
    system2::QSystem
    strength::Float64
end
function hamiltonian(scd::RotatingSemiClassicalDipole, t::Float64)
    saved_phase = scd.field.control.phase
    H = 0.5 * scd.strength * amplitude(scd.field, t) * X(scd.system2)
    scd.field.control.phase += pi/2
    H += 0.5 * scd.strength * amplitude(scd.field, t) * Y(scd.system2)
    scd.field.control.phase = saved_phase
    return H
end

# ParametricInteractions are interactions which modify
# parameters of systems rather than linearly add new Hamiltonians

#Flux control
type FluxTransmon <: ParametricInteraction
    flux::Field
    transmon::AbstractString #refer to transmon by label in parent CompositeQSystem to avoid broken reference on copy
    strength::Float64
end
FluxTransmon(flux::Field, transmon::TunableTransmon, strength::Float64) = FluxTransmon(flux, label(transmon), strength)

#Set the flux in the associated transmon
function update_params(sys::CompositeQSystem, ft::FluxTransmon, t::Float64)
    myTransmon = sys[ft.transmon]
    myTransmon.flux = ft.strength*amplitude(ft.flux,t)
end
