module QSimulator

export Resonator,
        Transmon,
        Qubit,
        CompositeQSystem,

        QuadratureControl,
        load_sequence!,

        Field,

        FlipFlop,
        SemiClassicalDipole,

        hamiltonian,
        expand,

       unitary_propagator,
       evolution_unitary,
       parallel_evolution_unitary,
       unitary_trajectory,
       run_sim

using NumericExtensions
using Grid
using Debug

abstract QSystem

label(q::QSystem) = q.label
dim(q::QSystem) = q.dim
raising(q::QSystem) = diagm(sqrt(1:(dim(q)-1)), -1)
lowering(q::QSystem) = diagm(sqrt(1:(dim(q)-1)), 1)
number(q::QSystem) = raising(q) * lowering(q)
X(q::QSystem) = raising(q) + lowering(q)

#Resonator 
type Resonator <: QSystem
    label::String
    freq::Float64
    dim::Int
end 
#Transmon
#Duffing approximation to a transmon 
type Transmon <: QSystem
    label::String
    E_J::Float64
    E_C::Float64
    d::Float64 = 0 #asymmetry 
    dim::Int
end 

#Basic two-level qubit
type Qubit <: QSystem
    label::String
    freq::Float64
end
dim(q::Qubit) = 2

#System hamiltonians 
hamiltonian(q::Qubit) = 2*pi*q.freq*number(q)
hamiltonian(q::Qubit, t::Float64) = hamiltonian(q)
hamiltonian(r::Resonator) = 2*pi*r.freq*number(r)
hamiltonian(r::Resonator, t::Float64) = hamiltonian(r)

#AWG channels are controls
abstract Control
label(c::Control) = c.label

#A pair of AWG channels driving an IQ mixer with a microwave source at a given frequency
type QuadratureControl <: Control
    label::String
    freq::Float64
    phase::Float64
    timeStep::Float64
    sequence_I::Union(InterpGrid, Nothing)
    sequence_Q::Union(InterpGrid, Nothing)
end
QuadratureControl(label, freq, phase=0., timeStep=1/1.2, sequence_I=nothing, sequence_Q=nothing) = QuadratureControl(label, freq, phase, timeStep, sequence_I, sequence_Q)

#Create the interpolated object
function load_sequence!(qc::QuadratureControl, seqDict::Dict, n::Int)
    #For quadrature controls we need two channels
    #Hack around off the expected BBNAPSx-xx style
    chan_I = label(qc)[end-1]
    chan_Q = label(qc)[end]
    APSLabel = label(qc)[1:7]
    qc.sequence_I = InterpGrid(seqDict[APSLabel]["chan_"*string(chan_I)][n], BCnearest, InterpNearest)
    qc.sequence_Q = InterpGrid(seqDict[APSLabel]["chan_"*string(chan_Q)][n], BCnearest, InterpNearest)
    return nothing
end
function amplitude(qc::QuadratureControl, t::Float64)
    #Create the complex I/Q pair
    #Scale time by the timestep
    phasor = qc.sequence_I[t/qc.timeStep] + 1im*qc.sequence_Q[t/qc.timeStep]
    return abs(phasor)*cos(2*pi*qc.freq*t + qc.phase + angle(phasor))
end

#Flux control 
#Complex non-linear control of E_J
type FluxControl <: Control
end

#AWG controls get connected, possibly through a transfer function, to fields
type Field
    control::Control
end
#Pull in amplitude from controls
#TODO: allow for transfer function
amplitude(f::Field, t::Float64) = amplitude(f.control, t)

abstract Interaction

type FlipFlop <: Interaction
    system1::QSystem
    system2::QSystem
    strength::Float64
end
function hamiltonian(f::FlipFlop, t::Float64=0.0)
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

type CompositeQSystem
    #An OrderedDict would be ideal here
    subSystems::Vector{QSystem}
    interactions::Vector{Interaction}
    subSystemExpansions::Vector{Vector{Vector{Int}}}
    interactionExpansions::Vector{Vector{Vector{Int}}}
end
CompositeQSystem() = CompositeQSystem(QSystem[], Interaction[], Vector{Vector{Int}}[], Vector{Vector{Int}}[])

function +(c::CompositeQSystem, q::QSystem)
    append!(c.subSystems, [q])
    update_expansion_indices!(c)
    return c
end

function +(c::CompositeQSystem, i::Interaction)
    append!(c.interactions, [i])
    update_expansion_indices!(c)
    return c
end

function update_expansion_indices!(c::CompositeQSystem)
    c.subSystemExpansions = [Vector{Int}[] for _ = 1:length(c.subSystems)]
    for (ct, sys) in enumerate(c.subSystems)
        c.subSystemExpansions[ct] = expand_indices([ct], dims(c))
    end

    c.interactionExpansions = [Vector{Int}[] for _ = 1:length(c.interactions)]
    for (ct, i) in enumerate(c.interactions)
        c.interactionExpansions[ct] = expand_indices(find_subsystem_pos(c, i), dims(c))
    end
end

labels(c::CompositeQSystem) = [label(s) for s in c.subSystems]
dims(c::CompositeQSystem) = [dim(s) for s in c.subSystems]
dim(c::CompositeQSystem) = prod([dim(s) for s in c.subSystems])

function find_subsystem_pos(c::CompositeQSystem, s::QSystem)
    @assert s in c.subSystems "Oops! Subsystem not found."
    findin(c.subSystems, [s])
end

function find_subsystem_pos(c::CompositeQSystem, i::Interaction)
    #Field-system interactions are one-body terms
    if isa(i, SemiClassicalDipole)
        return find_subsystem_pos(c, i.system2)
    else
        return [find_subsystem_pos(c, i.system1), find_subsystem_pos(c, i.system2)]
    end
end

function hamiltonian(c::CompositeQSystem, t::Float64)

    #Initialize Hamiltonian
    Htot = zeros(Complex128, dim(c), dim(c))

    #Add in all the terms    
    hamiltonian_add!(Htot, c, t)

    return Htot
end

function hamiltonian_add!(Ham::Matrix{Complex128}, c::CompositeQSystem, t::Float64)
    #Fast system hamiltonian calculator with total Hamiltonian preallocated
    
    #Zero the Hamiltonian memory
    Ham[:] = zero(Complex128)

    #Add together subsystem Hamiltonians
    for (ct, s) in enumerate(c.subSystems)
        expand_add!(Ham, hamiltonian(s, t), c.subSystemExpansions[ct])
    end

    #Add interactions
    for (ct, i) in enumerate(c.interactions)
        expand_add!(Ham, hamiltonian(i,t), c.interactionExpansions[ct])
    end
end

function expand(m::Matrix, actingOn::Vector, dims::Vector)
    #Expand an operator onto a larger Hilbert space
    # m: matrix form of  operator
    # actingOn: array of which subsystem index the operator should be acting on 
    # dims: array of dimensions of all the subsystems

    @assert size(m, 1) == prod(dims[actingOn]) "Oops! Dimensions of matrix do not match dims argument."

    #Create the large matrix by tensoring on identity
    l = length(dims)
    eyeIndices = filter(x->!(x in actingOn), 1:l)
    M = kron(m, eye(eltype(m), prod(dims[eyeIndices])))  

    #Reshape into multi-dimensional array given by subsystem dimensions
    #Since we have a matrix we repeat for rows then columns
    M = reshape(M, tuple([dims, dims]...))

    #Permute magic 
    forwardPerm = [actingOn, eyeIndices]
    reversePerm = invperm(forwardPerm)
    #Handle the way tensor product indices work (last subsystem is fastest)
    reversePerm = reverse(l+1-reversePerm)
    M = permutedims(M, tuple([reversePerm, reversePerm+l]...))

    #Reshape back
    return reshape(M, prod(dims), prod(dims))
end

function expand(m::Matrix, indices::Vector, sizeM::Int )
    M = zeros(eltype(m), (sizeM, sizeM))
    for (ct, inds) in enumerate(indices)
        M[inds] = m[ct]
    end
    return M
end

function expand_add!(M::Matrix, m::Matrix, indices::Vector )
    #Add to certain indices of M with terms from m according to expansion indices.
    for (ct, inds) in enumerate(indices)
        M[inds] += m[ct]
    end
end

function expand_indices(actingOn::Vector, dims::Vector)
    #Calculate the indices for expansion
    actingOnDim = prod(dims[actingOn])
    sm = (actingOnDim, actingOnDim)
    lenm = actingOnDim^2;
    M = expand(reshape(1:lenm, sm), actingOn, dims)
    return [find(M .== x) for x in 1:lenm]
end


function expm_eigen(A::Matrix, t)
    #Calculates expm(t*A) via eigenvalue decomposition and assuming Hermitian matrix
    F = eigfact(Hermitian(A))

    # V * diagm(exp(t*D)) * V'
    return scale(F[:vectors], exp(t*F[:values])) * F[:vectors]'
end

function evolution_unitary(Hnat::Matrix{Complex128}, 
                           controlHams, 
                           controlFields::Matrix{Float64}, 
                           controlFreqs::Vector{Float64})

    const timeStep = 0.01
    Uprop = eye(Complex128, size(Hnat,1))
    tmpH = similar(Hnat)
    localH = similar(Hnat)

    for timect = 1:size(controlFields,2)
        tmpH[:] = Hnat
        for controlct = 1:size(controlFields,1)
            localH[:] = controlHams[controlct]
            multiply!(localH, controlFields[controlct, timect]*cos(2*pi*timeStep*timect*controlFreqs[controlct]))
            add!(tmpH, localH)
        end
        Uprop *= expm_eigen(tmpH, 1im*2*pi*timeStep)
    end

    return Uprop'
end


function parallel_evolution_unitary(Hnat::Matrix{Complex128}, 
                                    controlHams, 
                                    controlFields::Matrix{Float64}, 
                                    controlFreqs::Vector{Float64})

    const timeStep = 0.01
    tmpH = similar(Hnat)
    localH = similar(Hnat)

    Uprop = @parallel (*) for timect = 1:size(controlFields,2)
        tmpH[:] = Hnat
        for controlct = 1:size(controlFields,1)
            localH[:] = controlHams[controlct]
            multiply!(localH, controlFields[controlct, timect]*cos(2*pi*timeStep*timect*controlFreqs[controlct]))
            add!(tmpH, localH)
        end
        expm_eigen(tmpH, 1im*2*pi*timeStep)
    end

    return Uprop'
end

function unitary_propagator(sys::CompositeQSystem, timeStep::Float64, endTime::Float64)

    #Preallocate Hamiltonian memory
    Ham = zeros(Complex128, (dim(sys), dim(sys)))
    Uprop = @parallel (*) for ct = 1:fld(endTime, timeStep)
        #a *= b expands to a = a*b
        expm_eigen(hamiltonian(sys, ct*timeStep), 1im*timeStep)
    end
    return Uprop'
end

function evolution_unitary(controlI, controlQ, calScale)

    const timeStep = 1.0
    const sx = 0.5*Complex128[0 1; 1 0]
    const sy = 0.5*Complex128[0 -1im; 1im 0]

    Uprop = @parallel (*) for timect = 1:length(controlI)
        expm_eigen(controlI[timect]*sx + controlQ[timect]*sy, 1im*2*pi*timeStep*calScale)
    end

    return Uprop'
end

function unitary_trajectory(controlI, controlQ, calScale)
    const timeStep = 1.0
    const sx = 0.5*Complex128[0 1; 1 0]
    const sy = 0.5*Complex128[0 -1im; 1im 0]
    const inputZ = Complex128[1, 0]
    const inputX = Complex128[1, 1]/sqrt(2)

    Uprop = eye(Complex128, 2)
    finalZ = Array(Vector{Complex128}, length(controlI))
    finalX = Array(Vector{Complex128}, length(controlI))

    for timect = 1:length(controlI)
        Uprop *= expm_eigen(controlI[timect]*sx + controlQ[timect]*sy, 1im*2*pi*timeStep*calScale)
        finalZ[timect] = Uprop' * inputZ
        finalX[timect] = Uprop' * inputX
    end

    return finalZ, finalX
end

function run_sim(awgdata, chI, chQ, calScale)
    results = Float64[]
    initialRho = [1. 0.; 0. 0.]
    const measOp = [1. 0.; 0. -1.]
    for ct = 1:length(awgdata[chI])
        U = evolution_unitary(awgdata[chI][ct], awgdata[chQ][ct], calScale)
        finalRho = U * initialRho * U'
        push!(results, real(trace(finalRho * measOp)))
    end
    return results
end

end # QSimulator
