module QSimulator

export Resonator,
        Transmon,
        Qubit,
        QuadratureControl,
        load_sequence!,

        Field,

        FlipFlop,
        SemiClassicalDipole,

        hamiltonian,
        expand,

       unitary_propagator,
       parallel_evolution_unitary,
       unitary_trajectory,
       run_sim

using NumericExtensions
using Grid

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
end
CompositeQSystem() = CompositeQSystem(QSystem[], Interaction[])

function +(c::CompositeQSystem, q::QSystem)
    append!(c.subSystems, [q])
    return c
end

function +(c::CompositeQSystem, i::Interaction)
    append!(c.interactions, [i])
    return c
end

labels(c::CompositeQSystem) = [label(s) for s in c.subSystems]
dims(c::CompositeQSystem) = [dim(s) for s in c.subSystems]
dim(c::CompositeQSystem) = prod([dim(s) for s in c.subSystems])

function find_subsystem_pos(c::CompositeQSystem, label::String)
    @assert label in labels(c) "Oops! subsystem not found."
    findin(labels(c), [label])
end


function hamiltonian(c::CompositeQSystem, t::Float64)
    #Put together subsystem Hamiltonians
    Htot = zeros(Complex128, dim(c), dim(c))

    for s in c.subSystems
        add!(Htot, expand(hamiltonian(s, t), [find_subsystem_pos(c, label(s))], dims(c)))
    end

    #Add interactions
    for i in c.interactions
        #Field-system interactions are one-body terms
        if isa(i, SemiClassicalDipole)
            add!(Htot, expand(hamiltonian(i,t), [find_subsystem_pos(c, label(i.system2))], dims(c)))
        #Other interactions are two body terms
        else
            add!(Htot, expand(hamiltonian(i,t), [find_subsystem_pos(c, label(i.system1)), find_subsystem_pos(c, label(i.system2))], dims(c)))
        end
    end

    return Htot
end

function expand(m::Matrix, actingOn::Vector, dims::Vector)
    #Expand an operator onto a larger Hilbert space
    # m: matrix form of  operator
    # actingOn: array of which subsystem index the operator should be acting on 
    # dims: array of dimensions of all the subsystems

    @assert size(m, 1) == prod(dims[actingOn]) "Oops! Dimensions of matrix do not match dims argument."

    #Create the large matrix by tensoring on identity
    l = length(dims)
    indices = [1:l]
    eyeIndices = filter(x->!(x in actingOn), indices)
    M = kron(m, eye(prod(dims[eyeIndices])))  

    #Reshape into multi-dimensional array given by subsystem dimensions
    #Since we have a matrix we repeat for rows then columns
    M = reshape(M, tuple([dims, dims]...))

    #Permute magic 
    forwardPerm = [actingOn, eyeIndices]
    reversePerm = invperm(forwardPerm)
    #Handle the way tensor product indices work (last subsystem in fastest)
    reversePerm = reverse(l+1-reversePerm)
    M = permutedims(M, tuple([reversePerm, reversePerm+l]...))

    #Reshape back
    return reshape(M, prod(dims), prod(dims))
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
