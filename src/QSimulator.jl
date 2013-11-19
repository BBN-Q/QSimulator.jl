# module QSimulator

# export evolution_unitary,
#        parallel_evolution_unitary,
#        unitary_trajectory

using NumericExtensions

abstract QSystem

#Resonator 
type Resonator <: QSystem
    name::String
    freq::Float64
    dim::Int
end 
#Transmon
#Duffing approximation to a transmon 
type Transmon <: QSystem
    name::String
    E_J::Float64
    E_C::Float64
    d::Float64 = 0 #asymmetry 
    dim::Int
end 

for t = [Resonator, Transmon]
    @eval begin
    raising(x::$t) = diagm(sqrt(1:(x.dim-1)), -1)
    lowering(x::$t) = diagm(sqrt(1:(x.dim-1)), 1)
    number(x::$t) = raising(x) * lowering(x)
    X(x::$t) = raising(x) + lowering(x)
    end
end

#Basic two-level qubit
type Qubit <: QSystem
    name::String
    freq::Float64
end
Qubit(name, freq) = Qubit(name, freq, dim)
raising(q::Qubit) = Float64[0 0;1 0]
lowering(q::Qubit) = Float64[0 1;0 0]

for t = [Resonator, Transmon, Qubit]
    @eval begin
        name(x::$t) = x.name
    end
end

for t = [Resonator, Transmon]
    @eval begin
        dim(x::$t) = x.dim
    end
end

type Field
end
#TODO: somehow pull in amplitude from controls
amplitude(f::Field, t::Float64) = 0.0

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
    return scd.strength*amplitude(system1, t)*X(system2)
end

type CompositeQSystem
    subSystems::Dict{String, QSystem}
    interactions::Vector
end

function +(c::CompositeQSystem, q::QSystem)
    c.subSystems[name(q)] = q
    return c
end

function +(c::CompositeQSystem, i::Interaction)
    append!(c.interactions, i)
    return c
end

dim(c::CompositeQSystem) = prod([dim(s) for s in values(c.subSystems)])

function hamiltonian(c::CompositeQSystem, t::Float64)
    #Put together subsystem Hamiltonians
    Htot = zeros(Complex128, dim(c), dim(c))

    for s in c.subSystems
        add!(Htot, expand(hamiltonian(s, t), [])
    end

    #Add interactions

end

function expand(m::Matrix, actingOn::Vector, dims::Vector)
    #Expand a bipartite operator onto a larger Hilbert space
    # m: matrix form of bipartite operator
    # actingOn: tuple of which subsystem index the operator should be acting on 
    # dims: tuple of dimensions of all the subsystems

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
    M = permutedims(M, tuple([reversePerm, reversePerm+l]...))

    #Reshape back
    return reshape(M, prod(dims), prod(dims))
end

#Abstract type for control.   
abstract QControl

#Microwave control 
#Of the form e(t)*cos(2pi*f*t)*Ham
type uWControl <: QControl
    freq::Float64
    Ham::Matrix{Complex128}
end 

#Flux control 
#Complex non-linear control of E_J
type FluxControl <: QControl
end






function expm_eigen(A::Matrix, t)
    #Calculates expm(t*A) via eigenvalue decomposition and assuming Hermitian matrix
    F = eigfact(Hermitian(A))

    # V * diagm(exp(t*D)) * V'
    return scale(F[:vectors], exp(t*F[:values])) * F[:vectors]'
end

function evolution_unitary(Hnat::Matrix{Complex128}, 
                           ampControlHams, 
                           ampControlFields::Matrix{Float64}, 
                           ampControlFreqs::Vector{Float64})

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

function evolution_unitary(controlI, controlQ, calScale)

    const timeStep = 1.0
    Hnat = eye(Complex128, 2)
    tmpH = similar(Hnat)
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

# module QSimulatorTest

# export run_sim, sim_setup, setup_test, run_parallel_sim

# function run_sim(awgdata, chI, chQ, calScale)
#     results = Float64[]
#     initialRho = [1. 0.; 0. 0.]
#     const measOp = [1. 0.; 0. -1.]
#     for ct = 1:length(awgdata[chI])
#         U = evolution_unitary(awgdata[chI][ct], awgdata[chQ][ct], calScale)
#         finalRho = U * initialRho * U'
#         push!(results, real(trace(finalRho * measOp)))
#     end
#     return results
# end

# function sim_setup(dimension, 
#                    numTimeSteps, 
#                    numControls)
#     #Create a random natural hamiltonian 
#     tmpMat = randn(dimension, dimension) + 1im*randn(dimension, dimension)
#     Hnat = tmpMat+tmpMat'

#     #Create random control Hamiltonians
#     controlHams = Array(Matrix{Complex128}, numControls)
#     for ct = 1:numControls
#         tmpMat[:] = randn(dimension, dimension) + 1im*randn(dimension, dimension)
#         controlHams[ct] = tmpMat+tmpMat'
#     end
#     #Create random controlfields
#     controlFields = randn(numControls, numTimeSteps)
        
#     #Control frequencies
#     controlFreqs = randn(numControls)

#     return Hnat, controlHams, controlFields, controlFreqs

# end

# function run_sim(Hnat, controlHams, controlFields, controlFreqs)
#     evolution_unitary(Hnat, controlHams, controlFields, controlFreqs)
# end

# function run_parallel_sim(Hnat, controlHams, controlFields, controlFreqs)
#     parallel_evolution_unitary(Hnat, controlHams, controlFields, controlFreqs)
# end

# function setup_test()
#   sim_setup(16,2000,4)
# end

# end # QSimulator.Test

# end # QSimulator
