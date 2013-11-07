module JuliaSimulator

export evolution_unitary,
       parallel_evolution_unitary,
       unitary_trajectory

using NumericExtensions

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

module JuliaSimulatorTest

export run_sim, sim_setup, setup_test, run_parallel_sim

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

function sim_setup(dimension, 
                   numTimeSteps, 
                   numControls)
    #Create a random natural hamiltonian 
    tmpMat = randn(dimension, dimension) + 1im*randn(dimension, dimension)
    Hnat = tmpMat+tmpMat'

    #Create random control Hamiltonians
    controlHams = Array(Matrix{Complex128}, numControls)
    for ct = 1:numControls
        tmpMat[:] = randn(dimension, dimension) + 1im*randn(dimension, dimension)
        controlHams[ct] = tmpMat+tmpMat'
    end
    #Create random controlfields
    controlFields = randn(numControls, numTimeSteps)
        
    #Control frequencies
    controlFreqs = randn(numControls)

    return Hnat, controlHams, controlFields, controlFreqs

end

function run_sim(Hnat, controlHams, controlFields, controlFreqs)
    evolution_unitary(Hnat, controlHams, controlFields, controlFreqs)
end

function run_parallel_sim(Hnat, controlHams, controlFields, controlFreqs)
    parallel_evolution_unitary(Hnat, controlHams, controlFields, controlFreqs)
end

function setup_test()
  sim_setup(16,2000,4)
end

end # JuliaSimulator.Test

end # JuliaSimulator
