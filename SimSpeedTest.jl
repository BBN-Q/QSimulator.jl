using NumericExtensions

function expm_eigen(A::Matrix, t)
    #Calculates expm(t*A) via eigenvalue decomposition and assuming Hermitian matrix
    F = eigfact(Hermitian(A))

    # V * diagm(exp(t*D)) * V'
    return scale(F[:vectors], exp(t*F[:values])) * F[:vectors]'
end

function evolution_unitary(Hnat::Matrix{Complex128}, controlHams, controlFields::Matrix{Float64}, controlFreqs::Vector{Float64})

    const timeStep = 0.01
    Uprop = eye(Complex128, size(Hnat,1))
    tmpH = similar(Hnat)
    localH = similar(Hnat)

    for timect = 1:size(controlFields,2)
        tmpH[:] = Hnat
        for controlct = 1:size(controlFields,1)
            # tmpH += controlFields[controlct, timect]*cos(2*pi*timeStep*timect*controlFreqs[controlct])*controlHams[:, :, controlct]
            localH[:] = controlHams[controlct]
            multiply!(localH, controlFields[controlct, timect]*cos(2*pi*timeStep*timect*controlFreqs[controlct]))
            add!(tmpH, localH)
        end
        Uprop *= expm_eigen(tmpH, -1im*2*pi*timeStep)
    end

    return Uprop
end

function parallel_evolution_unitary(Hnat::Matrix{Complex128}, controlHams, controlFields::Matrix{Float64}, controlFreqs::Vector{Float64})

    const timeStep = 0.01
    tmpH = similar(Hnat)
    localH = similar(Hnat)

    Uprop = @parallel (*) for timect = 1:size(controlFields,2)
        tmpH[:] = Hnat
        for controlct = 1:size(controlFields,1)
            # tmpH += controlFields[controlct, timect]*cos(2*pi*timeStep*timect*controlFreqs[controlct])*controlHams[:, :, controlct]
            localH[:] = controlHams[controlct]
            multiply!(localH, controlFields[controlct, timect]*cos(2*pi*timeStep*timect*controlFreqs[controlct]))
            add!(tmpH, localH)
        end
        expm_eigen(tmpH, -1im*2*pi*timeStep)
    end

    return Uprop
end

function sim_setup(dimension, numTimeSteps, numControls)
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

function run_sim()
    evolution_unitary(Hnat, controlHams, controlFields, controlFreqs)
end

function run_parallel_sim()
    parallel_evolution_unitary(Hnat, controlHams, controlFields, controlFreqs)
end

Hnat, controlHams, controlFields, controlFreqs = sim_setup(16, 2000, 4)
