export ## Methods
       unitary_propagator,
       unitary_evolution,
       liouvillian_propagator,
       liouvillian_evolution

using ExpmV

function expm_eigen(A::Matrix, t)
    #Calculates exp(t*A) via eigenvalue decomposition and assuming Hermitian matrix
    F = eigfact(Hermitian(A))

    return F[:vectors] * Diagonal(exp(t*F[:values])) * F[:vectors]'
end

function unitary_propagator(sys::CompositeQSystem,
                            timeStep::Float64,
                            startTime::Float64,
                            endTime::Float64;
                            parallelize=true)

    #Preallocate Hamiltonian memory
    Ham = zeros(Complex128, (dim(sys), dim(sys)))

    times = startTime:timeStep:(endTime-timeStep)

    if length(times) < 1
        error("Time step is too small")
    end

    if parallelize
        Uprop = @parallel (*) for time = times
            #a *= b expands to a = a*b
            hamiltonian_add!(Ham, sys, time)
            expm_eigen(Ham, 1im*2pi*timeStep)
        end
    else
        Uprop = eye(dim(sys))
        for time = times
            #a *= b expands to a = a*b
            hamiltonian_add!(Ham, sys, time)
            Uprop *= expm_eigen(Ham, 1im*2pi*timeStep)
        end
    end

    if (endTime-times[end]) > timeStep
        hamiltonian_add!(Ham, sys, times[end]+timeStep)
        Uprop *= expm_eigen(Ham, 1im*2pi*(endTime-times[end]-timeStep))
    end

    return Uprop'
end

function unitary_evolution{T<:Number}(state::Vector{T},
                                      sys::CompositeQSystem,
                                      timeStep::Float64,
                                      startTime::Float64,
                                      endTime::Float64)

    state_cp = copy(state)

    #Preallocate Hamiltonian memory
    Ham = spzeros(Complex128, dim(sys), dim(sys))

    times = startTime:timeStep:(endTime-timeStep)

    if length(times) < 1
        error("Time step is too small")
    end

    Uprop = eye(dim(sys))
    for time = times
        hamiltonian_add!(Ham, sys, time)
        state_cp = expmv(-1im*2pi*timeStep,Ham,state_cp)
    end

    if (endTime-times[end]) > timeStep
        hamiltonian_add!(Ham, sys, times[end]+timeStep)
        state_cp = expmv(-1im*2pi*(endTime-times[end]-timeStep),Ham,state_cp)
    end

    return state_cp
end

function liouvillian_propagator(sys::CompositeQSystem,
                                timeStep::Float64,
                                startTime::Float64,
                                endTime::Float64;
                                parallelize=true)

    #Preallocate memory
    liouv = zeros(Complex128, dim(sys)^2, dim(sys)^2)

    times = startTime:timeStep:(endTime-timeStep)

    if length(times) < 1
        error("Time step is too small")
    end

    if parallelize
        Lprop = @parallel (*) for time = times
            liouvillian_dual_add!(liouv, sys, time)
            expm(2pi*timeStep*liouv)
        end
    else
        Lprop = eye(dim(sys))
        for time = times
            liouvillian_dual_add!(liouv, sys, time)
            Lprop *= expm(2pi*timeStep*liouv)
        end
    end

    if (endTime-times[end]) > timeStep
        liouvillian_dual_add!(liouv, sys, times[end]+timeStep)
        Lprop *= expm(2pi*(endTime-times[end]-timeStep)*liouv)
    end

    return Lprop'
end

function liouvillian_evolution{T<:Number}(state::Matrix{T},
                                          sys::CompositeQSystem,
                                          timeStep::Float64,
                                          startTime::Float64,
                                          endTime::Float64)

    state_cp = copy(state[:])

    #Preallocate memory
    liouv = spzeros(Complex128, dim(sys)^2, dim(sys)^2)

    times = startTime:timeStep:(endTime-timeStep)

    if length(times) < 1
        error("Time step is too small")
    end

    Lprop = speye(dim(sys))
    for time = times
        liouvillian_add!(liouv, sys, time)
        state_cp = expmv(2pi*timeStep,liouv,state_cp)
    end

    if (endTime-times[end]) > timeStep
        liouvillian_add!(liouv, sys, times[end]+timeStep)
        state_cp = expmv(2pi*(endTime-times[end]-timeStep),liouv,state_cp)
    end

    return reshape(state_cp,size(state,1),size(state,2))
end
