export ## Methods
       unitary_propagator


function expm_eigen(A::Matrix, t)
    #Calculates expm(t*A) via eigenvalue decomposition and assuming Hermitian matrix
    F = eigfact(Hermitian(A))

    # V * diagm(exp(t*D)) * V'
    return scale(F[:vectors], exp(t*F[:values])) * F[:vectors]'
end

function unitary_propagator(sys::CompositeQSystem, timeStep::Float64, startTime::Float64, endTime::Float64)

    #Preallocate Hamiltonian memory
    Ham = zeros(Complex128, (dim(sys), dim(sys)))
    Uprop = @parallel (*) for time = startTime:timeStep:endTime
        #a *= b expands to a = a*b
        hamiltonian_add!(Ham, sys, time)
        expm_eigen(Ham, 1im*timeStep)
    end
    return Uprop'
end


