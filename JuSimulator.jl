using NumericExtensions

function expm_eigen(A::Matrix, t)
    #Calculates expm(t*A) via eigenvalue decomposition and assuming Hermitian matrix
    F = eigfact(Hermitian(A))

    # V * diagm(exp(t*D)) * V'
    return scale(F[:vectors], exp(t*F[:values])) * F[:vectors]'
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
