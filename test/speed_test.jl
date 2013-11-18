
function speed_test(repeats)
	Hnat, controlHams, controlFields, controlFreqs = sim_setup(16, 2000, 4)

	# run once for JIT
	evolution_unitary(Hnat, controlHams, controlFields, controlFreqs)
	parallel_evolution_unitary(Hnat, controlHams, controlFields, controlFreqs)

	# then time each
	tserial = @elapsed for ct = 1:repeats
		evolution_unitary(Hnat, controlHams, controlFields, controlFreqs)
	end

	tparallel = @elapsed for ct = 1:repeats
		parallel_evolution_unitary(Hnat, controlHams, controlFields, controlFreqs)
	end

	println("Serial execution in $(tserial/repeats) seconds")
	println("Parallel execution in $(tparallel/repeats) seconds")

	tserial, tparallel
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
