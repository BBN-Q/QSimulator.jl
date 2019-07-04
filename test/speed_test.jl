using QSimulator

function speed_test(repeats)
	T_f = 1000.0
	sys = sim_setup(4, T_f)
	dt = 1.0

	# run once for JIT
	unitary_propagator(sys, dt, 0.0, T_f, parallelize=false)
	unitary_propagator(sys, dt, 0.0, T_f, parallelize=true)

	# then time each
	tserial = @elapsed for ct = 1:repeats
		unitary_propagator(sys, dt, 0.0, T_f, parallelize=false)
	end

	tparallel = @elapsed for ct = 1:repeats
		unitary_propagator(sys, dt, 0.0, T_f, parallelize=true)
	end

	println("Serial execution in $(tserial/repeats) seconds")
	println("Parallel execution in $(tparallel/repeats) seconds")

	@test isapprox(abs(tserial/repeats - tparallel/repeats), 0.0, atol = 1e-1)
end

function sim_setup(dimension, numTimeSteps)
	α = -0.350
	J = 0.002
	Ω = 0.010
	Q1 = Duffing("q1", 0.0, α, dimension)
	Q2 = Duffing("q2", 1.0, α, dimension)

	drive = MicrowaveControl("CR", 1.0, timeStep=1.0)
	sys = Q1 + Q2 + FlipFlop(Q1, Q2, J) + RotatingSemiClassicalDipole(Field(drive), Q2, Ω)

	seq = ones(round(Int, numTimeSteps))
	load_sequence!(drive, seq)

	return sys
end

@testset "speed test" begin
    speed_test(5)
end
