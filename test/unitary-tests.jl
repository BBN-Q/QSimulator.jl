using Test
using QSimulator
import QSimulator: hamiltonian

function test_1()

  for i=1:100
      r1 = Resonator("r1", 1.0, 10)

      system = CompositeQSystem()

      system += r1

      state = randn(10) + 1im*randn(10)
      state = state/LinearAlgebra.norm(state,2)
      t = 1.0+10*abs(randn())
      #t = abs(randn())

      ans = state .* exp.(-1im*2*pi*t*Float64[k for k in 0:9])
      ev_state = unitary_evolution(state,system,.1,0.,t)

      @test isapprox(dot(ans,ev_state), 1.0, atol = 1e-12)
      #@test isapprox(norm(state, 2) 1.0 1e-12)
  end
end

function test_drive()

    # create the qubits
    q1 = Qubit("q1", 4.0)
    q2 = Qubit("q2", 4.5)

    # create a time dependant drive
    drive = MicrowaveControl("drive1", 4.0, timeStep = 1.0)

    # create the Hamiltonian with 'objects' and 'interactions'
    sys = q1 + q2 + FlipFlop(q1, q2, 0.005) + SemiClassicalDipole(Field(drive), q1, 0.01);

    seq = ones(200)
    load_sequence!(drive, seq);

    # test starting eigenvals
    E_test = LinearAlgebra.eigvals(hamiltonian(sys, 0.0))
    E = [-2.4999878473470395e-5,
            3.999975000433993,
            4.500024999566007,
            8.500024999878473]
    @test isapprox(E, E_test, atol = 1e-12)

    # define the initial state
    rhoIn = zeros(ComplexF64, dim(sys), dim(sys))
    rhoIn[1] = 1;
    #define the measurement projector
    evecs = LinearAlgebra.eigvecs(hamiltonian(sys))
    measOp = evecs[:,2] * evecs[:,2]'
    real(measOp)

    function run_simulation(sys, rhoIn, measOp, times)
    populations = zeros(length(times)-1)
    U = QSimulator.eye(ComplexF64, dim(sys))
    dt = 0.01
        for ct = 1:length(times)-1
            U = unitary_propagator(sys, dt, 0.0, times[ct+1])
            populations[ct] = real(LinearAlgebra.tr(U*rhoIn*U'*measOp))
        end
    return populations
    end

    times = 0.0:20.0
    populations_test = run_simulation(sys, rhoIn, measOp, times);
    populations = [0.001007229516897708,
                     0.003966998572908745,
                     0.008873938406129277,
                     0.01570878588050283,
                     0.024444712566132513,
                     0.03504742273249326,
                     0.047475299645262155,
                     0.06167955185968114,
                     0.0776044282644603,
                     0.09518740596921914,
                     0.11435947513329113,
                     0.1350453610888861,
                     0.15716387889776454,
                     0.1806281818148819,
                     0.20534618442196084,
                     0.23122082912331615,
                     0.2581505756295072,
                     0.28602967689427616,
                     0.31474873168055584,
                     0.3441949611994505]

    @test isapprox(populations, populations_test, atol = 1e-12)
end

@testset "unitary" begin
    test_1()
end

@testset "driven" begin
    test_drive()
end
