using Test
using QSimulator
import QuantumInfo

function test_1()
  q1 = Qubit("q1", 5.0)

  system = CompositeQSystem()

  system += q1

  ss_exp = system.subSystemExpansions
  in_exp = system.interactionExpansions

  system += Cooling("T1", 0.0, q1)

  @test ss_exp == system.subSystemExpansions
  @test in_exp == system.interactionExpansions
  @test system.dissipatorExpansions[1] == (IndexSet[Int[1,11],Int[2,12],Int[5,15],Int[6,16]],
                                           IndexSet[Int[1,6],Int[3,8],Int[9,14],Int[11,16]],
                                           IndexSet[Int[i] for i in 1:16])
end

function test_2()
  q1 = Qubit("q1", 5.0)
  q2 = Qubit("q2", 6.0)

  system = CompositeQSystem()

  system += q1
  system += q2
  system += FlipFlop(q1, q2, 0.0)

  ss_exp = system.subSystemExpansions
  in_exp = system.interactionExpansions

  system += Cooling("T1", 0.0, q1)

  @test ss_exp == system.subSystemExpansions
  @test in_exp == system.interactionExpansions

end

function test_3()
  system = CompositeQSystem()

  @test_throws ErrorException hamiltonian(system,0.0)
end

function test_4()
  sz = Float64[0 0; 0 1];
  sm = Float64[0 1; 0 0];

  q1 = Qubit("q1", 1.0)

  system = CompositeQSystem()

  L = zeros(ComplexF64, 4, 4)

  system += q1

  QSimulator.liouvillian_dual_add!(L, system, 0.0)
  @test isapprox(norm(L-QuantumInfo.hamiltonian(sz)', 2), 0.0)

  system += Cooling("T1", 0.1, q1)
  QSimulator.liouvillian_dual_add!(L, system, 0.0)
  @test isapprox(norm(L-(QuantumInfo.hamiltonian(sz)+.1*QuantumInfo.dissipator(sm))',2), 0.0)
end

function test_5()
  q1 = Qubit("q1", 0.0)

  system = CompositeQSystem()

  system += q1

  system += Cooling("T1", 1.0, q1)

  Lp = liouvillian_propagator(system,100.,0.,1000.)

  @test isapprox(norm(Lp-[1 0 0 1; 0 0 0 0; 0 0 0 0; 0 0 0 0], 2), 0.0, atol=1e-12)
end

function test_6()
  q1 = Qubit("q1", 0.0)

  system = CompositeQSystem()

  system += q1

  system += Cooling("T1", 1.0, q1)

  state = [0. 0.; 0. 1.]

  ev_state = liouvillian_evolution(state,system,100.,0.,1000.)

  @test isapprox(norm(ev_state - [1 0; 0 0], 1), 0.0, atol=1e-12)
end

function test_7()
  q1 = Qubit("q1", 0.0)

  system = CompositeQSystem()

  system += q1

  system += Cooling("T1", 1.0/2pi, q1)

  state = [0. 0.; 0. 1.]

  for i=1:10
      t = 1.0 + 10*abs(randn())
      ev_state = liouvillian_evolution(state,system,.01,0.,t)
      ans = diagm(0 => [1-exp(-t),exp(-t)])

      @test isapprox(norm(ev_state - ans, 1), 0.0, atol = 1e-12)
  end
end


@testset "dissipators 1" begin
    test_1()
end

@testset "dissipators 2" begin
    test_2()
end

@testset "dissipators 3" begin
    test_3()
end

@testset "dissipators 4" begin
    test_4()
end

@testset "dissipators 5" begin
    test_5()
end

@testset "dissipators 6" begin
    test_6()
end

@testset "dissipators 7" begin
    test_7()
end
