using Base.Test
using QSimulator
using QIP

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

  @test_throws hamiltonian(system,0.0)
end

function test_4()
  sz = [0 0; 0 1];
  sm = [0 1; 0 0];

  q1 = Qubit("q1", 1.0)

  system = CompositeQSystem()

  L = zeros(Complex128, 4, 4)

  system += q1
  QSimulator.liouvillian_add!(L, system, 0.0)
  @test_approx_eq norm(L-QIP.hamiltonian(sz)',2) 0.0

  system += Cooling("T1", 0.1, q1)
  QSimulator.liouvillian_add!(L, system, 0.0)
  @test_approx_eq norm(L-(QIP.hamiltonian(sz)+.1*QIP.dissipator(sm))',2) 0.0
end

function test_5()
  sz = [0 0; 0 1];
  sm = [0 1; 0 0];

  q1 = Qubit("q1", 0.0)

  system = CompositeQSystem()

  L = zeros(Complex128, 4, 4)

  system += q1

  system += Cooling("T1", 1.0, q1)
  QSimulator.liouvillian_add!(L, system, 0.0)

  @test_approx_eq_eps norm(liouvillian_propagator(system, 1000., 0., 1000.)-[1 0 0 1; 0 0 0 0; 0 0 0 0; 0 0 0 0], 2) 0.0 1e-12
end

test_1()
test_2()
test_3()
test_4()
test_5()