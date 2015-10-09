using Base.Test
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

  L = zeros(Complex128, 4, 4)

  system += q1
  @test_approx_eq norm(L-QuantumInfo.hamiltonian(sz)',2) 0.0
  QSimulator.liouvillian_dual_add!(L, system, 0.0)

  system += Cooling("T1", 0.1, q1)
  @test_approx_eq norm(L-(QuantumInfo.hamiltonian(sz)+.1*QuantumInfo.dissipator(sm))',2) 0.0
  QSimulator.liouvillian_dual_add!(L, system, 0.0)
end

function test_5()
  q1 = Qubit("q1", 0.0)

  system = CompositeQSystem()

  system += q1

  system += Cooling("T1", 1.0, q1)

  Lp = liouvillian_propagator(system,100.,0.,1000.)

  @test_approx_eq_eps norm(Lp-[1 0 0 1; 0 0 0 0; 0 0 0 0; 0 0 0 0], 2) 0.0 1e-12
end

function test_6()
  q1 = Qubit("q1", 0.0)

  system = CompositeQSystem()

  system += q1

  system += Cooling("T1", 1.0, q1)

  state = [0. 0.; 0. 1.]

  ev_state = liouvillian_evolution(state,system,100.,0.,1000.)

  @test_approx_eq_eps norm(ev_state - [1 0; 0 0], 1) 0.0 1e-12
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
      ans = diagm([1-exp(-t),exp(-t)])

      @test_approx_eq_eps norm(ev_state - ans, 1) 0.0 1e-12
  end
end


test_1()
test_2()
test_3()
test_4()
test_5()
test_6()
test_7()
