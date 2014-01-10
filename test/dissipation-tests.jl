using Base.Test
using QSimulator

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

test_1()
test_2()