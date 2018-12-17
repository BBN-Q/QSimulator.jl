using Test
using QSimulator

function test_1()

  for i=1:100
      r1 = Resonator("r1", 1.0, 10)

      system = CompositeQSystem()

      system += r1

      state = randn(10) + 1im*randn(10)
      state = state/norm(state,2)
      t = 1.0+10*abs(randn())
      #t = abs(randn())

      ans = state .* exp.(-1im*2*pi*t*Float64[k for k in 0:9])
      ev_state = unitary_evolution(state,system,.1,0.,t)

      @test isapprox(dot(ans,ev_state), 1.0, atol = 1e-12)
      #@test isapprox(norm(state, 2) 1.0 1e-12)
  end
end

@testset "unitary" begin
    test_1()
end
