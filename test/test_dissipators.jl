using Base.Test,
      QSimulator

let
    sys = CompositeQSystem();	
    
    qA = Qubit("qA",8.0);
    qB = Qubit("qB",8.0);
    
    g=.1
    
    sys += qA;
    sys += qB;
    sys += FlipFlop(qA,qB,g);

    q0 = [1 0; 0 0];
    q1 = [0 0; 0 1];
    
    obs = kron(q1,q0);
    
    Rabi_rate = 2g;
    max_time = 3/Rabi_rate;
    times = linspace(.0,max_time,1000);

    result = 
      pmap(times) do time
          local u
          local out
          out = zeros(3)
          u = unitary_propagator(sys, .01, .0, time);
          out[1] = real(trace(obs*u*kron(q0,q0)*u'));
          out[2] = real(trace(obs*u*kron(q0,q1)*u'));
          out[3] = real(trace(obs*u*kron(q1,q0)*u'));  
          out
      end;

    pred1 = zeros(length(times))
    pred2 = map(t->1/2*cos(2pi*2g*t+pi)+1/2,times);
    pred3 = map(t->1/2*cos(2pi*2g*t)+1/2,times);

    r1 = Float64[r[1] for r in result];
    r2 = Float64[r[2] for r in result];
    r3 = Float64[r[3] for r in result];

    @test_approx_eq_eps(maximum(abs(r1-pred1)),0.0,1e-2)
    @test_approx_eq_eps(maximum(abs(r2-pred2)),0.0,1e-2)
    @test_approx_eq_eps(maximum(abs(r3-pred3)),0.0,1e-2)

end

let
    sys3 = CompositeQSystem();
    
    qA = Qubit("qA",8.0);

    rrate = .1
    relaxA = Cooling("T1",.1,qA)
    
    sys3 += qA
    sys3 += relaxA

    max_time = 3/rrate
    times = linspace(.0,max_time,1000);

    rhos = [expm(time*generator(sys3))*vec([0 0; 0 1]) for time in times];

    p = Float64[ (vec([0 0; 0 1])'*rho)[1,1] for rho in rhos ];
    pred = map(t->exp(-2pi*t*rate(relaxA)),times);

    @test_approx_eq_eps(maximum(abs(p-pred)),0.0,1e-2)
end
