using BenchmarkTools

using QSimulator

const suite = BenchmarkGroup()

suite["unitary"] = BenchmarkGroup()
suite["unitary"]["propagator"] = BenchmarkGroup()
suite["unitary"]["pure state"] = BenchmarkGroup()

# Free evolution of an n level transmon
prop_suite = suite["unitary"]["propagator"]["free evolution"] = BenchmarkGroup()
state_suite = suite["unitary"]["pure state"]["free evolution"] = BenchmarkGroup()
for n = 2:4
    q0 = FixedDuffingTransmon("q0", 5, -0.2, n)
    cqs = CompositeQSystem([q0])
    add_hamiltonian!(cqs, q0)
    times = collect(linspace(0,1.0,201))
    # superposition of all levels
    ψ0 = (1/sqrt(dim(cqs))) * ones(Complex128, dim(cqs))
    prop_suite["single transmon ($n levels)"] = @benchmarkable unitary_propagator($cqs, $times);
    state_suite["single transmon ($n levels)"] = @benchmarkable unitary_state($cqs, $times, $ψ0);
end

# Dipole interaction between a chain of fixed transmons
for n = 2:4
    qs = [FixedDuffingTransmon("q$ct", 4.0 + 0.1*ct, -(0.2 + 0.01*ct), 3) for ct = 0:(n-1)]
    cqs = CompositeQSystem(qs)
    for q = qs
        add_hamiltonian!(cqs, hamiltonian(q), label(q))
    end
    for (qa,qb) = zip(qs[1:end-1], qs[2:end])
        add_hamiltonian!(cqs, 0.005*Dipole(qa, qb), [qa, qb])
    end
    times = collect(linspace(0,1.0,201))
    # superposition of all levels
    ψ0 = (1/sqrt(dim(cqs))) * ones(Complex128, dim(cqs))
    prop_suite["dipole chain ($n transmons)"] = @benchmarkable unitary_propagator($cqs, $times);
    state_suite["dipole chain ($n transmons)"] = @benchmarkable unitary_state($cqs, $times, $ψ0);
end

# Lab frame Rabi oscillations of an n level transmon
prop_suite = suite["unitary"]["propagator"]["rabi flops"] = BenchmarkGroup()
state_suite = suite["unitary"]["pure state"]["rabi flops"] = BenchmarkGroup()

qubit_freq = 5.0
nutation_freq = 0.02
for n = 2:4
    q0 = FixedDuffingTransmon("q0", qubit_freq, -0.2, n)
    cqs = CompositeQSystem([q0]);
    add_hamiltonian!(cqs, q0)
    add_hamiltonian!(cqs, microwave_drive(q0, t -> nutation_freq*cos(2π*qubit_freq * t)), q0);
    ψ0 = zeros(Complex128, n)
    ψ0[1] = 1
    times = collect(linspace(0,100,101))
    prop_suite["$n level transmon"] = @benchmarkable unitary_propagator($cqs, $times);
    state_suite["$n level transmon"] = @benchmarkable unitary_state($cqs, $times, $ψ0);
end

# Lab frame parametric interaction between two transmons with spectators
prop_suite = suite["unitary"]["propagator"]["parametric 2Q gate"] = BenchmarkGroup()
state_suite = suite["unitary"]["pure state"]["parametric 2Q gate"] = BenchmarkGroup()

# helper function to add flux drive
for n = 2:4
    q0 = FixedDuffingTransmon("q0", 3.94015, -0.1807,  3)
    q1 = TunableDuffingTransmon("q1",  0.172, 16.4, 0.55, 3)

    # add fixed frequency spectators
    spectator_qs = [FixedDuffingTransmon("q$ct", 4.0 + 0.1*ct, -(0.2 + 0.01*ct), 3) for ct = 2:(n-1)]

    # should get an iSWAP interaction at ≈ 122 MHz with drive amplitude 0.323 Φ₀
    mod_freq = 122.1/1e3
    amp = 0.323
    all_qs = [q0, q1, spectator_qs...]
    cqs = CompositeQSystem(all_qs)
    add_hamiltonian!(cqs, q0)
    add_hamiltonian!(cqs, flux_drive(q1, t -> sin(2π*mod_freq*t)), q1)
    add_hamiltonian!(cqs, 0.006*Dipole(q0, q1), [q0,q1])

    # add hamiltoians for spectators coupled to tunable transmon
    for q = all_qs[3:end]
        add_hamiltonian!(cqs, q)
        add_hamiltonian!(cqs, 0.006*Dipole(q, q1), [q,q1])
    end

    times = collect(0:200)
    ψ0 = Complex128[0.0; 1.0; 0.0] ⊗ Complex128[1.0; 0.0; 0.0] # start in 10 state
    for ct = 1:n-2
        ψ0 = ψ0 ⊗ Complex128[1.0; 0.0; 0.0]
    end

    prop_suite["$n transmons"] = @benchmarkable unitary_propagator($cqs, $times)
    state_suite["$n transmons"] = @benchmarkable unitary_state($cqs, $times, $ψ0)
end


# Rotating frame parametric interaction between two transmons with spectators
