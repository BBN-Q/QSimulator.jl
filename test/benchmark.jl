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
    q0 = DuffingTransmon("q0", n, DuffingSpec(5, -0.2))
    cqs = CompositeQSystem([q0])
    add_hamiltonian!(cqs, q0)
    times = range(0, stop=1.0, length=201)
    # superposition of all levels
    ψ₀ = (1/sqrt(dimension(cqs))) * ones(ComplexF64, dimension(cqs))
    prop_suite["single transmon ($n levels)"] = @benchmarkable unitary_propagator($cqs, $times);
    state_suite["single transmon ($n levels)"] = @benchmarkable unitary_state($cqs, $times, $ψ₀);
end

# Dipole interaction between a chain of fixed transmons
for n = 2:4
    qs = [DuffingTransmon("q$ct", 3, DuffingSpec(4.0 + 0.1*ct, -(0.2 + 0.01*ct))) for ct = 0:(n-1)]
    cqs = CompositeQSystem(qs)
    for q = qs
        add_hamiltonian!(cqs, hamiltonian(q), q)
    end
    for (qa,qb) = zip(qs[1:end-1], qs[2:end])
        add_hamiltonian!(cqs, 0.005*X([qa, qb]), [qa, qb])
    end
    times = range(0, stop=1.0, length=201)
    # superposition of all levels
    ψ₀ = (1/sqrt(dimension(cqs))) * ones(ComplexF64, dimension(cqs))
    prop_suite["dipole chain ($n transmons)"] = @benchmarkable unitary_propagator($cqs, $times);
    state_suite["dipole chain ($n transmons)"] = @benchmarkable unitary_state($cqs, $times, $ψ₀);
end

# Lab frame Rabi oscillations of an n level transmon
prop_suite = suite["unitary"]["propagator"]["rabi flops"] = BenchmarkGroup()
state_suite = suite["unitary"]["pure state"]["rabi flops"] = BenchmarkGroup()

qubit_freq = 5.0
nutation_freq = 0.02
for n = 2:4
    q0 = DuffingTransmon("q0", n, DuffingSpec(qubit_freq, -0.2))
    cqs = CompositeQSystem([q0]);
    add_hamiltonian!(cqs, q0)
    add_hamiltonian!(cqs, dipole_drive(q0, t -> nutation_freq*cos(2π*qubit_freq * t)), q0);
    # start in ground state
    ψ₀ = zeros(ComplexF64, n)
    ψ₀[1] = 1
    times = range(0, stop=100, length=101)
    prop_suite["$n level transmon"] = @benchmarkable unitary_propagator($cqs, $times);
    state_suite["$n level transmon"] = @benchmarkable unitary_state($cqs, $times, $ψ₀);
end

# Lab frame parametric interaction between two transmons with spectators
prop_suite = suite["unitary"]["propagator"]["lab frame parametric 2Q gate"] = BenchmarkGroup()
state_suite = suite["unitary"]["pure state"]["lab frame parametric 2Q gate"] = BenchmarkGroup()

# helper function to add flux drive
for n = 2:3
    q0 = DuffingTransmon("q0", 3, DuffingSpec(3.94015, -0.1807))
    q1 = PerturbativeTransmon("q1",  3, TransmonSpec(0.172, 16.4, 0.55))

    # add fixed frequency spectators
    spectator_qs = [DuffingTransmon("q$ct", 3, DuffingSpec(4.0 + 0.1*ct, -(0.2 + 0.01*ct))) for ct = 2:(n-1)]

    # should get an iSWAP interaction at ≈ 122 MHz with drive amplitude 0.323 Φ₀
    mod_freq = 122.1/1e3
    amp = 0.323
    all_qs = [q0, q1, spectator_qs...]
    cqs = CompositeQSystem(all_qs)
    add_hamiltonian!(cqs, q0)
    add_hamiltonian!(cqs, parametric_drive(q1, t -> sin(2π*mod_freq*t)), q1)
    add_hamiltonian!(cqs, 0.006*X([q0, q1]), [q0, q1])

    # add hamiltoians for spectators coupled to tunable transmon
    for q = all_qs[3:end]
        add_hamiltonian!(cqs, q)
        add_hamiltonian!(cqs, 0.006*X([q, q1]), [q, q1])
    end

    times = 0:200
    ψ₀ = ComplexF64[0.0; 1.0; 0.0] ⊗ ComplexF64[1.0; 0.0; 0.0] # start in 10 state
    for ct = 1:n-2
        ψ₀ = ψ₀ ⊗ ComplexF64[1.0; 0.0; 0.0]
    end

    prop_suite["$n transmons"] = @benchmarkable unitary_propagator($cqs, $times)
    state_suite["$n transmons"] = @benchmarkable unitary_state($cqs, $times, $ψ₀)
end


# Rotating frame parametric interaction between two transmons with spectators
prop_suite = suite["unitary"]["propagator"]["rotating frame parametric 2Q gate"] = BenchmarkGroup()
state_suite = suite["unitary"]["pure state"]["rotating frame parametric 2Q gate"] = BenchmarkGroup()

for n = 2:4

    q0 = DuffingTransmon("q0", 3, DuffingSpec(3.94015, -0.1807))
    q1 = PerturbativeTransmon("q1",  3, TransmonSpec(0.172, 16.4, 0.55))

    # add fixed frequency spectators
    spectator_qs = [DuffingTransmon("q$ct", 3, DuffingSpec(4.0 + 0.1*ct, -(0.2 + 0.01*ct))) for ct = 2:(n-1)]

    # should get an CZ02 interaction at ≈ 115 MHz
    freq = 115.5/1e3
    amp = 0.245

    all_qs = [q0, q1, spectator_qs...]
    cqs = CompositeQSystem(all_qs)
    add_hamiltonian!(cqs, q0)
    # add rotating frame Hamiltonian shifts
    add_hamiltonian!(cqs, -spec(q0).frequency*number(q0), q0)
    q1_freq = hamiltonian(q1, 0.0)[2,2]
    add_hamiltonian!(cqs, -q1_freq*number(q1), q1)

    diff_freq = spec(q0).frequency - q1_freq
    # time dependent flip flop interaction
    add_hamiltonian!(cqs, t -> 0.006*.5*X_Y([q0, q1], [diff_freq*t, 0.0]), [q0,q1])
    add_hamiltonian!(cqs, parametric_drive(q1, t -> amp*sin(2π*freq*t)), q1)

    # add hamiltoians for spectators coupled to tunable transmon
    for q = all_qs[3:end]
        add_hamiltonian!(cqs, q)
        # add rotating frame Hamiltonian shifts
        add_hamiltonian!(cqs, -spec(q).frequency*number(q), q)
        diff_freq = q1_freq - spec(q).frequency
        add_hamiltonian!(cqs, t -> 0.006*.5*X_Y([q1, q], [diff_freq*t, 0.0]), [q1, q])
    end

    times = 0:200
    ψ₀ = ComplexF64[0.0; 1.0; 0.0] ⊗ ComplexF64[0.0; 1.0; 0.0] # start in 11 state
    for ct = 1:n-2
        ψ₀ = ψ₀ ⊗ ComplexF64[1.0; 0.0; 0.0]
    end

    prop_suite["$n transmons"] = @benchmarkable unitary_propagator($cqs, $times)
    state_suite["$n transmons"] = @benchmarkable unitary_state($cqs, $times, $ψ₀);
end
