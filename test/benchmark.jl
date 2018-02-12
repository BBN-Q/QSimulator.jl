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
    add_hamiltonian!(cqs, hamiltonian(q0), "q0")
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

# Lab frame parmetric interaction between two transmons with spectators

# Rotating frame parametric interaction between two trasmons with spectators
