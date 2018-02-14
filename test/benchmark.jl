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
prop_suite = suite["unitary"]["propagator"]["rabi flops"] = BenchmarkGroup()
state_suite = suite["unitary"]["pure state"]["rabi flops"] = BenchmarkGroup()

# TODO replace when available in QSimulator
function qubit_drive(q::QSimulator.QSystem, drive::Function)
    function add_drive_ham!(ham, idxs, t)
        pulse = 2π * drive(t)
        drive_ham = real(pulse) * X(q) + imag(pulse) * Y(q)
        QSimulator.expand_add!(ham, drive_ham, idxs)
    end
end

qubit_freq = 5.0
nutation_freq = 0.02
for n = 2:4
    q0 = FixedDuffingTransmon("q0", qubit_freq, -0.2, n)
    cqs = CompositeQSystem([q0]);
    add_hamiltonian!(cqs, hamiltonian(q0), "q0")
    add_hamiltonian!(cqs, qubit_drive(q0, t -> nutation_freq*cos(2π*qubit_freq * t)), q0);
    ψ0 = zeros(Complex128, n)
    ψ0[1] = 1
    times = collect(linspace(0,100,101))
    prop_suite["$n level transmon"] = @benchmarkable unitary_propagator($cqs, $times);
    state_suite["$n level transmon"] = @benchmarkable unitary_state($cqs, $times, $ψ0);
end


# Lab frame parmetric interaction between two transmons with spectators

# Rotating frame parametric interaction between two trasmons with spectators
